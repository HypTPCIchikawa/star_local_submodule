#!/usr/bin/env python3
"""Validate QA YAML and emit a ROOT5-readable TSV manifest without overwriting.

Run from the repository root:
  python3 common/macro/prepareLambdaKfImp5Scan.py config/qa/lambda_kf_imp5_cut_scan_20260913.yaml /tmp/scan.tsv
Then (inside the established STAR environment, compiled with ACLiC):
  root4star -b -q 'common/macro/scanLambdaKfImp5Cuts.C+("/tmp/scan.tsv","new-output-directory")'

Input ROOT files are hashed read-only. This utility does not submit or run events.
"""
import argparse
import hashlib
import math
from pathlib import Path
import sys

import yaml


class UniqueKeyLoader(yaml.SafeLoader):
    pass


def unique_mapping(loader, node, deep=False):
    result = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if key in result:
            raise ValueError("duplicate YAML key: " + str(key))
        result[key] = loader.construct_object(value_node, deep=deep)
    return result


UniqueKeyLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, unique_mapping)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def finite(value):
    return not isinstance(value, bool) and isinstance(value, (int, float)) and math.isfinite(value)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[2]
    data = yaml.load(args.config.read_text(), Loader=UniqueKeyLoader)
    required = {"schemaVersion", "helixRoot", "kfRoot", "inputEvents", "massMin", "massMax",
                "massBinWidth", "window", "narrowWindow", "leftSideband", "rightSideband",
                "nearLeftSideband", "nearRightSideband", "farLeftSideband", "farRightSideband",
                "maxDistanceToPv", "maxDaughterDistance", "minCosPointing",
                "includeThreeDimensionalGrid", "representativeRule"}
    if not isinstance(data, dict) or set(data) != required:
        raise ValueError("missing or unknown QA YAML keys")
    if data["schemaVersion"] != 1:
        raise ValueError("unsupported schemaVersion")
    if not isinstance(data["inputEvents"], int) or isinstance(data["inputEvents"], bool) or data["inputEvents"] <= 0:
        raise ValueError("inputEvents must be a positive integer")
    for key in ("massMin", "massMax", "massBinWidth"):
        if not finite(data[key]):
            raise ValueError(key + " must be finite numeric")
    low, high, width = (data[k] for k in ("massMin", "massMax", "massBinWidth"))
    if width <= 0 or high <= low or not math.isclose((high-low)/width, round((high-low)/width), abs_tol=1.e-8):
        raise ValueError("invalid common histogram range/bin width")
    for key in ("window", "narrowWindow", "leftSideband", "rightSideband",
                "nearLeftSideband", "nearRightSideband", "farLeftSideband", "farRightSideband"):
        bounds = data[key]
        if not isinstance(bounds, list) or len(bounds) != 2 or not all(finite(x) for x in bounds):
            raise ValueError(key + " must have two finite bounds")
        if not low <= bounds[0] < bounds[1] <= high:
            raise ValueError(key + " lies outside the mass domain")
        if any(not math.isclose((x-low)/width, round((x-low)/width), abs_tol=1.e-8) for x in bounds):
            raise ValueError(key + " must align with common histogram bin edges")
    if not (data["leftSideband"][1] <= data["window"][0] <= data["narrowWindow"][0]
            < data["narrowWindow"][1] <= data["window"][1] <= data["rightSideband"][0]):
        raise ValueError("sidebands overlap the signal window, or narrow window is not nested")
    for prefix in ("near", "far"):
        if not (data[prefix + "LeftSideband"][1] <= data["window"][0]
                and data["window"][1] <= data[prefix + "RightSideband"][0]):
            raise ValueError(prefix + " sidebands overlap the signal window")
    for key in ("maxDistanceToPv", "maxDaughterDistance", "minCosPointing"):
        values = data[key]
        if not isinstance(values, list) or not values or not all(finite(x) for x in values):
            raise ValueError(key + " must be a nonempty numeric list")
        if len(set(values)) != len(values):
            raise ValueError(key + " has duplicate points")
        if key == "minCosPointing" and any(not -1 <= x <= 1 for x in values):
            raise ValueError("pointing cosine outside [-1,1]")
        if key == "maxDaughterDistance" and any(x < 0 for x in values):
            raise ValueError("negative daughter-distance threshold is not supported")
        if key == "maxDistanceToPv" and sum(x < 0 for x in values) > 1:
            raise ValueError("at most one disabled PV-distance point is permitted")
    if type(data["includeThreeDimensionalGrid"]) is not bool:
        raise ValueError("includeThreeDimensionalGrid must be boolean")
    if data["representativeRule"] != "nearest_window_count_then_fewest_changed_parameters":
        raise ValueError("unsupported representative rule")
    for key in ("helixRoot", "kfRoot"):
        source = Path(data[key])
        if not source.is_absolute():
            source = repo / source
        source = source.resolve(strict=True)
        data[key] = str(source)
        data[key + "Sha256"] = sha256(source)
    data["configPath"] = str(args.config.resolve(strict=True))
    data["configSha256"] = sha256(args.config)
    data["generatorSha256"] = sha256(Path(__file__))
    macro = Path(__file__).with_name("scanLambdaKfImp5Cuts.C")
    data["macroSha256"] = sha256(macro)
    rows = []
    for key in sorted(data):
        value = data[key]
        if isinstance(value, list):
            value = " ".join(format(x, ".17g") for x in value)
        elif isinstance(value, bool):
            value = "1" if value else "0"
        elif finite(value):
            value = format(value, ".17g")
        else:
            value = str(value)
        if any(c in value for c in "\t\r\n"):
            raise ValueError("manifest values may not contain tabs or newlines")
        rows.append(key + "\t" + value + "\n")
    # Exclusive creation: never replace previous manifests or analysis artifacts.
    with args.manifest.open("x") as stream:
        stream.writelines(rows)
    print("Created validated manifest:", args.manifest)
    print("Input hashes recorded; inputs were opened read-only.")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, TypeError, yaml.YAMLError) as error:
        print("ERROR:", error, file=sys.stderr)
        sys.exit(1)
