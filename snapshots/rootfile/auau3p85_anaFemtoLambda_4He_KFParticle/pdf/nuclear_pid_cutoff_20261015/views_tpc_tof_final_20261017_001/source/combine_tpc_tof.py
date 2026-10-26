#!/usr/bin/env python3
"""Combine already aggregated purity tables; never average purities or refit."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys

import yaml

SPECIES = ("d", "t", "he3", "he4")
VALUE_FIELDS = ("species", "rLo", "rHi", "purityNominal", "purityLow", "purityHigh", "rawSelected", "quality")


def sha(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_table(path):
    with Path(path).open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not set(VALUE_FIELDS).issubset(reader.fieldnames or []):
            raise ValueError("Missing required purity columns")
        return list(reader)


def quality(value):
    if value in ("1", "1.0", "True", "true"):
        return True
    if value in ("0", "0.0", "False", "false"):
        return False
    raise ValueError("Quality must be explicit")


def nullable_number(value):
    if value in ("", "None", "null", "nan"):
        return math.nan
    return float(value)


def checked_rows(rows, regime, minimum, maximum, targets):
    seen = set()
    for row in rows:
        species = row["species"]
        if species not in targets:
            raise ValueError("Unknown target")
        lo, hi = float(row["rLo"]), float(row["rHi"])
        key = species, round(lo, 9), round(hi, 9)
        if key in seen or not (hi > lo and lo >= minimum - 1e-8 and hi <= maximum + 1e-8):
            raise ValueError("Duplicate or out-of-range interval")
        seen.add(key)
        if regime == "TOF" and hi > float(targets[species]["displayMax"]) + 1e-8:
            raise ValueError("Point exceeds target domain")
        supported = quality(row["quality"])
        values = [nullable_number(row[k]) for k in ("purityNominal", "purityLow", "purityHigh")]
        if supported and not all(math.isfinite(x) for x in values):
            raise ValueError("Supported interval lacks finite purity")
        if any(math.isfinite(x) and not 0 <= x <= 1 for x in values):
            raise ValueError("Unphysical purity")
        if math.isfinite(values[1]) and math.isfinite(values[2]) and values[1] > values[2]:
            raise ValueError("Inverted model envelope")
        raw = float(row["rawSelected"])
        if not math.isfinite(raw) or raw < 0:
            raise ValueError("Invalid raw selected count")
        out = {k: row[k] for k in VALUE_FIELDS}
        out["quality"] = "1" if supported else "0"
        out["regime"] = regime
        yield out
    for species in targets:
        if not any(x[0] == species for x in seen):
            raise ValueError("Missing target coverage")


def check_coverage(rows, regime, minimum, maximum, targets, step):
    if not step > 0:
        raise ValueError("Explicit aggregate width must be positive")
    for species in targets:
        endpoint = maximum if regime == "TPC" else float(targets[species]["displayMax"])
        count = round((endpoint - minimum) / step)
        if abs(count * step - (endpoint - minimum)) > 1e-8:
            raise ValueError("Domain is not aligned to aggregate width")
        subset = sorted((float(row["rLo"]), float(row["rHi"])) for row in rows if row["species"] == species)
        if len(subset) != count:
            raise ValueError("Missing or extra aggregate intervals")
        for index, (lo, hi) in enumerate(subset):
            if abs(lo - (minimum + index * step)) > 1e-8 or abs(hi - lo - step) > 1e-8:
                raise ValueError("Gap or inconsistent aggregate interval")


def write_table(path, rows, fields):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(command, log, cwd, timeout):
    with Path(log).open("w") as stream:
        subprocess.run(command, cwd=cwd, stdout=stream, stderr=subprocess.STDOUT, check=True, timeout=timeout)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--settings", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    project = Path(__file__).resolve().parents[2]
    settings_path = Path(args.settings).resolve()
    cfg = yaml.safe_load(settings_path.read_text())
    if cfg["schemaVersion"] != 1 or set(cfg["targets"]) != set(SPECIES):
        raise ValueError("Unsupported schema or target set")
    paths = {key: (project / cfg[key]).resolve() for key in (
        "runtimeMainconf", "lowTable", "lowRanges", "selectionRanges", "lowPdf", "highTable", "highCompleted", "inputAudit", "plotSource", "lowPlotSource", "runtimeWrapper")}
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    complete = json.loads(paths["highCompleted"].read_text())
    if complete["executionComplete"] is not True or complete["protectedInputsVerified"] is not True:
        raise ValueError("High-R study not complete or inputs not verified")
    product = complete["products"][paths["highTable"].name]
    expected_sha = product["sha256"] if isinstance(product, dict) else product
    if sha(paths["highTable"]) != expected_sha:
        raise ValueError("High-R table differs from completion manifest")
    rmin, split, rmax = [float(cfg[k]) for k in ("rigidityMin", "rigiditySplit", "rigidityMax")]
    threshold = float(cfg["purityThreshold"])
    if not (rmin < split < rmax and 0 < threshold < 1):
        raise ValueError("Invalid explicit plotting range/threshold")
    low = list(checked_rows(read_table(paths["lowTable"]), "TPC", rmin, split, cfg["targets"]))
    high = list(checked_rows(read_table(paths["highTable"]), "TOF", split, rmax, cfg["targets"]))
    check_coverage(low, "TPC", rmin, split, cfg["targets"], float(cfg["aggregateWidth"]))
    check_coverage(high, "TOF", split, rmax, cfg["targets"], float(cfg["aggregateWidth"]))
    if int(cfg["lowConditionCount"]) != 20:
        raise ValueError("Renderer low-R condition label differs from frozen study")
    audit = json.loads(paths["inputAudit"].read_text())
    windows = {row["cacheKey"]: row for row in audit["cache"]["tagWindows"]}
    metadata = []
    for species in SPECIES:
        target = cfg["targets"][species]
        window = windows[target["sample"]]
        if window["strictBounds"] is not True:
            raise ValueError("Unexpected window boundary")
        condition = "{:.6g} < m^{{2}} < {:.6g}; R >= {:.3g}".format(window["savedMass2Min"], window["savedMass2Max"], split)
        if species == "he4" and target["sample"] != cfg["targets"]["d"]["sample"]:
            raise ValueError("He4 study must use the user-selected D window")
        metadata.append(dict(species=species, hasHigh=1, exactStudyWindow=1, displayMax=target["displayMax"],
                             conditionLabel=condition, detailLabel=target["detail"], missingLabel="High-R purity not evaluated",
                             curveLabel="TPC+TOF ({} fit conditions)".format(cfg["highConditionCount"])))
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=False)
    source = output / "source"
    source.mkdir()
    protected = {str(path): sha(path) for path in list(paths.values()) + [settings_path]}
    for key, path in paths.items():
        if key in ("plotSource", "lowPlotSource", "runtimeWrapper", "runtimeMainconf"):
            shutil.copy2(path, source / path.name)
    shutil.copy2(settings_path, source / settings_path.name)
    shutil.copy2(Path(__file__), source / Path(__file__).name)
    shutil.copy2(paths["highCompleted"], source / "high_completed.json")
    table, meta = output / "combined.tsv", output / "metadata.tsv"
    write_table(table, low + high, VALUE_FIELDS + ("regime",))
    write_table(meta, metadata, tuple(metadata[0]))
    shutil.copy2(paths["lowRanges"], output / "low_adopted_ranges.tsv")
    ranges = output / "production_adopted_ranges.tsv"
    shutil.copy2(paths["selectionRanges"], ranges)
    overview = output / "purity_tpc_tof_overview.pdf"
    call_args = [str(table), str(meta), str(ranges), str(overview)]
    call = str(source / paths["plotSource"].name) + "+(" + ",".join(json.dumps(x) for x in call_args)
    call += "," + ",".join(format(x, ".17g") for x in (rmin, split, rmax, threshold)) + ")"
    command = ["bash", str(paths["runtimeWrapper"]), str(paths["runtimeMainconf"]), "root4star", "-l", "-b", "-q", call]
    (output / "command.json").write_text(json.dumps(command, indent=2) + "\n")
    timeout = float(cfg["executionTimeoutSeconds"])
    run(command, output / "render.log", project, timeout)
    combined = output / "purity_tpc_tof_combined.pdf"
    run(["gs", "-dBATCH", "-dNOPAUSE", "-sDEVICE=pdfwrite", "-sOutputFile=" + str(combined), str(overview), str(paths["lowPdf"])], output / "merge.log", project, timeout)
    pages = {}
    for pdf, expected in ((overview, 3), (combined, 6)):
        log = output / (pdf.name + ".readback.log")
        run(["gs", "-dBATCH", "-dNOPAUSE", "-sDEVICE=nullpage", str(pdf)], log, project, timeout)
        count = len(re.findall(r"^Page [0-9]+\s*$", log.read_text(), re.MULTILINE))
        if count != expected:
            raise ValueError("Unexpected PDF page count")
        pages[pdf.name] = count
    if any(sha(Path(path)) != value for path, value in protected.items()):
        raise ValueError("Protected input changed during rendering")
    report = dict(schemaVersion=1, executionComplete=True, presentationOnly=True, refitted=False,
                  lowRows=len(low), highRows=len(high), pageCounts=pages, protectedInputs=protected,
                  note="He4 uses D mass2 window by user choice; old 3-4 He4 window is not evaluated",
                  products={str(p.relative_to(output)): sha(p) for p in output.rglob("*") if p.is_file()})
    (output / "completed.json").write_text(json.dumps(report, indent=2) + "\n")
    print("TPC_TOF_COMBINED_COMPLETE", output)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print("TPC_TOF_COMBINED_FAILED:", error, file=sys.stderr)
        raise
