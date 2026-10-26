#!/usr/bin/env python3
"""Reproducible, train-only, bounded post-KF Lambda purity scan.

Inputs are read-only CSV exports verified by exportLambdaPurityCandidates.C.
Run with host Python (STAR PYTHONPATH can select incompatible NumPy binaries):
  env -u PYTHONPATH /usr/bin/python3 common/macro/studyLambdaKfPurity.py \
    --config config/qa/lambda_kf_purity_study_20260927.yaml \
    --datasets datasets.tsv --output NEW_OUTPUT_DIRECTORY

This is candidate-output QA, not a replacement PicoDst/StChain event loop.
No mass constraint or mass-window optimization is performed. An optimum over
this finite, same-data exploratory grid is not a physics-optimal selection.
"""
import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
import sys

import numpy as np
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

# This mapping documents the existing Maker/interface API, not cut values.
CUT_FIELDS = {
    "nSigmaProton": ("protonPidPull", "abs_le"),
    "nSigmaPion": ("pionPidPull", "abs_le"),
    "imp5MinDCAProton": ("protonDcaToPv", "ge"),
    "imp5MinDCAPion": ("pionDcaToPv", "ge"),
    "maxDaughterDistance": ("daughterDistance", "le"),
    "maxDistanceToPv": ("distanceToPv", "le"),
    "minCosPointing": ("cosPointing", "ge"),
    "maxChi2Ndf": ("chi2Ndf", "le"),
    "maxTopoChi2Ndf": ("topoChi2Ndf", "le"),
    "minDecayLength": ("decayLength", "ge"),
    "minDecayLengthSignificance": ("decayLengthSignificance", "ge"),
    "minVertexLineSignificance": ("vertexLineLengthSignificance", "ge"),
}
UPSTREAM_DIRECTIONS = {
    "interfaceChiPrimaryCut": "ge", "finderChiPrimary2D": "ge",
    "finderLCut": "ge", "finderLdL2D": "ge",
    "finderChi2Ndf2D": "le", "topoChi2NdfCut": "le",
}
METRIC_FIELDS = [
    "Nwin", "Nleft", "Nright", "alpha", "B", "BError", "S", "SError",
    "purity", "purityError", "SoverB", "SoverBError", "SoverSqrtN",
    "signalRetention", "NwinRetention", "baselineS", "baselineNwin",
]
ROW_FIELDS = ["dataset_id", "point_id", "split", "sideband"] + METRIC_FIELDS + [
    "sourceChangedParameters", "addedCutCount", "changedParameters", "cuts_json", "score",
    "selectionEligible", "eligibilityReason",
]


def require(test, message):
    if not test:
        raise ValueError(message)


def finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def sha256(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def load_config(path):
    with Path(path).open() as source:
        c = yaml.load(source, Loader=UniqueKeyLoader)
    required = {
        "schemaVersion", "inputEvents", "reconstructedEvents", "eventIdsFile", "massMin", "massMax",
        "massBinWidth", "window", "sidebands", "displayMin", "displayMax", "title",
        "splitSeed", "trainFraction", "selectionScore", "signalRetentionFloors",
        "minWindowCountForSelection", "minSidebandCountForSelection", "maxPointsPerDataset",
        "bootstrapReplicates", "bootstrapSeed", "bootstrapMinValidFraction", "bootstrapQuantiles",
        "allowedUpstreamDifferences", "cuts", "includeAllPairs", "tripleFamilies",
        "fitRanges", "fitSigmaMin", "fitSigmaMax", "fitMeanMin", "fitMeanMax", "fitMinEntries",
    }
    require(isinstance(c, dict) and set(c) == required, "missing/unknown QA configuration keys")
    require(c["schemaVersion"] == 1, "unsupported schemaVersion")
    for key in ("inputEvents", "reconstructedEvents", "minWindowCountForSelection",
                "minSidebandCountForSelection", "maxPointsPerDataset", "bootstrapReplicates",
                "fitMinEntries"):
        require(type(c[key]) is int and c[key] > 0, key + " must be a positive integer")
    require(c["reconstructedEvents"] <= c["inputEvents"], "invalid event counts")
    require(isinstance(c["eventIdsFile"], str) and c["eventIdsFile"], "missing eventIdsFile")
    for key in ("splitSeed", "bootstrapSeed"):
        require(type(c[key]) is int and 0 <= c[key] < 2**32, key + " must fit uint32")
    require(c["bootstrapReplicates"] >= 200, "at least 200 bootstrap replicates required")
    require(finite(c["trainFraction"]) and 0 < c["trainFraction"] < 1, "invalid trainFraction")
    require(c["selectionScore"] in ("purity_minus_one_statistical_error", "purity"), "unknown ranking rule")
    require(finite(c["bootstrapMinValidFraction"]) and 0 < c["bootstrapMinValidFraction"] <= 1,
            "invalid bootstrapMinValidFraction")
    require(isinstance(c["title"], str) and c["title"] and "\n" not in c["title"] and "\t" not in c["title"],
            "invalid plot title")
    for key in ("massMin", "massMax", "massBinWidth", "displayMin", "displayMax",
                "fitSigmaMin", "fitSigmaMax", "fitMeanMin", "fitMeanMax"):
        require(finite(c[key]), key + " must be finite numeric")
    lo, hi, width = (c[k] for k in ("massMin", "massMax", "massBinWidth"))
    require(lo < hi and width > 0 and math.isclose((hi - lo) / width, round((hi - lo) / width), abs_tol=1e-8),
            "invalid common histogram axis")
    require(lo <= c["displayMin"] < c["displayMax"] <= hi, "display outside mass axis")
    require(0 < c["fitSigmaMin"] < c["fitSigmaMax"], "invalid fit sigma range")
    require(lo < c["fitMeanMin"] < c["fitMeanMax"] < hi, "invalid fit mean range")
    def bounds(value, name):
        require(isinstance(value, list) and len(value) == 2 and all(finite(x) for x in value), name + " needs two bounds")
        require(lo <= value[0] < value[1] <= hi, name + " outside axis")
        for x in value:
            require(math.isclose((x-lo)/width, round((x-lo)/width), abs_tol=1e-8), name + " not bin-aligned")
    bounds(c["window"], "window")
    require(isinstance(c["sidebands"], dict) and set(c["sidebands"]) == {"default", "near", "far"}, "need default/near/far sidebands")
    for name, sb in c["sidebands"].items():
        require(isinstance(sb, dict) and set(sb) == {"left", "right"}, "invalid sideband fields")
        bounds(sb["left"], name + " left"); bounds(sb["right"], name + " right")
        require(sb["left"][1] <= c["window"][0] and c["window"][1] <= sb["right"][0], "overlapping sidebands")
    require(isinstance(c["fitRanges"], list) and c["fitRanges"], "empty fit ranges")
    for item in c["fitRanges"]:
        bounds(item, "fitRange")
        require(item[0] < c["window"][0] and item[1] > c["window"][1], "fit range must span mass window")
    floors = c["signalRetentionFloors"]
    require(isinstance(floors, list) and floors and all(finite(x) and 0 < x <= 1 for x in floors), "invalid retention floors")
    require(len(set(floors)) == len(floors), "duplicate retention floors")
    quantiles = c["bootstrapQuantiles"]
    require(isinstance(quantiles, list) and len(quantiles) == 2 and all(finite(x) for x in quantiles)
            and 0 < quantiles[0] < quantiles[1] < 1, "invalid bootstrap quantiles")
    require(isinstance(c["allowedUpstreamDifferences"], list) and
            set(c["allowedUpstreamDifferences"]).issubset(UPSTREAM_DIRECTIONS), "unsupported upstream differences")
    require(isinstance(c["cuts"], dict) and c["cuts"] and set(c["cuts"]).issubset(CUT_FIELDS), "unknown/empty scan axes")
    for name, spec in c["cuts"].items():
        require(isinstance(spec, dict) and set(spec) == {"field", "operator", "values"}, "invalid axis " + name)
        require((spec["field"], spec["operator"]) == CUT_FIELDS[name], "wrong field/operator for " + name)
        values = spec["values"]
        require(isinstance(values, list) and values and all(finite(v) for v in values), "invalid thresholds for " + name)
        require(len(set(values)) == len(values), "duplicate threshold for " + name)
        require(all(v > 0 for v in values), "nonpositive study threshold for " + name)
        if name == "minCosPointing":
            require(all(v <= 1 for v in values), "cosine threshold >1")
    require(type(c["includeAllPairs"]) is bool, "includeAllPairs must be boolean")
    require(isinstance(c["tripleFamilies"], list), "invalid triple families")
    for family in c["tripleFamilies"]:
        require(isinstance(family, list) and len(family) == 3 and len(set(family)) == 3
                and set(family).issubset(c["cuts"]), "invalid triple family")
    return c


def enumerate_points(c):
    points = [()]
    seen = {()}
    axes = list(c["cuts"])
    families = [(key,) for key in axes]
    if c["includeAllPairs"]:
        families.extend(itertools.combinations(axes, 2))
    families.extend(c["tripleFamilies"])
    for family in families:
        for levels in itertools.product(*(c["cuts"][key]["values"] for key in family)):
            point = tuple(sorted(zip(family, map(float, levels))))
            if point not in seen:
                seen.add(point); points.append(point)
    require(len(points) <= c["maxPointsPerDataset"], "grid exceeds configured cap")
    return points


def event_train(run_id, event_id, seed, fraction):
    # Canonical decimal integer encoding; independent of Python randomized hash.
    key = "{}:{}:{}".format(seed, int(run_id), int(event_id)).encode("ascii")
    value = int.from_bytes(hashlib.sha256(key).digest()[:8], "big")
    return value < int(fraction * (1 << 64))


def read_numeric_csv(path):
    with Path(path).open(newline="") as source:
        reader = csv.reader(source)
        header = next(reader, None)
        require(header and len(header) == len(set(header)), "missing/duplicate CSV columns: " + str(path))
        rows = []
        for line_no, row in enumerate(reader, 2):
            require(len(row) == len(header), "wrong column count in {}:{}".format(path, line_no))
            rows.append([float(value) for value in row])
    require(rows, "empty CSV: " + str(path))
    values = np.asarray(rows, dtype=np.float64)
    require(np.isfinite(values).all(), "nonfinite CSV value: " + str(path))
    return {name: values[:, index] for index, name in enumerate(header)}


def read_metadata(path):
    text = Path(path).read_text()
    chunks = text.split("\nKFParticleEffectiveConfiguration\n")
    require(len(chunks) == 2, "missing/duplicate KF effective metadata")
    tail = chunks[1].split("\nKFEventSelectionConfiguration\n")
    require(len(tail) == 2, "missing/duplicate event metadata")
    sections = [yaml.load(x, Loader=UniqueKeyLoader) for x in (chunks[0], tail[0], tail[1])]
    require(all(isinstance(x, dict) for x in sections), "malformed export metadata")
    return sections


def mass_bins(mass, c):
    low, high = c["massMin"], c["massMax"]
    bins = int(round((high - low) / c["massBinWidth"]))
    # ROOT fixed-axis convention, explicitly retaining under/overflow.
    index = np.floor(bins * ((mass - low) / (high - low))).astype(np.int64) + 1
    index[mass < low] = 0
    index[mass >= high] = bins + 1
    return np.clip(index, 0, bins + 1)


def replay_selected(data, cuts):
    required = set(CUT_FIELDS) | {"imp5MaxPathLength", "minMass", "maxMass", "maxMassError"}
    require(required.issubset(cuts), "missing replay cut metadata")
    mask = (data["imp5PathValid"] == 1)
    mask &= np.abs(data["protonHelixPathLength"]) <= cuts["imp5MaxPathLength"]
    mask &= np.abs(data["pionHelixPathLength"]) <= cuts["imp5MaxPathLength"]
    mask &= (data["mass"] >= cuts["minMass"]) & (data["mass"] <= cuts["maxMass"])
    if cuts["maxMassError"] >= 0:
        mask &= data["massError"] <= cuts["maxMassError"]
    for key, (field, op) in CUT_FIELDS.items():
        # Original input PID/DCA already enforced; do not reapply a rounded
        # Float_t PID here. Validate those input gates separately below.
        if key in ("nSigmaProton", "nSigmaPion", "imp5MinDCAProton", "imp5MinDCAPion"):
            continue
        threshold = cuts[key]
        active = threshold > -1 if key == "minCosPointing" else threshold >= 0
        if active:
            mask &= data[field] <= threshold if op == "le" else data[field] >= threshold
    return mask


def load_dataset(dataset_id, folder, c, reference=None):
    paths = {key: folder / filename for key, filename in
             (("candidates", "candidates.csv"), ("histogram", "selected_histogram.csv"), ("metadata", "metadata.txt"))}
    require(all(path.is_file() for path in paths.values()), "missing export files: " + str(folder))
    data = read_numeric_csv(paths["candidates"])
    meta, cuts, event_meta = read_metadata(paths["metadata"])
    required_fields = {x[0] for x in CUT_FIELDS.values()} | {
        "runId", "eventId", "selected", "mass", "massError", "pdg", "imp5PathValid",
        "protonHelixPathLength", "pionHelixPathLength", "protonId", "pionId", "protonHasTof", "pionHasTof",
    }
    require(required_fields.issubset(data), "missing required candidate branches: " + dataset_id)
    n = len(data["mass"])
    require(meta.get("inputEvents") == c["inputEvents"] and meta.get("reconstructedEvents") == c["reconstructedEvents"],
            "wrong input/reconstructed counts: " + dataset_id)
    require(meta.get("rawEntries") == n, "wrong rawEntries: " + dataset_id)
    for name in ("selected", "imp5PathValid", "protonHasTof", "pionHasTof"):
        require(np.isin(data[name], (0, 1)).all(), "nonboolean " + name)
    for name in ("runId", "eventId", "protonId", "pionId"):
        require((data[name] == np.floor(data[name])).all() and (np.abs(data[name]) < 2**53).all(), "invalid integer " + name)
    require((data["pdg"] == 3122).all(), "study requires Lambda only")
    require(cuts.get("selectionProfile") == "lambda_imp5" and cuts.get("pidProfile") == "pico_nsigma"
            and cuts.get("backend") == "finder_topo" and cuts.get("schemaVersion") == 2,
            "not a full-KF TPC Imp5 source")
    for name in ("useTof", "strictTofPid", "cleanKaonsWithTof", "useHftTracksOnly", "reconstructAntiLambda"):
        require(cuts.get(name) in (False, 0), "unexpected source profile switch: " + name)
    require(not np.any(data["protonHasTof"]) and not np.any(data["pionHasTof"]), "unexpected Imp5 TOF values")
    require(cuts.get("applyLambdaGeometryCuts") in (True, 1), "purity study requires original Finder geometry enabled")
    for species in ("Proton", "Pion"):
        field = species.lower() + "PidPull"
        require((np.abs(data[field]) <= cuts["nSigma" + species] + 1e-6).all(), "upstream PID violation")
        require((data[species.lower() + "DcaToPv"] >= cuts["imp5MinDCA" + species]).all(), "upstream DCA violation")
    selected = data["selected"] == 1
    replay = replay_selected(data, cuts)
    require(np.array_equal(selected, replay), "selected bits do not exactly replay: " + dataset_id)
    require(selected.any(), "no selected baseline candidates: " + dataset_id)
    native = read_numeric_csv(paths["histogram"])
    require({"bin", "low", "high", "count", "error"}.issubset(native), "invalid histogram CSV columns")
    nbins = int(round((c["massMax"] - c["massMin"]) / c["massBinWidth"]))
    require(np.array_equal(native["bin"], np.arange(nbins+2)), "histogram needs contiguous under/regular/overflow bins")
    expected_low = c["massMin"] + (np.arange(nbins+2)-1) * c["massBinWidth"]
    require(np.allclose(native["low"], expected_low, rtol=0, atol=1e-12)
            and np.allclose(native["high"], expected_low+c["massBinWidth"], rtol=0, atol=1e-12), "different source histogram axis")
    bins = mass_bins(data["mass"], c)
    counts = np.bincount(bins[selected], minlength=nbins+2)
    require(np.array_equal(counts, native["count"]), "source mass histogram does not exactly replay: " + dataset_id)
    require(np.allclose(np.sqrt(counts), native["error"], rtol=1e-12, atol=1e-12), "source is not an unweighted count histogram")
    for key, spec in c["cuts"].items():
        original = cuts[key]
        for value in spec["values"]:
            enabled = original > -1 if key == "minCosPointing" else original >= 0
            if enabled:
                require(value <= original if spec["operator"] in ("le", "abs_le") else value >= original,
                        "scan must only tighten existing source bound: " + key)
    changes = {}
    if reference is not None:
        require(set(cuts) == set(reference["cuts"]), "different metadata key set across datasets")
        for key, value in cuts.items():
            if value != reference["cuts"][key]:
                require(key in c["allowedUpstreamDifferences"], "unapproved effective difference: " + key)
                old = reference["cuts"][key]
                require(value >= old if UPSTREAM_DIRECTIONS[key] == "ge" else value <= old,
                        "upstream variant must tighten: " + key)
                changes[key] = {"baseline": old, "value": value}
        path_keys = {"KFEventSelection.mainconf", "KFEventSelection.analysisInfo", "KFEventSelection.eventConfig"}
        require({k:v for k,v in event_meta.items() if k not in path_keys} ==
                {k:v for k,v in reference["event_meta"].items() if k not in path_keys}, "different effective event policies")
    keys = [(int(run), int(event)) for run, event in zip(data["runId"], data["eventId"])]
    flags = {key: event_train(*key, c["splitSeed"], c["trainFraction"]) for key in set(keys)}
    train = np.array([flags[key] for key in keys], dtype=bool)
    require(np.any(selected & train) and np.any(selected & ~train), "empty train/validation subset")
    return {"id": dataset_id, "folder": str(folder), "paths": paths, "meta": meta, "cuts": cuts,
            "event_meta": event_meta, "data": data, "selected": selected, "bins": bins,
            "splits": {"full": np.ones(n, dtype=bool), "train": train, "validation": ~train},
            "event_keys": keys, "upstream_changes": changes}


def make_region_masks(data, c):
    masks = {}
    for name, bounds in [("signal", c["window"])] + [
            (sb + "_" + side, c["sidebands"][sb][side]) for sb in ("default", "near", "far") for side in ("left", "right")]:
        masks[name] = (data["mass"] >= bounds[0]) & (data["mass"] < bounds[1])
    return masks


def cut_masks(dataset, c):
    cache = {}
    for key, spec in c["cuts"].items():
        values = dataset["data"][spec["field"]]
        if spec["operator"] == "abs_le":
            values = np.abs(values)
        for threshold in spec["values"]:
            cache[(key, float(threshold))] = values >= threshold if spec["operator"] == "ge" else values <= threshold
    return cache


def point_mask(dataset, point):
    mask = dataset["selected"].copy()
    for gate in point:
        mask &= dataset["cut_masks"][gate]
    return mask


def measure(n, left, right, alpha):
    n, left, right = float(n), float(left), float(right)
    b = alpha * (left + right)
    vb = alpha**2 * (left + right)
    signal = n - b
    result = {"Nwin": n, "Nleft": left, "Nright": right, "alpha": alpha,
              "B": b, "BError": math.sqrt(vb), "S": signal, "SError": math.sqrt(n+vb),
              "purity": float("nan"), "purityError": float("nan"), "SoverB": float("nan"),
              "SoverBError": float("nan"), "SoverSqrtN": float("nan")}
    if n > 0:
        result.update(purity=signal/n, purityError=math.sqrt(vb/n**2 + b*b/n**3), SoverSqrtN=signal/math.sqrt(n))
    if b > 0:
        result.update(SoverB=signal/b, SoverBError=math.sqrt(n/b**2 + n*n*vb/b**4))
    return result


def metrics(dataset, mask, split, sideband, c):
    use = mask & dataset["splits"][split]
    reg = dataset["regions"]
    sb = c["sidebands"][sideband]
    alpha = (c["window"][1]-c["window"][0]) / (sb["left"][1]-sb["left"][0]+sb["right"][1]-sb["right"][0])
    return measure(np.count_nonzero(use & reg["signal"]), np.count_nonzero(use & reg[sideband+"_left"]),
                   np.count_nonzero(use & reg[sideband+"_right"]), alpha)


def add_reference(result, baseline):
    result = dict(result)
    result["baselineS"], result["baselineNwin"] = baseline["S"], baseline["Nwin"]
    result["signalRetention"] = result["S"]/baseline["S"] if baseline["S"] > 0 else float("nan")
    result["NwinRetention"] = result["Nwin"]/baseline["Nwin"] if baseline["Nwin"] > 0 else float("nan")
    return result


def score(result, c):
    return result["purity"] - (result["purityError"] if c["selectionScore"] == "purity_minus_one_statistical_error" else 0)


def row_for(dataset, point_id, point, split, sb, result, c):
    nsource = len(dataset["upstream_changes"])
    row = {"dataset_id": dataset["id"], "point_id": point_id, "split": split, "sideband": sb,
           "sourceChangedParameters": nsource, "addedCutCount": len(point),
           "changedParameters": len(point)+nsource, "cuts_json": json.dumps(dict(point), sort_keys=True), "score": score(result, c)}
    row.update(result)
    reasons = []
    if split != "train" or sb != "default":
        reasons.append("not_training_default")
    if result["Nwin"] < c["minWindowCountForSelection"]:
        reasons.append("low_window_count")
    if result["Nleft"] + result["Nright"] < c["minSidebandCountForSelection"]:
        reasons.append("low_sideband_count")
    if result["S"] <= 0:
        reasons.append("nonpositive_signal")
    if not math.isfinite(row["score"]):
        reasons.append("nonfinite_score")
    row["selectionEligible"] = int(not reasons)
    row["eligibilityReason"] = ";".join(reasons) if reasons else "eligible_before_retention_floor"
    return row


def write_csv(path, fields, rows):
    with Path(path).open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def json_clean(value):
    if isinstance(value, dict):
        return {key: json_clean(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_clean(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_json(path, value):
    Path(path).write_text(json.dumps(json_clean(value), indent=2, sort_keys=True, allow_nan=False)+"\n")


def choose_points(train_rows, c):
    choices = []
    for floor in c["signalRetentionFloors"]:
        eligible = [row for row in train_rows if row["signalRetention"] >= floor and row["S"] > 0
                    and row["Nwin"] >= c["minWindowCountForSelection"]
                    and row["Nleft"]+row["Nright"] >= c["minSidebandCountForSelection"]
                    and math.isfinite(row["score"])]
        require(eligible, "no training candidate meets retention floor " + str(floor))
        # Numeric score first, exact-score ties prefer fewer changes, then
        # higher retained signal and deterministic dataset/point identity.
        best = min(eligible, key=lambda row: (-row["score"], row["changedParameters"], -row["S"],
                                             row["dataset_id"], row["point_id"]))
        choices.append((floor, best))
    return choices


def write_bootstrap(output, choices, datasets, masks, c, ordered_events):
    by_id = {d["id"]: d for d in datasets}
    keys = ordered_events
    event_lookup = {key: index for index, key in enumerate(keys)}
    event_ids = np.asarray(keys, dtype=np.int64)
    train_flags = np.array([event_train(*key, c["splitSeed"], c["trainFraction"]) for key in keys], dtype=bool)
    labels = ["signal"] + [sb+"_"+side for sb in ("default", "near", "far") for side in ("left", "right")]
    values = np.zeros((len(choices), len(keys), len(labels)), dtype=np.int64)
    for ci, choice in enumerate(choices):
        dataset = by_id[choice["dataset_id"]]
        event_index = np.array([event_lookup[key] for key in dataset["event_keys"]], dtype=np.int64)
        for ri, name in enumerate(labels):
            use = masks[choice["choice_id"]] & dataset["regions"][name]
            values[ci, :, ri] = np.bincount(event_index[use], minlength=len(keys))
    np.savez_compressed(output / "chosen_event_counts.npz", event_ids=event_ids,
                        is_train=train_flags, choice_ids=np.array([x["choice_id"] for x in choices]),
                        region_labels=np.array(labels), counts=values)
    rng = np.random.RandomState(c["bootstrapSeed"])
    # All 10,000 input events, including zero-candidate events, are retained.
    # Shared weights preserve baseline/selection covariance. This is a Poisson
    # event bootstrap, not candidate bootstrap or fixed-N Efron resampling.
    fields = ["choice_id", "split", "sideband", "quantity", "estimate", "bootstrapError",
              "quantileLow", "quantileHigh", "validReplicates", "replicates", "status"]
    rows = []
    all_draws = {}
    for split, use in (("full", np.ones(len(keys), dtype=bool)), ("train", train_flags), ("validation", ~train_flags)):
        local = values[:, use, :]
        weights = rng.poisson(1, size=(c["bootstrapReplicates"], int(use.sum())))
        # Flatten choices/regions to use a single small matrix multiplication.
        sums = weights.dot(local.transpose(1, 0, 2).reshape(int(use.sum()), -1))
        sums = sums.reshape(c["bootstrapReplicates"], len(choices), len(labels))
        actual = local.sum(axis=1)
        for si, sb_name in enumerate(("default", "near", "far")):
            sb = c["sidebands"][sb_name]
            alpha = (c["window"][1]-c["window"][0])/(sb["left"][1]-sb["left"][0]+sb["right"][1]-sb["right"][0])
            ns = sums[:, :, 0].astype(float)
            bs = alpha*(sums[:, :, 1+si*2]+sums[:, :, 2+si*2])
            ss = ns-bs
            actual_n = actual[:, 0].astype(float)
            actual_s = actual_n-alpha*(actual[:, 1+si*2]+actual[:, 2+si*2])
            with np.errstate(divide="ignore", invalid="ignore"):
                purities = ss/ns
                retention = ss/ss[:, :1]
                gain = purities-purities[:, :1]
                retention[ss[:, 0] <= 0, :] = np.nan
                actual_purity = actual_s/actual_n
                actual_retention = actual_s/actual_s[0]
            for ci, choice in enumerate(choices):
                for quantity, array, estimate in (
                        ("purity", purities[:, ci], actual_purity[ci]),
                        ("purityGain", gain[:, ci], actual_purity[ci]-actual_purity[0]),
                        ("signalRetention", retention[:, ci], actual_retention[ci])):
                    good = array[np.isfinite(array)]
                    valid = len(good) >= math.ceil(c["bootstrapMinValidFraction"]*c["bootstrapReplicates"])
                    quant = np.quantile(good, c["bootstrapQuantiles"]) if valid else [float("nan")]*2
                    rows.append({"choice_id": choice["choice_id"], "split": split, "sideband": sb_name,
                                 "quantity": quantity, "estimate": float(estimate),
                                 "bootstrapError": float(np.std(good, ddof=1)) if valid else float("nan"),
                                 "quantileLow": float(quant[0]), "quantileHigh": float(quant[1]),
                                 "validReplicates": len(good), "replicates": c["bootstrapReplicates"],
                                 "status": "ok" if valid else "insufficient_finite_replicates"})
                    all_draws["{}_{}_{}_{}".format(choice["choice_id"], split, sb_name, quantity)] = array
    write_csv(output / "paired_event_bootstrap.csv", fields, rows)
    np.savez_compressed(output / "paired_event_bootstrap_draws.npz", **all_draws)


def run(config_path, datasets_path, output):
    c = load_config(config_path)
    points = enumerate_points(c)
    require(not output.exists(), "output directory already exists; refusing overwrite")
    event_path = Path(c["eventIdsFile"])
    if not event_path.is_absolute():
        event_path = Path(__file__).resolve().parents[2] / event_path
    ordered_events = []
    with event_path.open() as source:
        for index, line in enumerate(source):
            fields = line.rstrip("\n").split("\t")
            require(len(fields) == 3, "eventIdsFile requires index/runId/eventId, headerless TSV")
            values = [int(field) for field in fields]
            require(values[0] == index, "eventIdsFile input index not contiguous from zero")
            ordered_events.append((values[1], values[2]))
    require(len(ordered_events) == c["inputEvents"] and len(set(ordered_events)) == len(ordered_events),
            "eventIdsFile must contain exactly inputEvents unique ordered event IDs")
    event_set = set(ordered_events)
    with datasets_path.open(newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        require(reader.fieldnames == ["id", "export_dir"], "datasets TSV needs header id<TAB>export_dir")
        source_rows = list(reader)
    require(source_rows, "no datasets")
    datasets, ids = [], set()
    for row in source_rows:
        name = row["id"]
        require(re.fullmatch(r"[A-Za-z0-9_-]+", name or "") and name not in ids, "invalid/duplicate dataset id")
        ids.add(name)
        folder = Path(row["export_dir"]).expanduser()
        if not folder.is_absolute():
            folder = datasets_path.parent / folder
        dataset = load_dataset(name, folder.resolve(), c, datasets[0] if datasets else None)
        require(set(dataset["event_keys"]).issubset(event_set), "candidate event absent from ordered first inputEvents list")
        dataset["cut_masks"] = cut_masks(dataset, c)
        dataset["regions"] = make_region_masks(dataset["data"], c)
        datasets.append(dataset)
    baseline = datasets[0]
    references = {(split, sb): metrics(baseline, baseline["selected"], split, sb, c)
                  for split in ("full", "train", "validation") for sb in ("default", "near", "far")}
    require(all(x["S"] > 0 and x["Nwin"] > 0 for x in references.values()), "baseline has nonpositive reference S/N")
    output.mkdir(parents=True, exist_ok=False)
    write_json(output / "status.json", {"complete": False, "phase": "scanning"})
    (output / "configuration.yaml").write_text(Path(config_path).read_text())
    source_manifest = {"config": {"path": str(config_path.resolve()), "sha256": sha256(config_path)},
                       "scanner": {"path": str(Path(__file__).resolve()), "sha256": sha256(__file__)},
                       "datasetsTsv": {"path": str(datasets_path.resolve()), "sha256": sha256(datasets_path)},
                       "eventIdsFile": {"path": str(event_path.resolve()), "sha256": sha256(event_path)},
                       "python": sys.version, "numpy": np.__version__, "yaml": yaml.__version__, "datasets": []}
    for d in datasets:
        source_manifest["datasets"].append({"id": d["id"], "upstreamChanges": d["upstream_changes"],
            "meta": d["meta"], "effectiveCuts": d["cuts"], "eventPolicy": d["event_meta"],
            "inputs": {key: {"path": str(path), "sha256": sha256(path)} for key, path in d["paths"].items()}})
    write_json(output / "sources.json", source_manifest)
    train_rows = []
    all_path = output / "all_metrics.csv"
    with all_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=ROW_FIELDS); writer.writeheader()
        for dataset in datasets:
            for point_id, point in enumerate(points):
                mask = point_mask(dataset, point)
                for split in ("train", "validation", "full"):
                    for sb in ("default", "near", "far"):
                        result = add_reference(metrics(dataset, mask, split, sb, c), references[(split, sb)])
                        row = row_for(dataset, point_id, point, split, sb, result, c)
                        writer.writerow(row)
                        if split == "train" and sb == "default":
                            train_rows.append(row)
            print("Scanned {}: {} bounded points; exact baseline replay passed".format(dataset["id"], len(points)), flush=True)
    selected = choose_points(train_rows, c)
    choices = [{"choice_id": "baseline", "label": "KF Imp5 baseline", "role": "baseline",
                "dataset_id": baseline["id"], "point_id": 0, "retention_floors": [], "cuts": {}}]
    index = {(baseline["id"], 0): choices[0]}
    for floor, row in selected:
        identity = (row["dataset_id"], row["point_id"])
        if identity not in index:
            choice = {"choice_id": "choice_{:02d}".format(len(choices)), "label": "", "role": "selected",
                      "dataset_id": row["dataset_id"], "point_id": row["point_id"], "retention_floors": [],
                      "cuts": dict(points[row["point_id"]])}
            choices.append(choice); index[identity] = choice
        index[identity]["retention_floors"].append(floor)
    for choice in choices[1:]:
        choice["label"] = "S retention >= " + "/".join("{:g}%".format(100*x) for x in choice["retention_floors"])
    write_json(output / "chosen_points.json", {"ranking": c["selectionScore"], "selectionSplit": "train",
        "selectionSideband": "default", "choices": choices,
        "floors": [{"floor": floor, "choice_id": index[(row["dataset_id"], row["point_id"])]["choice_id"],
                    "trainingMetrics": row} for floor, row in selected]})
    selection_rows = [{key: choice[key] for key in ("choice_id", "label", "role", "dataset_id", "point_id")}
                      for choice in choices]
    for row, choice in zip(selection_rows, choices):
        row["retention_floors"] = ";".join(str(x) for x in choice["retention_floors"])
        row["cuts_json"] = json.dumps(choice["cuts"], sort_keys=True)
    write_csv(output / "chosen_selections.csv", ["choice_id", "label", "role", "dataset_id", "point_id", "retention_floors", "cuts_json"], selection_rows)
    by_id = {d["id"]: d for d in datasets}
    chosen_metrics, histograms, native_histograms, memberships, chosen_masks = [], [], [], [], {}
    nbins = int(round((c["massMax"]-c["massMin"])/c["massBinWidth"]))
    for choice in choices:
        d = by_id[choice["dataset_id"]]
        point = points[choice["point_id"]]
        mask = point_mask(d, point); chosen_masks[choice["choice_id"]] = mask
        for split in ("full", "train", "validation"):
            use = mask & d["splits"][split]
            counts = np.bincount(d["bins"][use], minlength=nbins+2)
            norm = int(counts[1:-1].sum())
            if split == "full":
                for bin_id in range(nbins+2):
                    low = c["massMin"]+(bin_id-1)*c["massBinWidth"]
                    native_histograms.append({"selection": choice["choice_id"], "bin": bin_id,
                        "low": low, "high": low+c["massBinWidth"], "count": int(counts[bin_id]),
                        "error": math.sqrt(counts[bin_id])})
            for bin_id in range(1, nbins+1):
                low = c["massMin"]+(bin_id-1)*c["massBinWidth"]
                histograms.append({"choice_id": choice["choice_id"], "dataset_id": d["id"], "point_id": choice["point_id"],
                    "split": split, "bin": bin_id, "low": low, "high": low+c["massBinWidth"], "count": int(counts[bin_id]),
                    "error": math.sqrt(counts[bin_id]), "normalization_count": norm})
            for sb in ("default", "near", "far"):
                result = add_reference(metrics(d, mask, split, sb, c), references[(split, sb)])
                row = row_for(d, choice["point_id"], point, split, sb, result, c)
                row["choice_id"] = choice["choice_id"]; chosen_metrics.append(row)
        for raw_row in np.flatnonzero(mask):
            memberships.append({"choice_id": choice["choice_id"], "dataset_id": d["id"], "row": int(raw_row),
                "runId": int(d["data"]["runId"][raw_row]), "eventId": int(d["data"]["eventId"][raw_row]),
                "protonId": int(d["data"]["protonId"][raw_row]), "pionId": int(d["data"]["pionId"][raw_row]),
                "mass": float(d["data"]["mass"][raw_row])})
    write_csv(output / "chosen_metrics.csv", ["choice_id"]+ROW_FIELDS, chosen_metrics)
    write_csv(output / "histograms.csv", ["choice_id", "dataset_id", "point_id", "split", "bin", "low", "high", "count", "error", "normalization_count"], histograms)
    write_csv(output / "chosen_candidate_memberships.csv", ["choice_id", "dataset_id", "row", "runId", "eventId", "protonId", "pionId", "mass"], memberships)
    write_csv(output / "chosen_histograms.csv", ["selection", "bin", "low", "high", "count", "error"], native_histograms)
    with (output / "chosen_cuts.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["selection", "key", "operator", "value"])
        for choice in choices:
            writer.writerow([choice["choice_id"], "selected", "eq", 1])
            for key, value in choice["cuts"].items():
                writer.writerow([choice["choice_id"], key, c["cuts"][key]["operator"], value])
    write_bootstrap(output, choices, datasets, chosen_masks, c, ordered_events)
    manifest = {key: c[key] for key in ("massMin", "massMax", "massBinWidth", "displayMin", "displayMax", "title",
                "fitSigmaMin", "fitSigmaMax", "fitMeanMin", "fitMeanMax", "fitMinEntries")}
    manifest.update(windowLow=c["window"][0], windowHigh=c["window"][1],
                    fitRanges=";".join("{}:{}".format(*bounds) for bounds in c["fitRanges"]))
    with (output / "plot_manifest.tsv").open("w") as stream:
        for key, value in manifest.items():
            stream.write("{}\t{}\n".format(key, value))
    summary = {
        "baselineExactSelectedAndNativeHistogramReplay": True, "datasets": len(datasets),
        "pointsPerDataset": len(points), "totalPoints": len(points)*len(datasets), "choices": len(choices),
        "baseline": {split: references[(split, "default")] for split in ("full", "train", "validation")},
        "eventSplit": "sha256(seed:runId:eventId), first64bits; common event always same split",
        "inputEventSplitCounts": {"train": sum(event_train(*key, c["splitSeed"], c["trainFraction"]) for key in ordered_events),
                                  "validation": sum(not event_train(*key, c["splitSeed"], c["trainFraction"]) for key in ordered_events)},
        "bootstrap": "paired Poisson(1) EVENT cluster bootstrap with shared baseline/choice weights; fixed chosen cuts",
        "scope": "Exploratory finite-grid optimization on TRAIN/default only; validation/full and alternate sidebands never retune. Purity is constant-sideband S/Nwin, not absolute efficiency. No production configuration edits. No pT/eta/TOF/massError/mass-window tuning. Post-KF nSigma/DCA cuts are not asserted identical to upstream reruns.",
        "limitations": "Training score uses independent-Poisson window error; event bootstrap is separate and retains candidate correlations. Background-shape systematic, detector efficiency, feeddown, finite-sample selection bias and uncertainty of the selection procedure are not covered. A holdout from same10000 is an internal check, not independent-data confirmation. Metadata counters cannot prove ordered event identity; preserve source-run provenance.",
    }
    write_json(output / "summary.json", summary)
    write_json(output / "status.json", {"complete": True, "phase": "complete"})
    products = sorted(path for path in output.iterdir() if path.is_file())
    require(all(path.stat().st_size > 0 for path in products), "empty output product")
    with (output / "SHA256SUMS").open("w") as stream:
        for path in products:
            stream.write("{}  {}\n".format(sha256(path), path.name))
    print("PASS: training-only choices, holdout/full metrics, event bootstrap and provenance saved to " + str(output))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--datasets", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        run(args.config, args.datasets, args.output)
    except (ValueError, OSError, KeyError, csv.Error, yaml.YAMLError) as error:
        print("ERROR: " + str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

