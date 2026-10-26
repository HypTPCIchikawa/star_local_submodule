#!/usr/bin/env python3
"""Full-sample exploratory refinement near 70% retained Lambda signal.

This separate tool leaves the initial train-only study unchanged. Reuses its
validated export, exact baseline replay and paired event-bootstrap helpers.
The old train/validation partitions have already been inspected: they are
stability diagnostics ONLY and MUST NOT be described as a fresh holdout.

env -u PYTHONPATH /usr/bin/python3 common/macro/refineLambdaKfPurity.py \
  --config config/qa/lambda_kf_purity70_refinement_20260927.yaml \
  --datasets datasets.tsv --output NEW_DIRECTORY
"""
import argparse
import csv
from functools import lru_cache
import itertools
import json
import math
from pathlib import Path
import re
import sys

import numpy as np
import yaml
import studyLambdaKfPurity as base

SIDEBANDS = ("default", "near", "far")
SPLITS = ("full", "train", "validation")
EXTRA_FIELDS = ["poissonUpperMean", "backgroundUpper", "conservativePurity", "robustPurityScore",
                "minSignalRetentionAcrossSidebands", "splitRole"]
ROW_FIELDS = base.ROW_FIELDS + EXTRA_FIELDS


def load_config(path):
    policy = yaml.load(Path(path).read_text(), Loader=base.UniqueKeyLoader)
    required = {"schemaVersion", "baseConfig", "rankingSplit", "oldSplitsRole", "rankingMetric",
                "retentionMetric", "poissonUpperConfidenceLevel", "signalRetentionFloors",
                "nominalPurityRetentionFloor", "minimumSidebandCount", "minimumWindowCount",
                "maxPointsPerDataset", "includeAllPairs", "tripleAnchor", "grid", "fixedReference", "backgroundCheck"}
    base.require(isinstance(policy, dict) and set(policy) == required, "missing/unknown refinement keys")
    base.require(policy["schemaVersion"] == 1 and policy["rankingSplit"] == "full", "refinement must be full-sample schema1")
    base.require(policy["oldSplitsRole"] == "previously_inspected_stability_diagnostics_only", "old splits are NOT a fresh holdout")
    base.require(policy["rankingMetric"] == "worst_sideband_poisson_upper_background_purity", "unsupported ranking")
    base.require(policy["retentionMetric"] == "minimum_background_subtracted_signal_retention_across_sidebands", "unsupported retention definition")
    base.require(policy["minimumSidebandCount"] == 0, "refinement must not censor sparse sidebands")
    for key in ("minimumWindowCount", "maxPointsPerDataset"):
        base.require(type(policy[key]) is int and policy[key] > 0, key + " must be positive integer")
    cl = policy["poissonUpperConfidenceLevel"]
    base.require(base.finite(cl) and .5 < cl < 1, "invalid one-sided Poisson confidence level")
    floors = policy["signalRetentionFloors"]
    base.require(isinstance(floors, list) and floors and all(base.finite(v) and 0 < v <= 1 for v in floors)
                 and len(set(floors)) == len(floors), "invalid retention floors")
    base.require(base.finite(policy["nominalPurityRetentionFloor"]) and 0 < policy["nominalPurityRetentionFloor"] <= 1,
                 "invalid nominal purity retention floor")
    base_path = Path(policy["baseConfig"])
    if not base_path.is_absolute():
        base_path = Path(__file__).resolve().parents[2] / base_path
    c = base.load_config(base_path)
    grid = policy["grid"]
    base.require(isinstance(grid, dict) and grid and set(grid).issubset(base.CUT_FIELDS), "unknown/empty focused axes")
    cuts = {}
    for key, values in grid.items():
        base.require(isinstance(values, list) and values and all(base.finite(v) and v > 0 for v in values)
                     and len(set(values)) == len(values), "invalid focused levels: " + key)
        if key == "minCosPointing":
            base.require(all(v <= 1 for v in values), "pointing cosine >1")
        field, operator = base.CUT_FIELDS[key]
        cuts[key] = {"field": field, "operator": operator, "values": values}
    anchor = policy["tripleAnchor"]
    base.require(anchor in grid, "tripleAnchor absent from grid")
    base.require(type(policy["includeAllPairs"]) is bool, "includeAllPairs must be boolean")
    other = [key for key in grid if key != anchor]
    c.update(cuts=cuts, includeAllPairs=policy["includeAllPairs"],
             tripleFamilies=[[anchor, a, b] for a, b in itertools.combinations(other, 2)],
             maxPointsPerDataset=policy["maxPointsPerDataset"])
    reference = policy["fixedReference"]
    base.require(isinstance(reference, dict) and set(reference) == {"id", "label", "cuts"}, "invalid fixed reference")
    base.require(re.fullmatch(r"[A-Za-z0-9_-]+", reference["id"] or "") and reference["id"] != "baseline", "invalid reference id")
    base.require(isinstance(reference["label"], str) and reference["label"], "invalid reference label")
    base.require(isinstance(reference["cuts"], dict) and reference["cuts"], "empty reference cuts")
    for key, value in reference["cuts"].items():
        base.require(key in grid and value in grid[key], "reference not contained in grid")
    bg = policy["backgroundCheck"]
    bg_fields = {"backgroundDegrees", "backgroundProfileDeltaNll",
                 "backgroundProfileMaxFactor", "backgroundMinCounts"}
    base.require(isinstance(bg, dict) and set(bg) == bg_fields, "invalid background-check fields")
    base.require(bg["backgroundDegrees"] == "0;1;2", "unsupported background polynomial degrees")
    for key in bg_fields - {"backgroundDegrees"}:
        base.require(base.finite(bg[key]) and bg[key] >= 0, "invalid " + key)
    base.require(bg["backgroundProfileDeltaNll"] > 0 and bg["backgroundProfileMaxFactor"] > 1, "invalid background profile search")
    return policy, c, base_path


def poisson_cdf(k, mean):
    """P[N <= k | mean], evaluated in log space without SciPy."""
    base.require(type(k) is int and k >= 0 and math.isfinite(mean) and mean >= 0, "invalid Poisson parameters")
    if mean == 0:
        return 1.
    log_mean = math.log(mean)
    terms = [-mean + index*log_mean - math.lgamma(index+1) for index in range(k+1)]
    peak = max(terms)
    return min(1., math.exp(peak) * sum(math.exp(term-peak) for term in terms))


@lru_cache(maxsize=None)
def poisson_upper(k, confidence):
    """Exact one-sided Poisson mean upper limit, CDF(k; upper)=1-CL.

    Used as a sparse-background ranking proxy, NOT a confidence interval for
    purity: Nwin uncertainty, correlated candidates and model error are not
    covered by this single Poisson bound.
    """
    base.require(type(k) is int and k >= 0 and .5 < confidence < 1, "invalid Poisson upper-limit request")
    target = 1.-confidence
    low, high = 0., max(1., float(k+1))
    while poisson_cdf(k, high) > target:
        high *= 2.
    for iteration in range(80):
        midpoint = (low+high)/2
        if poisson_cdf(k, midpoint) > target:
            low = midpoint
        else:
            high = midpoint
    return high


def point_results(dataset, mask, references, c, policy):
    results, aggregates = {}, {}
    for split in SPLITS:
        scores, retentions = [], []
        for sb in SIDEBANDS:
            result = base.add_reference(base.metrics(dataset, mask, split, sb, c), references[(split, sb)])
            upper = poisson_upper(int(result["Nleft"]+result["Nright"]), policy["poissonUpperConfidenceLevel"])
            result["poissonUpperMean"] = upper
            result["backgroundUpper"] = result["alpha"]*upper
            result["conservativePurity"] = 1.-result["backgroundUpper"]/result["Nwin"] if result["Nwin"] > 0 else float("nan")
            results[(split, sb)] = result
            scores.append(result["conservativePurity"]); retentions.append(result["signalRetention"])
        robust = min(scores) if all(math.isfinite(x) for x in scores) else float("nan")
        retention = min(retentions) if all(math.isfinite(x) for x in retentions) else float("nan")
        aggregates[split] = {"robustPurityScore": robust, "minSignalRetentionAcrossSidebands": retention}
    return results, aggregates


def make_row(dataset, point_id, point, split, sb, result, aggregate, policy):
    row = dict(result)
    row.update(dataset_id=dataset["id"], point_id=point_id, split=split, sideband=sb,
               sourceChangedParameters=0, addedCutCount=len(point), changedParameters=len(point),
               cuts_json=json.dumps(dict(point), sort_keys=True), score=aggregate["robustPurityScore"], **aggregate)
    row["splitRole"] = "full_sample_exploratory_ranking" if split == "full" else policy["oldSplitsRole"]
    reasons = []
    if split != "full": reasons.append("previously_inspected_diagnostic_split")
    if result["Nwin"] < policy["minimumWindowCount"]: reasons.append("low_window_count")
    if result["S"] <= 0: reasons.append("nonpositive_signal")
    if not math.isfinite(row["score"]): reasons.append("nonfinite_score")
    if not math.isfinite(row["minSignalRetentionAcrossSidebands"]): reasons.append("nonfinite_retention")
    row["selectionEligible"] = int(not reasons)
    row["eligibilityReason"] = ";".join(reasons) if reasons else "eligible_before_all_sideband_retention_floor"
    return row


def select_row(rows, floor, nominal=False):
    eligible = [row for row in rows if row["selectionEligible"] and row["minSignalRetentionAcrossSidebands"] >= floor]
    base.require(eligible, "no eligible focused point for retention " + str(floor))
    key = "purity" if nominal else "robustPurityScore"
    return min(eligible, key=lambda row: (-row[key], row["changedParameters"],
                                          -row["S"], row["point_id"]))


def read_events(c):
    path = Path(c["eventIdsFile"])
    if not path.is_absolute(): path = Path(__file__).resolve().parents[2] / path
    events = []
    with path.open() as source:
        for index, line in enumerate(source):
            fields = line.rstrip("\n").split("\t")
            base.require(len(fields) == 3, "malformed input event TSV")
            values = [int(value) for value in fields]
            base.require(values[0] == index, "unordered input event indices")
            events.append((values[1], values[2]))
    base.require(len(events) == c["inputEvents"] and len(set(events)) == len(events), "wrong/duplicate input event IDs")
    return events, path


def load_input(datasets_path, c, events):
    with datasets_path.open(newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        base.require(reader.fieldnames == ["id", "export_dir"], "expected id/export_dir TSV")
        rows = list(reader)
    base.require(len(rows) == 1, "focused refinement requires the single preserved baseline source")
    row = rows[0]
    base.require(re.fullmatch(r"[A-Za-z0-9_-]+", row["id"] or ""), "invalid baseline dataset id")
    folder = Path(row["export_dir"])
    if not folder.is_absolute(): folder = datasets_path.parent / folder
    d = base.load_dataset(row["id"], folder.resolve(), c)
    base.require(set(d["event_keys"]).issubset(set(events)), "candidate event not in first inputEvents")
    d["cut_masks"] = base.cut_masks(d, c)
    d["regions"] = base.make_region_masks(d["data"], c)
    return d


def export_choices(output, choices, d, points, c, policy, references, events):
    selection_rows, metrics_rows, histograms, native, memberships, masks = [], [], [], [], [], {}
    nbins = int(round((c["massMax"]-c["massMin"])/c["massBinWidth"]))
    with (output / "chosen_cuts.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t"); writer.writerow(["selection", "key", "operator", "value"])
        for choice in choices:
            point = points[choice["point_id"]]
            mask = base.point_mask(d, point); masks[choice["choice_id"]] = mask
            results, aggregates = point_results(d, mask, references, c, policy)
            writer.writerow([choice["choice_id"], "selected", "eq", 1])
            for key, value in point: writer.writerow([choice["choice_id"], key, c["cuts"][key]["operator"], value])
            selection_rows.append({key: choice[key] for key in ("choice_id", "label", "role", "dataset_id", "point_id")})
            selection_rows[-1].update(retention_floors=";".join(str(x) for x in choice["retention_floors"]),
                                     cuts_json=json.dumps(dict(point), sort_keys=True))
            for split in SPLITS:
                use = mask & d["splits"][split]
                counts = np.bincount(d["bins"][use], minlength=nbins+2)
                norm = int(counts[1:-1].sum())
                for bin_id in range(nbins+2):
                    low = c["massMin"]+(bin_id-1)*c["massBinWidth"]
                    row = {"selection": choice["choice_id"], "bin": bin_id, "low": low,
                           "high": low+c["massBinWidth"], "count": int(counts[bin_id]), "error": math.sqrt(counts[bin_id])}
                    if split == "full": native.append(dict(row))
                    if 1 <= bin_id <= nbins:
                        row.pop("selection")
                        row.update(choice_id=choice["choice_id"], dataset_id=d["id"], point_id=choice["point_id"], split=split,
                                   normalization_count=norm)
                        histograms.append(row)
                for sb in SIDEBANDS:
                    row = make_row(d, choice["point_id"], point, split, sb, results[(split, sb)], aggregates[split], policy)
                    row["choice_id"] = choice["choice_id"]; metrics_rows.append(row)
            for index in np.flatnonzero(mask):
                memberships.append({"choice_id": choice["choice_id"], "dataset_id": d["id"], "row": int(index),
                    "runId": int(d["data"]["runId"][index]), "eventId": int(d["data"]["eventId"][index]),
                    "protonId": int(d["data"]["protonId"][index]), "pionId": int(d["data"]["pionId"][index]),
                    "mass": float(d["data"]["mass"][index])})
    base.write_csv(output / "chosen_selections.csv", ["choice_id", "label", "role", "dataset_id", "point_id", "retention_floors", "cuts_json"], selection_rows)
    base.write_csv(output / "chosen_metrics.csv", ["choice_id"]+ROW_FIELDS, metrics_rows)
    base.write_csv(output / "histograms.csv", ["choice_id", "dataset_id", "point_id", "split", "bin", "low", "high", "count", "error", "normalization_count"], histograms)
    base.write_csv(output / "chosen_histograms.csv", ["selection", "bin", "low", "high", "count", "error"], native)
    base.write_csv(output / "chosen_candidate_memberships.csv", ["choice_id", "dataset_id", "row", "runId", "eventId", "protonId", "pionId", "mass"], memberships)
    base.write_bootstrap(output, choices, [d], masks, c, events)


def run(config_path, datasets_path, output):
    policy, c, base_path = load_config(config_path)
    base.require(not output.exists(), "output directory exists; refusing overwrite")
    points = base.enumerate_points(c)
    ref_point = tuple(sorted((key, float(value)) for key, value in policy["fixedReference"]["cuts"].items()))
    base.require(ref_point in points, "fixed previous reference absent from enumeration")
    events, event_path = read_events(c)
    d = load_input(datasets_path, c, events)
    references = {(split, sb): base.metrics(d, d["selected"], split, sb, c) for split in SPLITS for sb in SIDEBANDS}
    base.require(all(r["S"] > 0 and r["Nwin"] > 0 for r in references.values()), "nonpositive baseline reference")
    output.mkdir(parents=True, exist_ok=False)
    base.write_json(output / "status.json", {"complete": False, "phase": "full_sample_refinement"})
    (output / "configuration.yaml").write_text(Path(config_path).read_text())
    (output / "resolved_base_configuration.yaml").write_text(yaml.safe_dump(c, sort_keys=False))
    sources = {"policy": {"path": str(config_path.resolve()), "sha256": base.sha256(config_path)},
        "baseConfig": {"path": str(base_path.resolve()), "sha256": base.sha256(base_path)},
        "scanner": {"path": str(Path(__file__).resolve()), "sha256": base.sha256(__file__)},
        "helpers": {"path": str(Path(base.__file__).resolve()), "sha256": base.sha256(base.__file__)},
        "eventIdsFile": {"path": str(event_path.resolve()), "sha256": base.sha256(event_path)},
        "datasetsTsv": {"path": str(datasets_path.resolve()), "sha256": base.sha256(datasets_path)},
        "dataset": {"id": d["id"], "meta": d["meta"], "effectiveCuts": d["cuts"], "eventPolicy": d["event_meta"],
            "inputs": {key: {"path": str(path), "sha256": base.sha256(path)} for key, path in d["paths"].items()}},
        "python": sys.version, "numpy": np.__version__, "yaml": yaml.__version__}
    base.write_json(output / "sources.json", sources)
    ranking_rows = []
    with (output / "all_metrics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=ROW_FIELDS); writer.writeheader()
        for point_id, point in enumerate(points):
            mask = base.point_mask(d, point)
            results, aggregates = point_results(d, mask, references, c, policy)
            for split in SPLITS:
                for sb in SIDEBANDS:
                    row = make_row(d, point_id, point, split, sb, results[(split, sb)], aggregates[split], policy)
                    writer.writerow(row)
                    if split == "full" and sb == "default": ranking_rows.append(row)
    base.write_csv(output / "ranking.csv", ROW_FIELDS, ranking_rows)
    choices, by_point, decisions = [], {}, []
    def add_choice(point_id, choice_id, label, role, floor=None):
        if point_id not in by_point:
            entry = {"choice_id": choice_id, "dataset_id": d["id"], "point_id": point_id,
                     "label": label, "role": role, "retention_floors": [], "cuts": dict(points[point_id]), "uses": []}
            choices.append(entry); by_point[point_id] = entry
        entry = by_point[point_id]
        entry["uses"].append(choice_id)
        if floor is not None: entry["retention_floors"].append(floor)
        return entry["choice_id"]
    add_choice(0, "baseline", "KF Imp5 baseline", "baseline")
    fixed = policy["fixedReference"]
    add_choice(points.index(ref_point), fixed["id"], fixed["label"], "fixed_previous_reference")
    for floor in policy["signalRetentionFloors"]:
        winner = select_row(ranking_rows, floor)
        name = "robust{:g}".format(100*floor)
        choice_id = add_choice(winner["point_id"], name, "Robust S retention >= {:g}%".format(100*floor), "full_ranked", floor)
        decisions.append({"criterion": "robust", "retentionFloor": floor, "choice_id": choice_id, "fullMetrics": winner})
    floor = policy["nominalPurityRetentionFloor"]
    winner = select_row(ranking_rows, floor, nominal=True)
    choice_id = add_choice(winner["point_id"], "nominal{:g}".format(100*floor), "Highest nominal purity (S >= {:g}%)".format(100*floor), "nominal_full_ranked", floor)
    decisions.append({"criterion": "nominal_default_purity", "retentionFloor": floor, "choice_id": choice_id, "fullMetrics": winner})
    base.write_json(output / "chosen_points.json", {"selectionSplit": "full", "oldSplitsRole": policy["oldSplitsRole"],
        "ranking": policy["rankingMetric"], "retentionDefinition": policy["retentionMetric"], "choices": choices, "decisions": decisions})
    export_choices(output, choices, d, points, c, policy, references, events)
    manifest = {key: c[key] for key in ("massMin", "massMax", "massBinWidth", "displayMin", "displayMax", "title",
        "fitSigmaMin", "fitSigmaMax", "fitMeanMin", "fitMeanMax", "fitMinEntries")}
    manifest.update(windowLow=c["window"][0], windowHigh=c["window"][1],
                    fitRanges=";".join("{}:{}".format(*bounds) for bounds in c["fitRanges"]))
    manifest.update(policy["backgroundCheck"])
    manifest["backgroundRegions"] = ";".join("{}:{}:{}:{}:{}".format(sb, *c["sidebands"][sb]["left"], *c["sidebands"][sb]["right"]) for sb in SIDEBANDS)
    with (output / "plot_manifest.tsv").open("w") as stream:
        for key, value in manifest.items(): stream.write("{}\t{}\n".format(key, value))
    summary = {"baselineExactSelectedAndNativeHistogramReplay": True, "pointsPerDataset": len(points), "choices": len(choices),
        "inputEvents": c["inputEvents"], "reconstructedEvents": c["reconstructedEvents"],
        "selectionSplit": "full", "oldSplitsRole": policy["oldSplitsRole"],
        "baseline": {split: references[(split, "default")] for split in SPLITS},
        "scope": "Exploratory full-sample refinement on previously inspected first10000 events. Old train/validation are descriptive stability partitions, NOT a fresh holdout. No production changes, mass/window/acceptance optimization or absolute efficiency estimate.",
        "ranking": "Minimum over default/near/far of 1-alpha*PoissonUpper(nSideband,CL)/Nwin; eligibility uses minimum same-sideband S/Sbaseline. No sparse-sideband censor; zero counts retain positive upper background. nominal default-purity champion uses the SAME eligibility.",
        "limitations": "Conservative score is a ranking proxy, NOT a confidence bound for purity: fixed Nwin, Poisson background approximation and mass-shape assumptions. Overlapping sideband scores are not independent evidence. Shared-event bootstrap keeps chosen cuts fixed and does not cover tuning bias, selection-process uncertainty, background-model systematics, or detector efficiency. Independent data/embedding needed for confirmation.",
        "poissonUpperConfidenceLevel": policy["poissonUpperConfidenceLevel"],
        "poissonZeroCountUpper": poisson_upper(0, policy["poissonUpperConfidenceLevel"])}
    base.write_json(output / "summary.json", summary)
    base.write_json(output / "status.json", {"complete": True, "phase": "complete"})
    products = sorted(path for path in output.iterdir() if path.is_file())
    base.require(all(path.stat().st_size > 0 for path in products), "empty refinement product")
    with (output / "SHA256SUMS").open("w") as stream:
        for path in products: stream.write("{}  {}\n".format(base.sha256(path), path.name))
    print("PASS full-sample EXPLORATORY refinement: {} points, {} unique choices; prior splits are diagnostics only".format(len(points), len(choices)))
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

