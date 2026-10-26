#!/usr/bin/env python3
"""Joint FULL-sample exploratory ranking of real upstream KF reruns.

The first dataset is the common denominator, never each variant's own yield.
All input events are used for discovery; there is no training/holdout claim.
Run independent, frozen-cut data afterwards. No PicoDst event loop is here.

env -u PYTHONPATH /usr/bin/python3 common/macro/scanLambdaKfUpstreamPurity.py \
  --config study.yaml --datasets datasets.tsv --output NEW_DIRECTORY
"""
import argparse
import csv
import itertools
import json
import math
from pathlib import Path
import re
import sys

import numpy as np
import yaml
import studyLambdaKfPurity as base
import refineLambdaKfPurity as refinement

SIDEBANDS = ("default", "near", "far")
UPSTREAM_DIRECTIONS = dict(base.UPSTREAM_DIRECTIONS, minNHitsFit="ge", minNHitsRatio="ge")
ROW_FIELDS = base.ROW_FIELDS + ["poissonUpperMean", "backgroundUpper", "conservativePurity",
                              "robustPurityScore", "minSignalRetentionAcrossSidebands", "splitRole"]
BACKGROUND_DEFAULTS = {"backgroundDegrees": "0;1;2", "backgroundProfileDeltaNll": 0.821187207574908,
                       "backgroundProfileMaxFactor": 100, "backgroundMinCounts": 8}


def repo_path(value):
    path = Path(value)
    return path if path.is_absolute() else Path(__file__).resolve().parents[2] / path


def load_config(path):
    policy = yaml.load(Path(path).read_text(), Loader=base.UniqueKeyLoader)
    required = {"schemaVersion", "baseConfig", "eventIdsFile", "inputEvents", "reconstructedEvents",
                "grid", "includeAllPairs", "tripleAnchor", "retentionFloors", "minimumWindowCount",
                "poissonUpperConfidenceLevel", "maxPointsPerDataset"}
    optional = {"title", "bootstrapReplicates", "bootstrapSeed", "backgroundCheck"}
    base.require(isinstance(policy, dict) and required.issubset(policy) and
                 set(policy).issubset(required | optional), "missing/unknown upstream-study keys")
    base.require(policy["schemaVersion"] == 1, "unsupported schemaVersion")
    for key in ("inputEvents", "reconstructedEvents", "minimumWindowCount", "maxPointsPerDataset"):
        base.require(type(policy[key]) is int and policy[key] > 0, key + " must be positive integer")
    base.require(policy["reconstructedEvents"] <= policy["inputEvents"], "invalid event counts")
    base.require(isinstance(policy["eventIdsFile"], str) and policy["eventIdsFile"], "missing event IDs path")
    floors = policy["retentionFloors"]
    base.require(isinstance(floors, list) and floors and all(base.finite(x) and 0 < x <= 1 for x in floors)
                 and len(set(floors)) == len(floors), "invalid retentionFloors")
    cl = policy["poissonUpperConfidenceLevel"]
    base.require(base.finite(cl) and .5 < cl < 1, "invalid Poisson confidence level")
    base_path = repo_path(policy["baseConfig"])
    c = base.load_config(base_path)
    grid = policy["grid"]
    base.require(isinstance(grid, dict) and grid and set(grid).issubset(base.CUT_FIELDS), "unknown/empty grid")
    cuts = {}
    for key, levels in grid.items():
        base.require(isinstance(levels, list) and levels and all(base.finite(x) and x > 0 for x in levels)
                     and len(set(levels)) == len(levels), "invalid levels for " + key)
        if key == "minCosPointing": base.require(all(x <= 1 for x in levels), "cosine exceeds one")
        field, op = base.CUT_FIELDS[key]
        cuts[key] = {"field": field, "operator": op, "values": levels}
    anchor = policy["tripleAnchor"]
    base.require(anchor in cuts, "tripleAnchor absent from grid")
    base.require(type(policy["includeAllPairs"]) is bool, "includeAllPairs must be boolean")
    others = [key for key in cuts if key != anchor]
    c.update(cuts=cuts, includeAllPairs=policy["includeAllPairs"],
             tripleFamilies=[[anchor, a, b] for a, b in itertools.combinations(others, 2)],
             maxPointsPerDataset=policy["maxPointsPerDataset"], eventIdsFile=policy["eventIdsFile"],
             inputEvents=policy["inputEvents"], reconstructedEvents=policy["reconstructedEvents"],
             bootstrapReplicates=policy.get("bootstrapReplicates", 400),
             bootstrapSeed=policy.get("bootstrapSeed", c["bootstrapSeed"]))
    base.require(type(c["bootstrapReplicates"]) is int and c["bootstrapReplicates"] >= 200,
                 "at least 200 bootstrap replicates required")
    base.require(type(c["bootstrapSeed"]) is int and 0 <= c["bootstrapSeed"] < 2**32, "invalid bootstrap seed")
    c["title"] = policy.get("title", c["title"])
    base.require(isinstance(c["title"], str) and c["title"] and not any(x in c["title"] for x in "\t\n\r"),
                 "invalid title")
    bg = policy.get("backgroundCheck", dict(BACKGROUND_DEFAULTS))
    base.require(isinstance(bg, dict) and set(bg) == set(BACKGROUND_DEFAULTS), "invalid backgroundCheck keys")
    base.require(bg["backgroundDegrees"] == "0;1;2", "unsupported background degrees")
    for key in set(bg) - {"backgroundDegrees"}:
        base.require(base.finite(bg[key]) and bg[key] >= 0, "invalid background setting " + key)
    base.require(bg["backgroundProfileDeltaNll"] > 0 and bg["backgroundProfileMaxFactor"] > 1,
                 "invalid background profile range")
    policy["backgroundCheck"] = bg
    return policy, c, base_path


def audit_source_differences(dataset, reference):
    """No helper-global monkeypatch: preserve old studies and audit explicitly."""
    changes = {}
    cuts, original = dataset["cuts"], reference["cuts"]
    base.require(set(cuts) == set(original), "different effective metadata key set")
    for key, value in cuts.items():
        if value == original[key]: continue
        base.require(key in UPSTREAM_DIRECTIONS, "unapproved effective difference: " + key)
        old = original[key]
        base.require(base.finite(value) and base.finite(old), "nonnumeric upstream difference: " + key)
        base.require(value >= old if UPSTREAM_DIRECTIONS[key] == "ge" else value <= old,
                     "upstream variant must tighten: " + key)
        if key == "minNHitsFit": base.require(value == int(value) and value >= 0, "invalid hit count")
        if key == "minNHitsRatio": base.require(0 <= value <= 1, "invalid hit ratio")
        changes[key] = {"baseline": old, "value": value}
    ignored = {"KFEventSelection.mainconf", "KFEventSelection.analysisInfo", "KFEventSelection.eventConfig"}
    event = lambda d: {k: v for k, v in d["event_meta"].items() if k not in ignored}
    base.require(event(dataset) == event(reference), "different effective event policies")
    return changes


def read_events(c):
    """Accept one explicit ledger header or the legacy headerless format."""
    path = repo_path(c["eventIdsFile"])
    events = []
    with path.open(newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        for line_number, fields in enumerate(reader):
            if line_number == 0 and fields == ["index", "runId", "eventId"]:
                continue
            base.require(len(fields) == 3, "malformed input event TSV")
            values = [int(field) for field in fields]
            base.require(values[0] == len(events), "unordered input event indices")
            events.append((values[1], values[2]))
    base.require(len(events) == c["inputEvents"] and len(set(events)) == len(events),
                 "wrong/duplicate input event IDs")
    return events, path


def load_inputs(path, c, events):
    with Path(path).open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        base.require(reader.fieldnames == ["id", "export_dir"], "datasets TSV needs id/export_dir header")
        rows = list(reader)
    base.require(rows, "empty datasets list")
    result, ids, event_set = [], set(), set(events)
    for row in rows:
        name = row["id"]
        base.require(re.fullmatch(r"[A-Za-z0-9_-]+", name or "") and name not in ids, "invalid/duplicate dataset id")
        ids.add(name)
        folder = Path(row["export_dir"])
        if not folder.is_absolute(): folder = Path(path).parent / folder
        # Reuse full source, selected-bit and exact native-histogram validation.
        # Cross-source audit below additionally permits two hit-quality gates.
        dataset = base.load_dataset(name, folder.resolve(), c)
        dataset["upstream_changes"] = audit_source_differences(dataset, result[0]) if result else {}
        base.require(set(dataset["event_keys"]).issubset(event_set), "candidate event absent from input event IDs")
        dataset["regions"] = base.make_region_masks(dataset["data"], c)
        dataset["cut_masks"] = base.cut_masks(dataset, c)
        result.append(dataset)
    return result


def point_results(dataset, mask, reference, c, policy):
    results = {}
    for sb in SIDEBANDS:
        value = base.add_reference(base.metrics(dataset, mask, "full", sb, c), reference[sb])
        upper = refinement.poisson_upper(int(value["Nleft"] + value["Nright"]), policy["poissonUpperConfidenceLevel"])
        value.update(poissonUpperMean=upper, backgroundUpper=value["alpha"] * upper)
        value["conservativePurity"] = 1 - value["backgroundUpper"] / value["Nwin"] if value["Nwin"] else float("nan")
        results[sb] = value
    def minimum(key):
        values = [results[sb][key] for sb in SIDEBANDS]
        return min(values) if all(math.isfinite(x) for x in values) else float("nan")
    aggregate = {"robustPurityScore": minimum("conservativePurity"),
                 "minSignalRetentionAcrossSidebands": minimum("signalRetention")}
    return results, aggregate


def make_row(dataset, point_id, point, sb, value, aggregate, policy):
    row = dict(value, **aggregate)
    nup = len(dataset["upstream_changes"])
    row.update(dataset_id=dataset["id"], point_id=point_id, split="full", sideband=sb,
               sourceChangedParameters=nup, addedCutCount=len(point), changedParameters=nup + len(point),
               cuts_json=json.dumps(dict(point), sort_keys=True), score=aggregate["robustPurityScore"],
               splitRole="full_sample_exploratory_ranking")
    reasons = []
    if value["Nwin"] < policy["minimumWindowCount"]: reasons.append("low_window_count")
    if aggregate["minSignalRetentionAcrossSidebands"] <= 0: reasons.append("nonpositive_signal")
    if not math.isfinite(row["score"]): reasons.append("nonfinite_score")
    if not math.isfinite(row["minSignalRetentionAcrossSidebands"]): reasons.append("nonfinite_retention")
    row.update(selectionEligible=int(not reasons),
               eligibilityReason=";".join(reasons) if reasons else "eligible_before_retention_floor")
    return row


def select_row(rows, floor, nominal=False):
    eligible = [r for r in rows if r["selectionEligible"] and r["minSignalRetentionAcrossSidebands"] >= floor]
    base.require(eligible, "no eligible point at retention floor " + str(floor))
    metric = "purity" if nominal else "robustPurityScore"
    return min(eligible, key=lambda r: (-r[metric], r["changedParameters"], -r["S"], r["dataset_id"], r["point_id"]))


def write_bootstrap(output, choices, datasets, masks, c, events):
    """Shared full-event Poisson weights, including zero-candidate events."""
    lookup = {event: i for i, event in enumerate(events)}
    by_id = {d["id"]: d for d in datasets}
    labels = ["signal"] + [sb + "_" + side for sb in SIDEBANDS for side in ("left", "right")]
    counts = np.zeros((len(choices), len(events), len(labels)), dtype=np.int64)
    for ci, choice in enumerate(choices):
        d = by_id[choice["dataset_id"]]
        indices = np.array([lookup[event] for event in d["event_keys"]], dtype=np.int64)
        for ri, region in enumerate(labels):
            use = masks[choice["choice_id"]] & d["regions"][region]
            counts[ci, :, ri] = np.bincount(indices[use], minlength=len(events))
    np.savez_compressed(output / "chosen_event_counts.npz", event_ids=np.asarray(events, dtype=np.int64),
                        choice_ids=np.array([x["choice_id"] for x in choices]), region_labels=np.array(labels), counts=counts)
    rng = np.random.RandomState(c["bootstrapSeed"])
    weights = rng.poisson(1, size=(c["bootstrapReplicates"], len(events)))
    sums = weights.dot(counts.transpose(1, 0, 2).reshape(len(events), -1))
    sums = sums.reshape(c["bootstrapReplicates"], len(choices), len(labels))
    actual = counts.sum(axis=1)
    rows, draws = [], {}
    for si, sb_name in enumerate(SIDEBANDS):
        sb = c["sidebands"][sb_name]
        alpha = (c["window"][1] - c["window"][0]) / (sb["left"][1] - sb["left"][0] + sb["right"][1] - sb["right"][0])
        n = sums[:, :, 0].astype(float)
        signal = n - alpha * (sums[:, :, 1 + 2 * si] + sums[:, :, 2 + 2 * si])
        actual_n = actual[:, 0].astype(float)
        actual_s = actual_n - alpha * (actual[:, 1 + 2 * si] + actual[:, 2 + 2 * si])
        with np.errstate(divide="ignore", invalid="ignore"):
            purity = signal / n
            retention = signal / signal[:, :1]
            retention[signal[:, 0] <= 0, :] = np.nan
            p_est, r_est = actual_s / actual_n, actual_s / actual_s[0]
        for ci, choice in enumerate(choices):
            for quantity, values, estimate in (("purity", purity[:, ci], p_est[ci]),
                    ("purityGain", purity[:, ci] - purity[:, 0], p_est[ci] - p_est[0]),
                    ("signalRetention", retention[:, ci], r_est[ci])):
                good = values[np.isfinite(values)]
                valid = len(good) >= math.ceil(c["bootstrapMinValidFraction"] * c["bootstrapReplicates"])
                quantiles = np.quantile(good, c["bootstrapQuantiles"]) if valid else [float("nan")] * 2
                rows.append(dict(choice_id=choice["choice_id"], split="full", sideband=sb_name, quantity=quantity,
                    estimate=float(estimate), bootstrapError=float(np.std(good, ddof=1)) if valid else float("nan"),
                    quantileLow=float(quantiles[0]), quantileHigh=float(quantiles[1]), validReplicates=len(good),
                    replicates=c["bootstrapReplicates"], status="ok" if valid else "insufficient_finite_replicates"))
                draws["{}_full_{}_{}".format(choice["choice_id"], sb_name, quantity)] = values
    base.write_csv(output / "paired_event_bootstrap.csv", list(rows[0]), rows)
    np.savez_compressed(output / "paired_event_bootstrap_draws.npz", **draws)


def export_choices(output, choices, datasets, points, c, policy, references, events):
    by_id = {d["id"]: d for d in datasets}
    selections, metrics, histograms, native, cuts_rows, source_rows, masks = [], [], [], [], [], [], {}
    nbins = int(round((c["massMax"] - c["massMin"]) / c["massBinWidth"]))
    for choice in choices:
        d = by_id[choice["dataset_id"]]
        point = points[choice["point_id"]]
        mask = base.point_mask(d, point)
        masks[choice["choice_id"]] = mask
        results, aggregate = point_results(d, mask, references, c, policy)
        counts = np.bincount(d["bins"][mask], minlength=nbins + 2)
        norm = int(counts[1:-1].sum())
        selections.append({key: choice[key] for key in ("choice_id", "label", "role", "dataset_id", "point_id")})
        selections[-1].update(retention_floors=";".join(map(str, choice["retention_floors"])),
                              cuts_json=json.dumps(dict(point), sort_keys=True))
        cuts_rows.append([choice["choice_id"], "selected", "eq", 1])
        for key, value in point: cuts_rows.append([choice["choice_id"], key, c["cuts"][key]["operator"], value])
        for sb in SIDEBANDS:
            row = make_row(d, choice["point_id"], point, sb, results[sb], aggregate, policy)
            row["choice_id"] = choice["choice_id"]
            metrics.append(row)
        for i in range(nbins + 2):
            low = c["massMin"] + (i - 1) * c["massBinWidth"]
            item = dict(selection=choice["choice_id"], bin=i, low=low, high=low + c["massBinWidth"],
                        count=int(counts[i]), error=math.sqrt(counts[i]))
            native.append(item)
            if 1 <= i <= nbins:
                item = dict(item); item.pop("selection")
                item.update(choice_id=choice["choice_id"], dataset_id=d["id"], point_id=choice["point_id"],
                            split="full", normalization_count=norm)
                histograms.append(item)
    base.write_csv(output / "chosen_selections.csv", list(selections[0]), selections)
    base.write_csv(output / "chosen_metrics.csv", ["choice_id"] + ROW_FIELDS, metrics)
    base.write_csv(output / "histograms.csv", ["choice_id", "dataset_id", "point_id", "split", "bin", "low", "high", "count", "error", "normalization_count"], histograms)
    base.write_csv(output / "chosen_histograms.csv", ["selection", "bin", "low", "high", "count", "error"], native)
    def write_cuts(path, rows):
        with path.open("w", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t")
            writer.writerow(["selection", "key", "operator", "value"]); writer.writerows(rows)
    write_cuts(output / "chosen_cuts.tsv", cuts_rows)
    # Every upstream source needs its OWN ROOT verification: a baseline ROOT
    # cannot reproduce candidates lost by an upstream hit/Finder gate.
    verify = output / "verification"; verify.mkdir()
    for dataset_id in dict.fromkeys(choice["dataset_id"] for choice in choices):
        wanted = {x["choice_id"] for x in choices if x["dataset_id"] == dataset_id}
        cut_path, hist_path = verify / (dataset_id + "_cuts.tsv"), verify / (dataset_id + "_histograms.csv")
        write_cuts(cut_path, [row for row in cuts_rows if row[0] in wanted])
        base.write_csv(hist_path, ["selection", "bin", "low", "high", "count", "error"],
                       [row for row in native if row["selection"] in wanted])
        source_rows.append(dict(dataset_id=dataset_id, inputRoot=by_id[dataset_id]["meta"]["inputRoot"],
                                cuts=str(cut_path.relative_to(output)), histograms=str(hist_path.relative_to(output))))
    base.write_csv(output / "verification_sources.csv", list(source_rows[0]), source_rows)
    write_bootstrap(output, choices, datasets, masks, c, events)


def run(config_path, datasets_path, output):
    config_path, datasets_path, output = Path(config_path), Path(datasets_path), Path(output)
    policy, c, base_path = load_config(config_path)
    base.require(not output.exists(), "output directory exists; refusing overwrite")
    points = base.enumerate_points(c)
    events, event_path = read_events(c)
    datasets = load_inputs(datasets_path, c, events)
    baseline = datasets[0]
    references = {sb: base.metrics(baseline, baseline["selected"], "full", sb, c) for sb in SIDEBANDS}
    base.require(all(r["S"] > 0 and r["Nwin"] > 0 for r in references.values()),
                 "baseline has empty/nonpositive sideband-subtracted signal; cannot rank retention")
    output.mkdir(parents=True, exist_ok=False)
    base.write_json(output / "status.json", {"complete": False, "phase": "full_sample_upstream_scan"})
    (output / "configuration.yaml").write_text(config_path.read_text())
    (output / "resolved_base_configuration.yaml").write_text(yaml.safe_dump(c, sort_keys=False))
    def file_record(path): return {"path": str(Path(path).resolve()), "sha256": base.sha256(path)}
    base.write_json(output / "sources.json", {
        "policy": file_record(config_path), "baseConfig": file_record(base_path), "scanner": file_record(__file__),
        "helpers": [file_record(base.__file__), file_record(refinement.__file__)],
        "datasetsTsv": file_record(datasets_path), "eventIdsFile": file_record(event_path),
        "python": sys.version, "numpy": np.__version__, "yaml": yaml.__version__,
        "datasets": [{"id": d["id"], "meta": d["meta"], "effectiveCuts": d["cuts"], "eventPolicy": d["event_meta"],
                      "upstreamChanges": d["upstream_changes"], "inputs": {k: file_record(p) for k, p in d["paths"].items()}}
                     for d in datasets]})
    ranking = []
    with (output / "all_metrics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=ROW_FIELDS); writer.writeheader()
        for dataset in datasets:
            for i, point in enumerate(points):
                results, aggregate = point_results(dataset, base.point_mask(dataset, point), references, c, policy)
                for sb in SIDEBANDS:
                    row = make_row(dataset, i, point, sb, results[sb], aggregate, policy)
                    writer.writerow(row)
                    if sb == "default": ranking.append(row)
            print("Scanned {}: {} points; common baseline {}".format(dataset["id"], len(points), baseline["id"]), flush=True)
    base.write_csv(output / "ranking.csv", ROW_FIELDS, ranking)
    by_id = {d["id"]: d for d in datasets}
    choices, decisions = [], []
    def choice(name, dataset_id, point_id, label, floor=None, nominal=False):
        d = by_id[dataset_id]
        point = dict(points[point_id])
        overrides = {key: value["value"] for key, value in d["upstream_changes"].items()}
        overrides.update(point)
        choices.append(dict(choice_id=name, dataset_id=dataset_id, point_id=point_id, label=label,
            role="baseline" if name == "baseline" else ("nominal_full_ranked" if nominal else "full_ranked"),
            retention_floors=[] if floor is None else [floor], cuts=point, downstreamCuts=point,
            upstreamChanges=d["upstream_changes"], combinedCutOverrides=overrides))
    choice("baseline", baseline["id"], 0, "KF Imp5 baseline")
    for floor in policy["retentionFloors"]:
        winner = select_row(ranking, floor)
        name = "best{:g}".format(100 * floor)
        choice(name, winner["dataset_id"], winner["point_id"], "Best robust purity (S >= {:g}%)".format(100 * floor), floor)
        decisions.append(dict(criterion="robust", retentionFloor=floor, choice_id=name, fullMetrics=winner))
    floor = min(policy["retentionFloors"])
    winner = select_row(ranking, floor, nominal=True)
    name = "nominal{:g}".format(100 * floor)
    choice(name, winner["dataset_id"], winner["point_id"], "Highest nominal purity (S >= {:g}%)".format(100 * floor), floor, True)
    decisions.append(dict(criterion="nominal_default_purity", retentionFloor=floor, choice_id=name, fullMetrics=winner))
    base.write_json(output / "chosen_points.json", dict(selectionSplit="full", choices=choices, decisions=decisions,
        ranking="worst_sideband_poisson_upper_background_purity", baselineDataset=baseline["id"],
        retentionDefinition="minimum same-sideband background-subtracted S / COMMON baseline S"))
    export_choices(output, choices, datasets, points, c, policy, references, events)
    manifest = {key: c[key] for key in ("massMin", "massMax", "massBinWidth", "displayMin", "displayMax", "title",
                "fitSigmaMin", "fitSigmaMax", "fitMeanMin", "fitMeanMax", "fitMinEntries")}
    manifest.update(windowLow=c["window"][0], windowHigh=c["window"][1],
                    fitRanges=";".join("{}:{}".format(*x) for x in c["fitRanges"]), **policy["backgroundCheck"])
    manifest["backgroundRegions"] = ";".join("{}:{}:{}:{}:{}".format(sb, *c["sidebands"][sb]["left"], *c["sidebands"][sb]["right"]) for sb in SIDEBANDS)
    with (output / "plot_manifest.tsv").open("w") as stream:
        for key, value in manifest.items(): stream.write("{}\t{}\n".format(key, value))
    summary = dict(baselineExactSelectedAndNativeHistogramReplay=True, datasets=len(datasets),
        pointsPerDataset=len(points), totalPoints=len(points) * len(datasets), selectionSplit="full",
        inputEvents=c["inputEvents"], reconstructedEvents=c["reconstructedEvents"], baseline=references,
        baselineDataset=baseline["id"], choices=len(choices),
        poissonZeroCountUpper=refinement.poisson_upper(0, policy["poissonUpperConfidenceLevel"]),
        scope="FULL-sample exploratory upstream/downstream scan; no held-out validation. Freeze cuts before independent data.",
        limitations="Ranking is a sparse-sideband proxy, NOT a purity confidence limit. Fixed-cut paired full-event bootstrap does not cover tuning bias or background-shape systematics. Retention is background-subtracted candidate yield, not detector efficiency. Zero sideband counts are not proof of zero background. Input counters alone cannot prove ordered event identity; retain actual run/input provenance.",
        plotting="Use checkLambdaPurityRefinement.C. Do not use drawLambdaPurityStudy.C: its hardcoded train/holdout annotations do not apply.",
        verification="Use verification_sources.csv and each source-specific cuts/histograms pair; combined chosen_cuts.tsv must NOT be replayed on one baseline ROOT.")
    base.write_json(output / "summary.json", summary)
    base.write_json(output / "status.json", {"complete": True, "phase": "complete"})
    products = sorted(p for p in output.rglob("*") if p.is_file())
    base.require(all(p.stat().st_size > 0 for p in products), "empty study output")
    with (output / "SHA256SUMS").open("w") as stream:
        for path in products: stream.write("{}  {}\n".format(base.sha256(path), path.relative_to(output)))
    print("PASS exploratory upstream scan: {} datasets, {} total points".format(len(datasets), len(points) * len(datasets)))
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
