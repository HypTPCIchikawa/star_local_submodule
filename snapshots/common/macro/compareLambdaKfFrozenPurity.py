#!/usr/bin/env python3
"""Evaluate actual default/frozen-best KF reruns on independent input events.

There is NO grid, extra candidate selection, or winner re-selection here.
Usage: --config TEN_K_SCORER.yaml --baseline EXPORT --best EXPORT
 --event-ids validation_event_ids.tsv --expected-events 100000
 --frozen-choice SCAN/chosen_points.json --output NEW_DIRECTORY
"""
import argparse
import csv
import json
import math
from pathlib import Path
import sys

import numpy as np
import yaml
import studyLambdaKfPurity as base
import scanLambdaKfUpstreamPurity as scan

SIDEBANDS = scan.SIDEBANDS


def audit_frozen_cuts(baseline, best, frozen, source_manifest):
    choices = [x for x in frozen.get("choices", []) if x.get("choice_id") == "best50"]
    base.require(len(choices) == 1, "frozen selection must contain exactly one best50")
    choice = choices[0]
    overrides = choice.get("combinedCutOverrides")
    base.require(isinstance(overrides, dict), "missing combinedCutOverrides in frozen best50")
    allowed = set(scan.UPSTREAM_DIRECTIONS) | set(base.CUT_FIELDS)
    base.require(set(overrides).issubset(allowed) and all(base.finite(v) for v in overrides.values()),
                 "unknown/nonnumeric frozen cut override")
    old, actual = baseline["cuts"], best["cuts"]
    base.require(set(old) == set(actual), "different default/best effective metadata keys")
    base.require(set(overrides).issubset(old), "frozen cut absent from actual effective configuration")
    expected = dict(old); expected.update(overrides)
    base.require(actual == expected, "actual best effective cuts do not exactly match frozen overrides")
    ignored = {"KFEventSelection.mainconf", "KFEventSelection.analysisInfo", "KFEventSelection.eventConfig"}
    policy = lambda d: {k: v for k, v in d["event_meta"].items() if k not in ignored}
    base.require(policy(baseline) == policy(best), "different default/best effective event policies")
    baseline_id = frozen.get("baselineDataset")
    sources = [d for d in source_manifest.get("datasets", []) if d.get("id") == baseline_id]
    base.require(len(sources) == 1, "frozen source manifest lacks unique discovery baseline")
    prior = sources[0]
    base.require(prior.get("effectiveCuts") == old, "100k default cuts differ from frozen discovery default")
    prior_policy = {k: v for k, v in prior.get("eventPolicy", {}).items() if k not in ignored}
    base.require(prior_policy == policy(baseline), "100k event policy differs from frozen discovery policy")
    return choice, prior


def bootstrap_sums(counts, replicates, seed, block_size=10):
    """At most block_size * nEvents Poisson weights exist at once."""
    base.require(type(block_size) is int and block_size > 0, "invalid bootstrap block size")
    choices, events, regions = counts.shape
    flattened = counts.transpose(1, 0, 2).reshape(events, choices * regions)
    result = np.empty((replicates, choices, regions), dtype=np.int64)
    rng = np.random.RandomState(seed)
    for start in range(0, replicates, block_size):
        end = min(start + block_size, replicates)
        weights = rng.poisson(1, size=(end - start, events))
        result[start:end] = weights.dot(flattened).reshape(end - start, choices, regions)
    return result


def write_bootstrap(output, datasets, c, events):
    lookup = {event: i for i, event in enumerate(events)}
    labels = ["signal"] + [sb + "_" + side for sb in SIDEBANDS for side in ("left", "right")]
    counts = np.zeros((2, len(events), len(labels)), dtype=np.int64)
    for ci, d in enumerate(datasets):
        indices = np.array([lookup[key] for key in d["event_keys"]], dtype=np.int64)
        for ri, region in enumerate(labels):
            use = d["selected"] & d["regions"][region]
            counts[ci, :, ri] = np.bincount(indices[use], minlength=len(events))
    ids = ["baseline", "best50"]
    np.savez_compressed(output / "chosen_event_counts.npz", event_ids=np.asarray(events, dtype=np.int64),
                        choice_ids=np.array(ids), region_labels=np.array(labels), counts=counts)
    sums = bootstrap_sums(counts, c["bootstrapReplicates"], c["bootstrapSeed"])
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
        for ci, name in enumerate(ids):
            for quantity, values, estimate in (("purity", purity[:, ci], p_est[ci]),
                    ("purityGain", purity[:, ci] - purity[:, 0], p_est[ci] - p_est[0]),
                    ("signalRetention", retention[:, ci], r_est[ci])):
                good = values[np.isfinite(values)]
                valid = len(good) >= math.ceil(c["bootstrapMinValidFraction"] * c["bootstrapReplicates"])
                qs = np.quantile(good, c["bootstrapQuantiles"]) if valid else [float("nan")] * 2
                rows.append(dict(choice_id=name, split="full", sideband=sb_name, quantity=quantity,
                    estimate=float(estimate), bootstrapError=float(np.std(good, ddof=1)) if valid else float("nan"),
                    quantileLow=float(qs[0]), quantileHigh=float(qs[1]), validReplicates=len(good),
                    replicates=c["bootstrapReplicates"], status="ok" if valid else "insufficient_finite_replicates"))
                draws["{}_full_{}_{}".format(name, sb_name, quantity)] = values
    base.write_csv(output / "paired_event_bootstrap.csv", list(rows[0]), rows)
    np.savez_compressed(output / "paired_event_bootstrap_draws.npz", **draws)


def run(config_path, baseline_dir, best_dir, event_ids, expected_events, frozen_path, output):
    config_path, frozen_path, output = Path(config_path), Path(frozen_path), Path(output)
    base.require(type(expected_events) is int and expected_events > 0, "expected-events must be positive")
    base.require(not output.exists(), "output directory exists; refusing overwrite")
    policy, c, base_path = scan.load_config(config_path)
    discovery_events, discovery_event_path = scan.read_events(c)
    frozen = json.loads(frozen_path.read_text())
    source_path = frozen_path.parent / "sources.json"
    source_manifest = json.loads(source_path.read_text())
    base.require(frozen.get("selectionSplit") == "full", "unexpected frozen discovery selectionSplit")
    headers = [base.read_metadata(Path(folder) / "metadata.txt")[0] for folder in (baseline_dir, best_dir)]
    base.require(all(h.get("inputEvents") == expected_events for h in headers), "wrong actual validation input-event count")
    reconstructed = headers[0].get("reconstructedEvents")
    base.require(type(reconstructed) is int and 0 < reconstructed <= expected_events and
                 headers[1].get("reconstructedEvents") == reconstructed, "default/best reconstructed counts differ or are invalid")
    c.update(inputEvents=expected_events, reconstructedEvents=reconstructed, eventIdsFile=str(event_ids),
             cuts={}, tripleFamilies=[])
    events, event_path = scan.read_events(c)
    overlap = set(events).intersection(discovery_events)
    base.require(not overlap, "validation events overlap frozen-cut discovery events")
    event_set = set(events)
    datasets = []
    for name, folder in (("baseline", baseline_dir), ("best50", best_dir)):
        d = base.load_dataset(name, Path(folder), c)
        base.require(set(d["event_keys"]).issubset(event_set), "candidate event absent from validation event ledger")
        d["regions"] = base.make_region_masks(d["data"], c)
        datasets.append(d)
    choice, discovery_baseline = audit_frozen_cuts(datasets[0], datasets[1], frozen, source_manifest)
    base.require(discovery_baseline.get("meta", {}).get("inputEvents") == len(discovery_events),
                 "discovery ledger count disagrees with frozen baseline provenance")
    references = {sb: base.metrics(datasets[0], datasets[0]["selected"], "full", sb, c) for sb in SIDEBANDS}
    base.require(all(r["Nwin"] > 0 and r["S"] > 0 for r in references.values()), "nonpositive validation baseline signal")
    output.mkdir(parents=True, exist_ok=False)
    base.write_json(output / "status.json", dict(complete=False, phase="frozen_validation"))
    (output / "configuration.yaml").write_text(config_path.read_text())
    (output / "resolved_validation_configuration.yaml").write_text(yaml.safe_dump(c, sort_keys=False))
    selections, histograms, native, metrics, metric_map = [], [], [], [], {}
    bins = int(round((c["massMax"] - c["massMin"]) / c["massBinWidth"]))
    for d in datasets:
        name = d["id"]
        selections.append(dict(choice_id=name, dataset_id=name, point_id=0,
            label="KF Imp5 default" if name == "baseline" else "Frozen best50 (independent validation)",
            role="baseline" if name == "baseline" else "frozen_validation", retention_floors="",
            cuts_json=json.dumps({} if name == "baseline" else choice["combinedCutOverrides"], sort_keys=True)))
        counts = np.bincount(d["bins"][d["selected"]], minlength=bins + 2)
        norm = int(counts[1:-1].sum())
        for i in range(bins + 2):
            low = c["massMin"] + (i - 1) * c["massBinWidth"]
            row = dict(selection=name, bin=i, low=low, high=low + c["massBinWidth"],
                       count=int(counts[i]), error=math.sqrt(counts[i]))
            native.append(row)
            if 1 <= i <= bins:
                row = dict(row); row.pop("selection")
                row.update(choice_id=name, dataset_id=name, point_id=0, split="full", normalization_count=norm)
                histograms.append(row)
        metric_map[name] = {}
        for sb in SIDEBANDS:
            value = base.add_reference(base.metrics(d, d["selected"], "full", sb, c), references[sb])
            upper = scan.refinement.poisson_upper(int(value["Nleft"] + value["Nright"]),
                                                   policy["poissonUpperConfidenceLevel"])
            value.update(poissonUpperMean=upper, backgroundUpper=value["alpha"] * upper)
            value["conservativePurity"] = 1 - value["backgroundUpper"] / value["Nwin"] if value["Nwin"] else float("nan")
            metric_map[name][sb] = value
            metrics.append(dict(choice_id=name, dataset_id=name, split="full", sideband=sb, **value))
    base.write_csv(output / "chosen_selections.csv", list(selections[0]), selections)
    base.write_csv(output / "chosen_metrics.csv", ["choice_id", "dataset_id", "split", "sideband"] + base.METRIC_FIELDS + ["poissonUpperMean", "backgroundUpper", "conservativePurity"], metrics)
    base.write_csv(output / "histograms.csv", ["choice_id", "dataset_id", "point_id", "split", "bin", "low", "high", "count", "error", "normalization_count"], histograms)
    base.write_csv(output / "chosen_histograms.csv", ["selection", "bin", "low", "high", "count", "error"], native)
    write_bootstrap(output, datasets, c, events)
    # Keep beam energy distinct from sqrt(s_NN); remove exploration wording
    # and long discovery titles from the independent frozen-cut spectrum.
    energy_token = c["title"].split()[0].lower()
    beam_energy = {"auau13p5": "13.5", "auau3p85": "3.85"}.get(energy_token)
    if beam_energy:
        title = "Au+Au FXT; {} GeV beam; {:,} input events; frozen cuts".format(beam_energy, expected_events)
    else:
        title = "KF Lambda; {:,} independent input events; frozen cuts".format(expected_events)
    manifest = {key: c[key] for key in ("massMin", "massMax", "massBinWidth", "displayMin", "displayMax")}
    manifest.update(title=title, windowLow=c["window"][0], windowHigh=c["window"][1],
                    analysisRole="independent_frozen_validation", inputEvents=expected_events, **policy["backgroundCheck"])
    manifest["backgroundRegions"] = ";".join("{}:{}:{}:{}:{}".format(sb, *c["sidebands"][sb]["left"], *c["sidebands"][sb]["right"]) for sb in SIDEBANDS)
    for filename in ("plot_manifest.tsv", "scan_config.tsv"):
        with (output / filename).open("w") as stream:
            for key, value in manifest.items(): stream.write("{}\t{}\n".format(key, value))
    def record(path): return dict(path=str(Path(path).resolve()), sha256=base.sha256(path))
    provenance = dict(config=record(config_path), frozenChoice=record(frozen_path), discoverySources=record(source_path),
        baseConfig=record(base_path), eventIdsFile=record(event_path), discoveryEventIdsFile=record(discovery_event_path),
        code=[record(__file__), record(base.__file__), record(scan.__file__)], sharedDiscoveryEvents=len(overlap),
        datasets=[dict(id=d["id"], meta=d["meta"], effectiveCuts=d["cuts"], eventPolicy=d["event_meta"],
                       inputs={k: record(p) for k, p in d["paths"].items()}) for d in datasets])
    base.write_json(output / "sources.json", provenance)
    previous = [x.get("fullMetrics", {}) for x in frozen.get("decisions", []) if x.get("choice_id") == "best50" and x.get("criterion") == "robust"]
    comparison = dict(discoveryInputEvents=len(discovery_events), validationInputEvents=expected_events,
                      discoveryBest50=previous[0] if len(previous) == 1 else None,
                      validationBest50=metric_map["best50"]["default"], role="descriptive_only_no_retuning")
    if len(previous) == 1 and previous[0].get("S", 0) > 0:
        comparison["signalPerInputEventRatio_validation_over_discovery"] = (
            metric_map["best50"]["default"]["S"] / expected_events) / (previous[0]["S"] / len(discovery_events))
    base.write_json(output / "discovery_validation_comparison.json", comparison)
    summary = dict(phase="independent_frozen_validation", inputEvents=expected_events,
        reconstructedEvents=reconstructed, sharedDiscoveryEvents=0, cutReselection=False,
        frozenChoice="best50", combinedCutOverrides=choice["combinedCutOverrides"], metrics=metric_map,
        bootstrapReplicates=c["bootstrapReplicates"], bootstrapBlockSize=10,
        baselineExactSelectedAndNativeHistogramReplay=True, bestExactSelectedAndNativeHistogramReplay=True,
        limitations="Both spectra are actual completed reconstruction selections, not post-hoc candidate masks. Purity uses the fixed discovery mass window and sidebands; model-dependent backgrounds and feeddown remain. Paired event bootstrap uses frozen cuts and does not estimate detector efficiency. Discovery/validation yield-per-event differences can reflect run/event populations; do not retune on this sample.",
        plotting="CSV schema supports checkLambdaPurityRefinement.C numeric background checks, but its hardcoded exploratory/previously-inspected annotations are inappropriate for this independent sample. Use a correctly labeled final overlay.")
    base.write_json(output / "summary.json", summary)
    base.write_json(output / "status.json", dict(complete=True, phase="complete"))
    products = sorted(p for p in output.iterdir() if p.is_file())
    base.require(all(p.stat().st_size > 0 for p in products), "empty frozen-validation output")
    with (output / "SHA256SUMS").open("w") as stream:
        for path in products: stream.write("{}  {}\n".format(base.sha256(path), path.name))
    print("PASS frozen comparison: {:,} independent input events; no cut reselection".format(expected_events))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--best", required=True, type=Path)
    parser.add_argument("--event-ids", required=True, type=Path)
    parser.add_argument("--expected-events", required=True, type=int)
    parser.add_argument("--frozen-choice", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        run(args.config, args.baseline, args.best, args.event_ids, args.expected_events, args.frozen_choice, args.output)
    except (ValueError, OSError, KeyError, csv.Error, yaml.YAMLError) as error:
        print("ERROR: " + str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
