#!/usr/bin/env python3
"""Verify a frozen best50 PicoDst rerun against its exploratory selection.

Read-only CSV/metadata verification; no ROOT or PicoDst event loop. Candidate
order is deliberately ignored, but daughter-pair multiplicities and all native
mass bins (including both flows) must agree exactly. Sub-bin mass differences
are reported without silently applying a numerical tolerance to bin counts.
"""
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import sys

import numpy as np
import yaml
import scanLambdaKfUpstreamPurity as upstream
import studyLambdaKfPurity as base

PAIR_FIELDS = ("runId", "eventId", "protonId", "pionId", "pdg")
EVENT_PATH_KEYS = {"KFEventSelection.mainconf", "KFEventSelection.analysisInfo",
                   "KFEventSelection.eventConfig"}


def record(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": base.sha256(path)}


def event_policy(dataset):
    return {k: v for k, v in dataset["event_meta"].items() if k not in EVENT_PATH_KEYS}


def load_recorded_dataset(source, c):
    """Use exactly the export paths and hashes retained by the scanner."""
    names = {"candidates": "candidates.csv", "histogram": "selected_histogram.csv",
             "metadata": "metadata.txt"}
    base.require(set(source["inputs"]) == set(names), "unexpected recorded export files")
    folder = Path(source["inputs"]["metadata"]["path"]).resolve().parent
    for key, name in names.items():
        item = source["inputs"][key]
        path = Path(item["path"]).resolve()
        base.require(path == folder / name, "recorded export files have inconsistent locations")
        base.require(base.sha256(path) == item["sha256"], "source export hash changed: " + str(path))
    dataset = base.load_dataset(source["id"], folder, c)
    base.require(dataset["cuts"] == source["effectiveCuts"], "source effective metadata changed")
    base.require(dataset["event_meta"] == source["eventPolicy"], "source event metadata changed")
    return dataset


def selected_pairs(dataset, mask):
    values = dataset["data"]
    return [tuple(int(values[key][index]) for key in PAIR_FIELDS)
            for index in np.flatnonzero(mask)]


def pair_mass_groups(dataset, mask):
    result = defaultdict(list)
    indices = np.flatnonzero(mask)
    for key, index in zip(selected_pairs(dataset, mask), indices):
        result[key].append(float(dataset["data"]["mass"][index]))
    return {key: sorted(values) for key, values in result.items()}


def differences(a, b):
    delta = a - b
    return [{"pair": list(pair), "multiplicity": int(count)}
            for pair, count in sorted(delta.items())[:20]]


def verify(scan, replay_export, config):
    scan, replay_export, config = Path(scan).resolve(), Path(replay_export).resolve(), Path(config).resolve()
    status = json.loads((scan / "status.json").read_text())
    base.require(status.get("complete") is True, "exploratory scan is incomplete")
    choices = json.loads((scan / "chosen_points.json").read_text())
    sources = json.loads((scan / "sources.json").read_text())
    base.require(choices.get("selectionSplit") == "full", "unexpected selection split")
    best = [x for x in choices["choices"] if x["choice_id"] == "best50"]
    base.require(len(best) == 1, "need exactly one frozen best50 choice")
    best = best[0]
    policy, c, _ = upstream.load_config(config)
    base.require(base.sha256(config) == sources["policy"]["sha256"], "score configuration changed after ranking")
    events, event_path = upstream.read_events(c)
    base.require(base.sha256(event_path) == sources["eventIdsFile"]["sha256"], "input event ledger changed")
    event_set = set(events)
    source_by_id = {item["id"]: item for item in sources["datasets"]}
    base.require(len(source_by_id) == len(sources["datasets"]), "duplicate source dataset IDs")
    baseline_id = choices["baselineDataset"]
    base.require(sources["datasets"][0]["id"] == baseline_id, "baseline identity mismatch")
    baseline = load_recorded_dataset(source_by_id[baseline_id], c)
    selected_source = best["dataset_id"]
    original = baseline if selected_source == baseline_id else load_recorded_dataset(source_by_id[selected_source], c)
    changes = upstream.audit_source_differences(original, baseline)
    points = base.enumerate_points(c)
    point_id = best["point_id"]
    base.require(type(point_id) is int and 0 <= point_id < len(points), "invalid selected point ID")
    point = points[point_id]
    base.require(dict(point) == best["downstreamCuts"], "selected point differs from downstreamCuts")
    overrides = {key: item["value"] for key, item in changes.items()}
    overrides.update(dict(point))
    base.require(overrides == best["combinedCutOverrides"], "combined frozen overrides do not reproduce selected point")
    expected_cuts = dict(baseline["cuts"], **overrides)
    original["cut_masks"] = base.cut_masks(original, c)
    expected_mask = base.point_mask(original, point)
    # The replay is already tightened; original grid levels may now be looser
    # than its bounds. Keep full source/selected/histogram validation, but do
    # not try to run the old exploratory grid on this finalized output.
    replay_config = dict(c, cuts={})
    actual = base.load_dataset("best50_replay", replay_export, replay_config)
    base.require(actual["cuts"] == expected_cuts, "replay effective cuts differ from baseline plus frozen overrides")
    base.require(event_policy(actual) == event_policy(baseline), "replay effective event policy changed")
    for dataset in (baseline, original, actual):
        base.require(set(dataset["event_keys"]).issubset(event_set), "candidate event absent from frozen input ledger")
    actual_mask = actual["selected"]
    expected_pairs = Counter(selected_pairs(original, expected_mask))
    actual_pairs = Counter(selected_pairs(actual, actual_mask))
    pair_match = expected_pairs == actual_pairs
    nbins = int(round((c["massMax"] - c["massMin"]) / c["massBinWidth"]))
    expected_counts = np.bincount(original["bins"][expected_mask], minlength=nbins + 2)
    actual_counts = np.bincount(actual["bins"][actual_mask], minlength=nbins + 2)
    # Verify the selected point against the scanner's independently saved
    # all-bin table too, not only a freshly recomputed mask.
    with (scan / "chosen_histograms.csv").open(newline="") as stream:
        frozen_rows = [row for row in csv.DictReader(stream) if row["selection"] == "best50"]
    base.require([int(row["bin"]) for row in frozen_rows] == list(range(nbins + 2)), "frozen histogram lacks ordered native bins")
    frozen_counts = np.asarray([float(row["count"]) for row in frozen_rows])
    base.require(np.array_equal(expected_counts, frozen_counts), "selected source mask differs from frozen histogram")
    mismatch_bins = np.flatnonzero(expected_counts != actual_counts)
    mass_difference = None
    changed_masses = None
    if pair_match:
        expected_groups = pair_mass_groups(original, expected_mask)
        actual_groups = pair_mass_groups(actual, actual_mask)
        delta = [abs(a - b) for key in expected_groups
                 for a, b in zip(expected_groups[key], actual_groups[key])]
        mass_difference = max(delta, default=0.)
        changed_masses = sum(x != 0. for x in delta)
    passed = pair_match and len(mismatch_bins) == 0
    return dict(schemaVersion=1, passed=bool(passed), choice_id="best50", sourceDataset=selected_source,
        checks=dict(sourceProvenanceHashes=True, sourceSelectedBitsAndNativeHistogramReplay=True,
            actualSelectedBitsAndNativeHistogramReplay=True, effectiveCutsMatch=True,
            effectiveEventPolicyMatch=True, eventCountsMatch=True, candidateEventsInFrozenLedger=True,
            frozenHistogramMatchesSourceMask=True, selectedPairMultisetsMatch=pair_match,
            allNativeMassBinsMatch=len(mismatch_bins) == 0),
        inputEvents=c["inputEvents"], reconstructedEvents=c["reconstructedEvents"],
        pairFields=list(PAIR_FIELDS), expectedSelected=int(expected_mask.sum()), actualSelected=int(actual_mask.sum()),
        expectedUniquePairs=len(expected_pairs), actualUniquePairs=len(actual_pairs),
        missingPairExamples=differences(expected_pairs, actual_pairs), extraPairExamples=differences(actual_pairs, expected_pairs),
        nativeBinsIncludingFlows=nbins + 2, massAxis=[c["massMin"], c["massMax"]],
        mismatchedBins=[dict(bin=int(i), expected=int(expected_counts[i]), actual=int(actual_counts[i])) for i in mismatch_bins],
        maxAbsMassDifference=mass_difference, unequalMatchedMassValues=changed_masses,
        massComparison="No mass tolerance used for acceptance: exact pair multiplicities and exact bin counts required. Per-pair sorted masses report sub-bin differences only.",
        expectedEffectiveCuts=expected_cuts, actualEffectiveCuts=actual["cuts"],
        provenance=dict(scanStatus=record(scan / "status.json"), chosenPoints=record(scan / "chosen_points.json"),
            scanSources=record(scan / "sources.json"), frozenHistogram=record(scan / "chosen_histograms.csv"),
            scoreConfiguration=record(config), inputEventIds=record(event_path),
            verifier=record(__file__), replayInputs={key: record(path) for key, path in actual["paths"].items()}),
        limitations="Validates replay exports and candidate event membership; the event ledger and reconstruction run/input provenance establish zero-candidate events and ordered event identity. This is not a purity or efficiency confidence interval.")


def run(scan, replay_export, config, output):
    output = Path(output)
    base.require(not output.exists(), "output already exists; refusing overwrite")
    try:
        report = verify(scan, replay_export, config)
    except (ValueError, OSError, KeyError, TypeError, csv.Error, yaml.YAMLError) as error:
        report = dict(schemaVersion=1, passed=False, error=str(error))
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan", required=True, type=Path)
    parser.add_argument("--replay-export", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        report = run(args.scan, args.replay_export, args.config, args.output)
    except (ValueError, OSError) as error:
        print("ERROR: " + str(error), file=sys.stderr)
        return 2
    print(("PASS" if report["passed"] else "FAIL") + " frozen best50 replay: " + str(args.output))
    if "error" in report: print(report["error"], file=sys.stderr)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
