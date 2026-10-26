#!/usr/bin/env python3
"""Independent, stdlib-only audit of frozen-comparison numeric products.

No imports from the production scan/comparison helpers are used by the audit.
Actual selected-candidate exports and all native histogram bins are replayed.
Retention below 50% on validation is reported, NEVER rejected or reoptimized.

--comparison DIR --frozen-manifest FROZEN_DECISION.json --output NEW_AUDIT.json
[--expected-events 100000]; or --self-test for synthetic positive/negative tests.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys


def require(value, message):
    if not value: raise ValueError(message)


def number(value):
    result = float(value)
    require(math.isfinite(result), "nonfinite numeric audit input")
    return result


def integer(value):
    result = number(value)
    require(result == int(result), "noninteger count or event identifier")
    return int(result)


def close(actual, expected, name):
    require(math.isclose(number(actual), number(expected), rel_tol=1e-10, abs_tol=1e-9),
            "numeric mismatch {}: {} != {}".format(name, actual, expected))


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for data in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(data)
    return digest.hexdigest()


def check_record(record):
    path = Path(record["path"])
    require(sha256(path) == record["sha256"], "source hash changed: " + str(path))
    return path


def ledger(path, expected=None):
    result = []
    with Path(path).open(newline="") as stream:
        for line, fields in enumerate(csv.reader(stream, delimiter="\t")):
            if line == 0 and fields == ["index", "runId", "eventId"]: continue
            require(len(fields) == 3, "invalid ledger column count")
            index, run, event = map(int, fields)
            require(index == len(result), "ledger index discontinuity")
            result.append((run, event))
    require(result and len(result) == len(set(result)), "empty/duplicate input events")
    if expected is not None: require(len(result) == expected, "ledger event count mismatch")
    return set(result), len(result)


def metadata_counts(path):
    result = {}
    with Path(path).open() as stream:
        for line in stream:
            if line.strip() == "KFParticleEffectiveConfiguration": break
            if ":" not in line: continue
            key, value = line.split(":", 1)
            if key in ("inputEvents", "reconstructedEvents", "rawEntries"):
                result[key] = int(value.strip())
    require(set(result) == {"inputEvents", "reconstructedEvents", "rawEntries"}, "incomplete export metadata")
    return result


def manifest(path):
    result = {}
    with Path(path).open() as stream:
        for line in stream:
            if not line.strip() or line.startswith("#"): continue
            key, value = line.rstrip("\n\r").split("\t", 1)
            require(key not in result, "duplicate plot manifest key")
            result[key] = value
    require(result["analysisRole"] == "independent_frozen_validation", "wrong manifest analysis role")
    return result


def sidebands(plot):
    result = {}
    for group in plot["backgroundRegions"].split(";"):
        fields = group.split(":")
        require(len(fields) == 5 and fields[0] not in result, "invalid background region")
        result[fields[0]] = list(map(number, fields[1:]))
    require(set(result) == {"default", "near", "far"}, "wrong sideband set")
    return result


def inspect_export(source, event_set, expected_events, plot, bands):
    inputs = {key: check_record(value) for key, value in source["inputs"].items()}
    metadata = metadata_counts(inputs["metadata"])
    require(metadata["inputEvents"] == expected_events, "actual export input-event count mismatch")
    require(0 < metadata["reconstructedEvents"] <= expected_events, "invalid reconstructed-event count")
    for key, value in metadata.items(): require(source["meta"][key] == value, "source/export metadata mismatch")
    low, high, width = [number(plot[k]) for k in ("massMin", "massMax", "massBinWidth")]
    bins = int(round((high - low) / width)); close(low + bins * width, high, "axis width")
    counts = [0] * (bins + 2)
    wl, wh = number(plot["windowLow"]), number(plot["windowHigh"])
    nwin, nraw, nselected = 0, 0, 0
    counts_sb = {name: [0, 0] for name in bands}
    with inputs["candidates"].open(newline="") as stream:
        reader = csv.DictReader(stream)
        required = {"runId", "eventId", "selected", "mass", "pdg"}
        require(required.issubset(reader.fieldnames or []), "missing candidate audit columns")
        for row in reader:
            nraw += 1
            event = integer(row["runId"]), integer(row["eventId"])
            require(event in event_set, "candidate event absent from input ledger")
            require(integer(row["pdg"]) == 3122, "unexpected non-Lambda candidate")
            selected = integer(row["selected"])
            require(selected in (0, 1), "nonboolean selected flag")
            if not selected: continue
            nselected += 1
            mass = number(row["mass"])
            index = 0 if mass < low else (bins + 1 if mass >= high else int(math.floor(bins * ((mass - low) / (high - low)))) + 1)
            require(0 <= index <= bins + 1, "invalid candidate mass bin")
            counts[index] += 1
            nwin += int(wl <= mass < wh)
            for name, (ll, lh, rl, rh) in bands.items():
                counts_sb[name][0] += int(ll <= mass < lh)
                counts_sb[name][1] += int(rl <= mass < rh)
    require(nraw == metadata["rawEntries"] and nselected > 0, "raw count mismatch or empty selected sample")
    with inputs["histogram"].open(newline="") as stream:
        native = list(csv.DictReader(stream))
    require(len(native) == bins + 2, "native histogram missing bins/flow cells")
    for index, row in enumerate(native):
        require(integer(row["bin"]) == index, "native bin order mismatch")
        require(integer(row["count"]) == counts[index], "native histogram differs from actual selected candidates")
        close(row["low"], low + (index - 1) * width, "native low edge")
        close(row["high"], low + index * width, "native high edge")
        close(number(row["error"]) ** 2, counts[index], "native Poisson error")
    measured = {}
    for name, (left, right) in counts_sb.items():
        ll, lh, rl, rh = bands[name]
        alpha = (wh - wl) / (lh - ll + rh - rl)
        background = alpha * (left + right)
        signal = nwin - background
        measured[name] = dict(Nwin=nwin, Nleft=left, Nright=right, alpha=alpha, B=background, S=signal,
            BError=alpha * math.sqrt(left + right), SError=math.sqrt(nwin + alpha * alpha * (left + right)),
            purity=signal / nwin if nwin else None,
            purityError=math.sqrt(alpha * alpha * (left + right) / (nwin * nwin) + background * background / (nwin ** 3)) if nwin else None)
    return dict(metadata=metadata, counts=counts, metrics=measured, normalizer=sum(counts[1:-1]), bins=bins)


def audit(comparison, frozen_manifest, expected_events=100000):
    comparison, frozen_manifest = Path(comparison), Path(frozen_manifest)
    require(type(expected_events) is int and expected_events > 0, "invalid expected event count")
    status = json.loads((comparison / "status.json").read_text())
    require(status.get("complete") is True, "comparison incomplete")
    summary = json.loads((comparison / "summary.json").read_text())
    require(summary["phase"] == "independent_frozen_validation" and summary["cutReselection"] is False,
            "comparison is not a frozen validation")
    require(summary["inputEvents"] == expected_events, "summary event count mismatch")
    sources = json.loads((comparison / "sources.json").read_text())
    for name in ("config", "baseConfig", "discoverySources"):
        check_record(sources[name])
    for record in sources["code"]: check_record(record)
    events, _ = ledger(check_record(sources["eventIdsFile"]), expected_events)
    discovery, discovery_count = ledger(check_record(sources["discoveryEventIdsFile"]))
    require(not events.intersection(discovery), "discovery/validation event overlap")
    require(summary["sharedDiscoveryEvents"] == sources["sharedDiscoveryEvents"] == 0, "overlap status mismatch")
    decision_path = check_record(sources["frozenChoice"])
    freeze = json.loads(frozen_manifest.read_text())
    require(freeze["noRetuningOnValidation"] is True and freeze["validationExaminedBeforeFreeze"] is False,
            "freeze manifest does not assert pre-validation fixed cuts")
    require(freeze["selectionEvents"] == discovery_count and freeze["validationEvents"] == expected_events,
            "freeze manifest event counts differ")
    require(freeze["decisionSha256"] == sha256(decision_path), "initial frozen-decision SHA does not match")
    require(Path(freeze["decisionPath"]).resolve() == decision_path.resolve(), "wrong initial frozen-decision path")
    decision = json.loads(decision_path.read_text())
    chosen = [x for x in decision["choices"] if x["choice_id"] == "best50"]
    require(len(chosen) == 1 and chosen[0] == freeze["choice"], "initial frozen best50 content changed")
    require(chosen[0]["combinedCutOverrides"] == summary["combinedCutOverrides"], "summary frozen cuts changed")
    plot = manifest(comparison / "plot_manifest.tsv")
    require(integer(plot["inputEvents"]) == expected_events, "plot event count mismatch")
    bands = sidebands(plot)
    actual = {}
    require(len(sources["datasets"]) == 2, "need exactly default and best datasets")
    for source in sources["datasets"]:
        name = source["id"]
        require(name in ("baseline", "best50") and name not in actual, "unknown/duplicate source dataset")
        actual[name] = inspect_export(source, events, expected_events, plot, bands)
    require(actual["baseline"]["metadata"]["reconstructedEvents"] == actual["best50"]["metadata"]["reconstructedEvents"]
            == summary["reconstructedEvents"], "reconstruction event counts differ")
    with (comparison / "chosen_metrics.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    require(len(rows) == 6, "expected two selections times three sideband metrics")
    seen = set()
    for row in rows:
        key = row["choice_id"], row["sideband"]
        require(key not in seen and key[0] in actual and key[1] in bands and row["split"] == "full", "duplicate/unknown metric")
        seen.add(key)
        measured = actual[key[0]]["metrics"][key[1]]
        baseline = actual["baseline"]["metrics"][key[1]]
        require(baseline["S"] > 0 and measured["Nwin"] > 0, "undefined validation reference metric")
        measured.update(signalRetention=measured["S"] / baseline["S"], NwinRetention=measured["Nwin"] / baseline["Nwin"],
                        baselineS=baseline["S"], baselineNwin=baseline["Nwin"])
        for field, expected in measured.items():
            close(row[field], expected, "CSV " + str(key) + " " + field)
            close(summary["metrics"][key[0]][key[1]][field], expected, "summary " + str(key) + " " + field)
    with (comparison / "histograms.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    nbins = actual["baseline"]["bins"]
    require(len(rows) == 2 * nbins, "comparison histogram missing bins")
    seen = set()
    low, width = number(plot["massMin"]), number(plot["massBinWidth"])
    for row in rows:
        key = row["choice_id"], integer(row["bin"])
        require(key not in seen and key[0] in actual and 1 <= key[1] <= nbins and row["split"] == "full", "invalid histogram bin")
        seen.add(key)
        count = actual[key[0]]["counts"][key[1]]
        close(row["count"], count, "comparison histogram count")
        close(number(row["error"]) ** 2, count, "comparison Poisson error")
        close(row["low"], low + (key[1] - 1) * width, "comparison low edge")
        close(row["high"], low + key[1] * width, "comparison high edge")
        close(row["normalization_count"], actual[key[0]]["normalizer"], "whole-range normalization count")
    for line in (comparison / "SHA256SUMS").read_text().splitlines():
        expected, name = line.split("  ", 1)
        require(sha256(comparison / name) == expected, "comparison product checksum changed: " + name)
    retentions = {sb: actual["best50"]["metrics"][sb]["signalRetention"] for sb in bands}
    return dict(passed=True, expectedInputEvents=expected_events, reconstructedEvents=summary["reconstructedEvents"],
        discoveryInputEvents=discovery_count, sharedEvents=0, initialFrozenDecisionSha256=freeze["decisionSha256"],
        actualSelectedCandidateAndNativeHistogramReplay=True, allNumericMetricsRecomputedIndependently=True,
        fullRangeNormalizationAudited=True, minimumValidationSignalRetention=min(retentions.values()),
        validationRetentionAtLeast50PercentAcrossSidebands=all(x >= .5 for x in retentions.values()),
        validationSignalRetentionBySideband=retentions,
        retentionPolicy="Observation only: a validation yield below 50% is not a QA failure or permission to retune.",
        metrics={name: value["metrics"] for name, value in actual.items()})


def self_test():
    import importlib.util
    import tempfile
    import unittest
    repo = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location("frozen_audit_fixture", repo / "tests/test_lambda_kf_frozen_purity.py")
    fixture = importlib.util.module_from_spec(spec); spec.loader.exec_module(fixture)

    class AuditTest(unittest.TestCase):
        def make(self, folder):
            config, baseline, best, events, choice = fixture.fixture(folder)
            comparison = folder / "comparison"
            fixture.frozen.run(config, baseline, best, events, 100, choice, comparison)
            chosen = json.loads(choice.read_text())["choices"][0]
            freeze = dict(noRetuningOnValidation=True, validationExaminedBeforeFreeze=False,
                selectionEvents=100, validationEvents=100, decisionPath=str(choice), decisionSha256=sha256(choice), choice=chosen)
            freeze_path = folder / "frozen_decision.json"
            freeze_path.write_text(json.dumps(freeze))
            return comparison, freeze_path

        def test_numeric_replay_and_frozen_hash(self):
            with tempfile.TemporaryDirectory(prefix="lambda-final-audit-") as tmp:
                comparison, frozen = self.make(Path(tmp))
                result = audit(comparison, frozen, 100)
                self.assertTrue(result["passed"])
                self.assertAlmostEqual(result["metrics"]["best50"]["default"]["signalRetention"], 360 / 585)
                obj = json.loads(frozen.read_text()); obj["decisionSha256"] = "0" * 64
                frozen.write_text(json.dumps(obj))
                with self.assertRaisesRegex(ValueError, "frozen-decision SHA"):
                    audit(comparison, frozen, 100)

        def test_tampered_numeric_summary_rejected(self):
            with tempfile.TemporaryDirectory(prefix="lambda-final-audit-") as tmp:
                comparison, frozen = self.make(Path(tmp))
                path = comparison / "summary.json"
                obj = json.loads(path.read_text()); obj["metrics"]["best50"]["default"]["S"] += 1
                path.write_text(json.dumps(obj))
                with self.assertRaisesRegex(ValueError, "numeric mismatch"):
                    audit(comparison, frozen, 100)

    return unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(AuditTest)).wasSuccessful()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", type=Path)
    parser.add_argument("--frozen-manifest", type=Path)
    parser.add_argument("--expected-events", type=int, default=100000)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test: return 0 if self_test() else 1
    if not (args.comparison and args.frozen_manifest and args.output):
        parser.error("--comparison, --frozen-manifest and --output are required")
    try:
        require(not args.output.exists(), "refusing existing audit output")
        result = audit(args.comparison, args.frozen_manifest, args.expected_events)
        with args.output.open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False); stream.write("\n")
        print("PASS independent numeric audit; minimum validation S retention = {:.4%}".format(result["minimumValidationSignalRetention"]))
    except (ValueError, KeyError, OSError, csv.Error) as error:
        print("ERROR: " + str(error), file=sys.stderr); return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
