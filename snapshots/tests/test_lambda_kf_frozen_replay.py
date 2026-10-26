#!/usr/bin/env python3
"""Synthetic, no-ROOT verification of frozen upstream/downstream replay."""
import copy
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import yaml

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "common/macro"))
import verifyLambdaKfFrozenReplay as verify

SPEC = importlib.util.spec_from_file_location("replay_fixture", REPO / "tests/test_lambda_kf_upstream_purity.py")
FIXTURE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(FIXTURE)


def rewrite_export(folder, table, cuts=None, event=None):
    head, old_cuts, old_event = verify.base.read_metadata(folder / "metadata.txt")
    cuts = old_cuts if cuts is None else cuts
    event = old_event if event is None else event
    rows = [dict(zip(table, values)) for values in zip(*table.values())]
    verify.base.write_csv(folder / "candidates.csv", list(table), rows)
    head["rawEntries"] = len(rows)
    FIXTURE.write_metadata(folder / "metadata.txt", head, cuts, event)
    counts = np.bincount(verify.base.mass_bins(table["mass"][table["selected"] == 1],
                        {"massMin": 1.05, "massMax": 1.25, "massBinWidth": .001}), minlength=202)
    hist = [dict(bin=i, low=1.05 + (i - 1) * .001, high=1.05 + i * .001,
                 count=int(counts[i]), error=float(np.sqrt(counts[i]))) for i in range(202)]
    verify.base.write_csv(folder / "selected_histogram.csv", list(hist[0]), hist)


def fixture(folder, optional_cut=False):
    config, datasets, export, source = FIXTURE.fixture(folder)
    scan = folder / "scan"
    verify.upstream.run(config, datasets, scan)
    chosen = json.loads((scan / "chosen_points.json").read_text())
    best = next(x for x in chosen["choices"] if x["choice_id"] == "best50")
    _, c, _ = verify.upstream.load_config(config)
    if optional_cut:
        point = {"minDecayLengthSignificance": 10.}
        best["point_id"] = next(i for i, p in enumerate(verify.base.enumerate_points(c)) if dict(p) == point)
        best["downstreamCuts"] = point
        best["combinedCutOverrides"].update(point)
        (scan / "chosen_points.json").write_text(json.dumps(chosen))
    head, cuts, event = verify.base.read_metadata(source / "metadata.txt")
    cuts.update(best["downstreamCuts"])
    replay = folder / "replay"; replay.mkdir()
    FIXTURE.write_metadata(replay / "metadata.txt", dict(head, inputRoot="synthetic_replay.root"), cuts, event)
    table = verify.base.read_numeric_csv(source / "candidates.csv")
    table["selected"] = verify.base.replay_selected(table, cuts).astype(float)
    rewrite_export(replay, table)
    return config, scan, replay, source


class FrozenReplayTest(unittest.TestCase):
    def test_exact_replay_and_no_overwrite(self):
        with tempfile.TemporaryDirectory(prefix="lambda-frozen-replay-") as tmp:
            folder = Path(tmp)
            config, scan, replay, source = fixture(folder)
            output = folder / "audit.json"
            report = verify.run(scan, replay, config, output)
            self.assertTrue(report["passed"], report)
            self.assertEqual(report["expectedSelected"], 360)
            self.assertEqual(report["nativeBinsIncludingFlows"], 202)
            self.assertEqual(report["maxAbsMassDifference"], 0.)
            before = output.read_bytes()
            with self.assertRaisesRegex(ValueError, "refusing overwrite"):
                verify.run(scan, replay, config, output)
            self.assertEqual(output.read_bytes(), before)

    def test_candidate_order_irrelevant_and_optional_final_cut_supported(self):
        with tempfile.TemporaryDirectory(prefix="lambda-frozen-replay-") as tmp:
            folder = Path(tmp)
            config, scan, replay, source = fixture(folder, optional_cut=True)
            table = verify.base.read_numeric_csv(replay / "candidates.csv")
            rewrite_export(replay, {key: values[::-1] for key, values in table.items()})
            report = verify.run(scan, replay, config, folder / "audit.json")
            self.assertTrue(report["passed"], report)
            self.assertEqual(report["actualEffectiveCuts"]["minDecayLengthSignificance"], 10.)

    def test_pair_identity_mismatch_fails_even_with_identical_histogram(self):
        with tempfile.TemporaryDirectory(prefix="lambda-frozen-replay-") as tmp:
            folder = Path(tmp)
            config, scan, replay, source = fixture(folder)
            table = verify.base.read_numeric_csv(replay / "candidates.csv")
            table["protonId"][0] = 99
            rewrite_export(replay, table)
            report = verify.run(scan, replay, config, folder / "audit.json")
            self.assertFalse(report["passed"])
            self.assertFalse(report["checks"]["selectedPairMultisetsMatch"])
            self.assertTrue(report["checks"]["allNativeMassBinsMatch"])
            self.assertEqual(len(report["missingPairExamples"]), 1)
            self.assertIsNone(report["maxAbsMassDifference"])

    def test_pair_multiplicity_is_not_collapsed_to_set(self):
        with tempfile.TemporaryDirectory(prefix="lambda-frozen-replay-") as tmp:
            folder = Path(tmp)
            config, scan, replay, source = fixture(folder)
            table = verify.base.read_numeric_csv(replay / "candidates.csv")
            selected_indices = np.flatnonzero(table["selected"] == 1)
            a, b = selected_indices[:2]
            for key in verify.PAIR_FIELDS: table[key][b] = table[key][a]
            rewrite_export(replay, table)
            report = verify.run(scan, replay, config, folder / "audit.json")
            self.assertFalse(report["passed"])
            self.assertEqual(report["actualSelected"], report["expectedSelected"])
            self.assertEqual(report["actualUniquePairs"], report["expectedUniquePairs"] - 1)

    def test_mass_bin_difference_fails_without_tolerance(self):
        with tempfile.TemporaryDirectory(prefix="lambda-frozen-replay-") as tmp:
            folder = Path(tmp)
            config, scan, replay, source = fixture(folder)
            table = verify.base.read_numeric_csv(replay / "candidates.csv")
            table["mass"][0] += .0011
            rewrite_export(replay, table)
            report = verify.run(scan, replay, config, folder / "audit.json")
            self.assertFalse(report["passed"])
            self.assertTrue(report["checks"]["selectedPairMultisetsMatch"])
            self.assertEqual(len(report["mismatchedBins"]), 2)
            self.assertAlmostEqual(report["maxAbsMassDifference"], .0011)

    def test_sub_bin_mass_change_reported_not_hidden(self):
        with tempfile.TemporaryDirectory(prefix="lambda-frozen-replay-") as tmp:
            folder = Path(tmp)
            config, scan, replay, source = fixture(folder)
            table = verify.base.read_numeric_csv(replay / "candidates.csv")
            table["mass"][0] += 1.e-7
            rewrite_export(replay, table)
            report = verify.run(scan, replay, config, folder / "audit.json")
            self.assertTrue(report["passed"], report)
            self.assertEqual(report["unequalMatchedMassValues"], 1)
            self.assertAlmostEqual(report["maxAbsMassDifference"], 1.e-7)

    def test_effective_cut_and_event_policy_mismatches_fail(self):
        for kind in ("cut", "policy"):
            with tempfile.TemporaryDirectory(prefix="lambda-frozen-replay-") as tmp:
                folder = Path(tmp)
                config, scan, replay, source = fixture(folder)
                table = verify.base.read_numeric_csv(replay / "candidates.csv")
                head, cuts, event = verify.base.read_metadata(replay / "metadata.txt")
                if kind == "cut": cuts["minNHitsFit"] += 1
                else: event["KFEventSelection.policy"] = "changed_policy"
                rewrite_export(replay, table, cuts, event)
                report = verify.run(scan, replay, config, folder / "audit.json")
                self.assertFalse(report["passed"])
                self.assertIn("effective " + ("cuts" if kind == "cut" else "event policy"), report["error"])

    def test_source_hash_tampering_and_event_ledger_changes_fail(self):
        for kind in ("source", "ledger"):
            with tempfile.TemporaryDirectory(prefix="lambda-frozen-replay-") as tmp:
                folder = Path(tmp)
                config, scan, replay, source = fixture(folder)
                if kind == "source":
                    target = source / "metadata.txt"
                    target.write_text(target.read_text() + "\n")
                else:
                    target = folder / "events.tsv"
                    target.write_text(target.read_text().replace("0\t1\t1\n", "0\t1\t9999\n"))
                report = verify.run(scan, replay, config, folder / "audit.json")
                self.assertFalse(report["passed"])
                self.assertIn("changed", report["error"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
