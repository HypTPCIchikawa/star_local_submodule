#!/usr/bin/env python3
"""Synthetic, no-ROOT regression checks for joint upstream purity ranking."""
import copy
import csv
import importlib.util
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import yaml

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "common/macro"))
import scanLambdaKfUpstreamPurity as scan

SPEC = importlib.util.spec_from_file_location("upstream_fixture", REPO / "tests/test_lambda_purity_scan.py")
FIXTURE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(FIXTURE)


def write_metadata(path, head, cuts, event):
    path.write_text(yaml.safe_dump(head) + "\nKFParticleEffectiveConfiguration\n" + yaml.safe_dump(cuts)
                    + "\nKFEventSelectionConfiguration\n" + yaml.safe_dump(event))


def fixture(folder):
    c, base_config, export, datasets = FIXTURE.fixture(folder)
    head, cuts, event = scan.base.read_metadata(export / "metadata.txt")
    cuts.update(minNHitsFit=15, minNHitsRatio=.52, finderLdL2D=3.)
    write_metadata(export / "metadata.txt", head, cuts, event)
    upstream = folder / "upstream"; upstream.mkdir()
    table = scan.base.read_numeric_csv(export / "candidates.csv")
    # Half the true signal, no SB candidates: independent upstream output,
    # including invalid final-cut rows for replay coverage.
    take = (table["protonId"] <= 7) | (table["selected"] == 0)
    rows = [dict(zip(table, values)) for values in zip(*(values[take] for values in table.values()))]
    scan.base.write_csv(upstream / "candidates.csv", list(table), rows)
    selected_mass = np.array([r["mass"] for r in rows if r["selected"]])
    counts = np.bincount(scan.base.mass_bins(selected_mass, c), minlength=202)
    hist = [dict(bin=i, low=1.05 + (i - 1) * .001, high=1.05 + i * .001,
                 count=int(counts[i]), error=math.sqrt(counts[i])) for i in range(202)]
    scan.base.write_csv(upstream / "selected_histogram.csv", list(hist[0]), hist)
    changed = dict(cuts, minNHitsFit=20, minNHitsRatio=.6)
    write_metadata(upstream / "metadata.txt", dict(head, rawEntries=len(rows), inputRoot="synthetic_upstream.root"), changed, event)
    datasets.write_text("id\texport_dir\nbaseline\texport\nhits\tupstream\n")
    policy = dict(schemaVersion=1, baseConfig=str(base_config), eventIdsFile=c["eventIdsFile"],
        inputEvents=100, reconstructedEvents=100, grid={"nSigmaProton": [2.5], "minDecayLengthSignificance": [10.]},
        includeAllPairs=True, tripleAnchor="minDecayLengthSignificance", retentionFloors=[.5, .7, .8],
        minimumWindowCount=20, poissonUpperConfidenceLevel=.9, maxPointsPerDataset=100,
        title="Synthetic full-data upstream scan", bootstrapReplicates=200)
    config = folder / "upstream.yaml"; config.write_text(yaml.safe_dump(policy, sort_keys=False))
    return config, datasets, export, upstream


class UpstreamPurityTest(unittest.TestCase):
    def test_sparse_background_upper_is_positive(self):
        self.assertAlmostEqual(scan.refinement.poisson_upper(0, .9), math.log(10), places=13)
        self.assertLess(1 - .75 * scan.refinement.poisson_upper(0, .9) / 20, 1.)

    def test_half_retention_and_stable_joint_tiebreak(self):
        common = dict(selectionEligible=1, robustPurityScore=.98, purity=.99, changedParameters=2, S=100.)
        rows = [dict(common, dataset_id="z", point_id=0, minSignalRetentionAcrossSidebands=.5),
                dict(common, dataset_id="a", point_id=0, minSignalRetentionAcrossSidebands=.5),
                dict(common, dataset_id="bad", point_id=0, minSignalRetentionAcrossSidebands=.499,
                     robustPurityScore=.9999)]
        self.assertEqual(scan.select_row(rows, .5)["dataset_id"], "a")
        with self.assertRaisesRegex(ValueError, "no eligible"):
            scan.select_row(rows, .7)
        rows[0]["changedParameters"] = 1
        self.assertEqual(scan.select_row(rows, .5)["dataset_id"], "z")

    def test_hit_whitelist_tightening_and_event_policy(self):
        reference = dict(cuts={"minNHitsFit": 15, "minNHitsRatio": .52, "nSigmaProton": 3},
                         event_meta={"mode": "fxtmult", "KFEventSelection.mainconf": "one"})
        candidate = copy.deepcopy(reference)
        candidate["cuts"].update(minNHitsFit=20, minNHitsRatio=.6)
        candidate["event_meta"]["KFEventSelection.mainconf"] = "two"
        changes = scan.audit_source_differences(candidate, reference)
        self.assertEqual(set(changes), {"minNHitsFit", "minNHitsRatio"})
        candidate["cuts"]["minNHitsFit"] = 10
        with self.assertRaisesRegex(ValueError, "must tighten"):
            scan.audit_source_differences(candidate, reference)
        candidate["cuts"]["minNHitsFit"] = 20
        candidate["cuts"]["nSigmaProton"] = 2.5
        with self.assertRaisesRegex(ValueError, "unapproved effective difference"):
            scan.audit_source_differences(candidate, reference)
        candidate["cuts"]["nSigmaProton"] = 3
        candidate["event_meta"]["mode"] = "refmult"
        with self.assertRaisesRegex(ValueError, "event policies"):
            scan.audit_source_differences(candidate, reference)

    def test_joint_pipeline_common_denominator_and_provenance(self):
        old_files = [REPO / "common/macro/studyLambdaKfPurity.py", REPO / "common/macro/refineLambdaKfPurity.py"]
        old_hashes = [scan.base.sha256(path) for path in old_files]
        with tempfile.TemporaryDirectory(prefix="lambda-upstream-test-") as tmp:
            folder = Path(tmp)
            config, datasets, export, upstream = fixture(folder)
            out = folder / "result"
            result = scan.run(config, datasets, out)
            self.assertEqual(result["datasets"], 2)
            self.assertEqual(result["totalPoints"], 8)
            self.assertEqual(result["selectionSplit"], "full")
            self.assertAlmostEqual(result["baseline"]["default"]["S"], 585.)
            choices = json.loads((out / "chosen_points.json").read_text())["choices"]
            best = next(x for x in choices if x["choice_id"] == "best50")
            self.assertEqual(best["dataset_id"], "hits")
            self.assertEqual(best["combinedCutOverrides"], {"minNHitsFit": 20, "minNHitsRatio": .6})
            with (out / "chosen_metrics.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            row = next(r for r in rows if r["choice_id"] == "best50" and r["sideband"] == "default")
            self.assertAlmostEqual(float(row["signalRetention"]), 360. / 585.)
            self.assertLess(float(row["conservativePurity"]), 1.)
            self.assertEqual({r["split"] for r in rows}, {"full"})
            self.assertTrue(json.loads((out / "status.json").read_text())["complete"])
            event_counts = np.load(out / "chosen_event_counts.npz")
            self.assertEqual(event_counts["event_ids"].shape, (100, 2))
            self.assertTrue((event_counts["counts"][:, 90:, :] == 0).all())
            self.assertTrue((out / "verification/hits_cuts.tsv").exists())
            self.assertNotIn("baseline\t", (out / "verification/hits_cuts.tsv").read_text())
            with (out / "paired_event_bootstrap.csv").open() as stream:
                boot = list(csv.DictReader(stream))
            self.assertEqual({r["split"] for r in boot}, {"full"})
            self.assertTrue(all(float(r["bootstrapError"]) == 0 for r in boot
                                if r["choice_id"] == "baseline" and r["quantity"] == "signalRetention"))
            for line in (out / "SHA256SUMS").read_text().splitlines():
                expected, name = line.split("  ", 1)
                self.assertEqual(expected, scan.base.sha256(out / name))
            with self.assertRaisesRegex(ValueError, "refusing overwrite"):
                scan.run(config, datasets, out)
        self.assertEqual(old_hashes, [scan.base.sha256(path) for path in old_files])

    def test_exact_histogram_replay_required_for_every_upstream_source(self):
        with tempfile.TemporaryDirectory(prefix="lambda-upstream-test-") as tmp:
            folder = Path(tmp)
            config, datasets, export, upstream = fixture(folder)
            hist = scan.base.read_numeric_csv(upstream / "selected_histogram.csv")
            hist["count"][65] += 1
            rows = [dict(zip(hist, values)) for values in zip(*hist.values())]
            scan.base.write_csv(upstream / "selected_histogram.csv", list(hist), rows)
            with self.assertRaisesRegex(ValueError, "histogram does not exactly replay"):
                scan.run(config, datasets, folder / "result")
            self.assertFalse((folder / "result").exists())

    def test_empty_candidate_csv_fails_clearly(self):
        with tempfile.TemporaryDirectory(prefix="lambda-upstream-test-") as tmp:
            folder = Path(tmp)
            config, datasets, export, upstream = fixture(folder)
            path = export / "candidates.csv"
            path.write_text(path.read_text().splitlines()[0] + "\n")
            with self.assertRaisesRegex(ValueError, "empty CSV"):
                scan.run(config, datasets, folder / "result")

    def test_nonpositive_baseline_signal_is_not_a_retention_denominator(self):
        with tempfile.TemporaryDirectory(prefix="lambda-upstream-test-") as tmp:
            folder = Path(tmp)
            config, datasets, export, upstream = fixture(folder)
            policy, c, base_path = scan.load_config(config)
            table = scan.base.read_numeric_csv(export / "candidates.csv")
            table["mass"][:] = 1.1025
            rows = [dict(zip(table, values)) for values in zip(*table.values())]
            scan.base.write_csv(export / "candidates.csv", list(table), rows)
            counts = np.bincount(scan.base.mass_bins(table["mass"][table["selected"] == 1], c), minlength=202)
            hist = [dict(bin=i, low=1.05 + (i - 1) * .001, high=1.05 + i * .001,
                         count=int(counts[i]), error=math.sqrt(counts[i])) for i in range(202)]
            scan.base.write_csv(export / "selected_histogram.csv", list(hist[0]), hist)
            with self.assertRaisesRegex(ValueError, "empty/nonpositive sideband-subtracted signal"):
                scan.run(config, datasets, folder / "result")
            self.assertFalse((folder / "result").exists())

    def test_headered_and_legacy_event_ledger_without_lost_first_event(self):
        with tempfile.TemporaryDirectory(prefix="lambda-upstream-test-") as tmp:
            path = Path(tmp) / "events.tsv"
            c = {"eventIdsFile": str(path), "inputEvents": 2}
            for header in ("", "index\trunId\teventId\n"):
                path.write_text(header + "0\t100\t321\n1\t100\t654\n")
                events, resolved = scan.read_events(c)
                self.assertEqual(events, [(100, 321), (100, 654)])
            path.write_text("index\trunId\teventId\n1\t100\t321\n2\t100\t654\n")
            with self.assertRaisesRegex(ValueError, "unordered"):
                scan.read_events(c)

    def test_no_hidden_global_whitelist_mutation(self):
        self.assertNotIn("minNHitsFit", scan.base.UPSTREAM_DIRECTIONS)
        self.assertNotIn("minNHitsRatio", scan.base.UPSTREAM_DIRECTIONS)


if __name__ == "__main__":
    unittest.main(verbosity=2)
