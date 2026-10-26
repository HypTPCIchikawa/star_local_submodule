#!/usr/bin/env python3
"""Synthetic tests for fixed-cut independent-statistics comparison."""
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
import compareLambdaKfFrozenPurity as frozen

SPEC = importlib.util.spec_from_file_location("upstream_test_fixture", REPO / "tests/test_lambda_kf_upstream_purity.py")
FIXTURE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(FIXTURE)


def fixture(folder):
    config, datasets, baseline, best = FIXTURE.fixture(folder)
    validation = folder / "validation_events.tsv"
    validation.write_text("index\trunId\teventId\n" + (folder / "events.tsv").read_text())
    # The discovery ledger is disjoint; the actual exported candidate events
    # are the independent validation input for this no-ROOT fixture.
    (folder / "events.tsv").write_text("".join("{}\t1\t{}\n".format(i, i + 1001) for i in range(100)))
    head, cuts, event = frozen.base.read_metadata(baseline / "metadata.txt")
    saved = folder / "discovery"; saved.mkdir()
    choice = dict(choice_id="best50", combinedCutOverrides={"minNHitsFit": 20, "minNHitsRatio": .6})
    choice_path = saved / "chosen_points.json"
    frozen.base.write_json(choice_path, dict(selectionSplit="full", baselineDataset="baseline", choices=[choice],
        decisions=[dict(choice_id="best50", criterion="robust", fullMetrics=dict(S=360., purity=1., signalRetention=360. / 585.))]))
    frozen.base.write_json(saved / "sources.json", dict(datasets=[dict(id="baseline", effectiveCuts=cuts, eventPolicy=event, meta=head)]))
    return config, baseline, best, validation, choice_path


class FrozenPurityTest(unittest.TestCase):
    def test_block_bootstrap_matches_dense_shared_weights(self):
        counts = np.arange(2 * 20 * 7, dtype=np.int64).reshape(2, 20, 7) % 4
        replicates, seed = 23, 1234
        weights = np.random.RandomState(seed).poisson(1, size=(replicates, 20))
        expected = weights.dot(counts.transpose(1, 0, 2).reshape(20, -1)).reshape(replicates, 2, 7)
        for block in (1, 7, 10, 23):
            self.assertTrue(np.array_equal(frozen.bootstrap_sums(counts, replicates, seed, block), expected))

    def test_complete_actual_selection_comparison_and_no_overwrite(self):
        with tempfile.TemporaryDirectory(prefix="lambda-frozen-test-") as tmp:
            folder = Path(tmp)
            config, baseline, best, events, choice = fixture(folder)
            out = folder / "comparison"
            result = frozen.run(config, baseline, best, events, 100, choice, out)
            self.assertEqual(result["phase"], "independent_frozen_validation")
            self.assertEqual(result["sharedDiscoveryEvents"], 0)
            self.assertFalse(result["cutReselection"])
            self.assertAlmostEqual(result["metrics"]["best50"]["default"]["signalRetention"], 360. / 585.)
            self.assertEqual(result["metrics"]["best50"]["default"]["Nwin"], 360.)
            self.assertTrue(json.loads((out / "status.json").read_text())["complete"])
            counts = np.load(out / "chosen_event_counts.npz")
            self.assertEqual(counts["counts"].shape, (2, 100, 7))
            self.assertTrue((counts["counts"][:, 90:, :] == 0).all())
            compare = json.loads((out / "discovery_validation_comparison.json").read_text())
            self.assertEqual(compare["signalPerInputEventRatio_validation_over_discovery"], 1.)
            self.assertIn("independent_frozen_validation", (out / "scan_config.tsv").read_text())
            self.assertFalse((out / "ranking.csv").exists())
            for line in (out / "SHA256SUMS").read_text().splitlines():
                expected, name = line.split("  ", 1)
                self.assertEqual(expected, frozen.base.sha256(out / name))
            with self.assertRaisesRegex(ValueError, "refusing overwrite"):
                frozen.run(config, baseline, best, events, 100, choice, out)

    def test_concise_beam_energy_titles_without_exploration_wording(self):
        for token, expected in (("auau13p5", "Au+Au FXT; 13.5 GeV beam; 100 input events; frozen cuts"),
                                ("auau3p85", "Au+Au FXT; 3.85 GeV beam; 100 input events; frozen cuts"),
                                ("synthetic", "KF Lambda; 100 independent input events; frozen cuts")):
            with tempfile.TemporaryDirectory(prefix="lambda-frozen-title-test-") as tmp:
                folder = Path(tmp)
                config, baseline, best, events, choice = fixture(folder)
                policy = yaml.safe_load(config.read_text())
                policy["title"] = token + " FXT: upstream/downstream KF Lambda purity exploration"
                config.write_text(yaml.safe_dump(policy, sort_keys=False))
                out = folder / "comparison"
                frozen.run(config, baseline, best, events, 100, choice, out)
                title = next(line.split("\t", 1)[1] for line in (out / "plot_manifest.tsv").read_text().splitlines()
                             if line.startswith("title\t"))
                self.assertEqual(title, expected)
                self.assertNotIn("exploration", title)

    def test_discovery_overlap_and_wrong_event_counts_rejected(self):
        with tempfile.TemporaryDirectory(prefix="lambda-frozen-test-") as tmp:
            folder = Path(tmp)
            config, baseline, best, events, choice = fixture(folder)
            with self.assertRaisesRegex(ValueError, "wrong actual validation input-event count"):
                frozen.run(config, baseline, best, events, 101, choice, folder / "wrong_count")
            (folder / "events.tsv").write_text(events.read_text())
            with self.assertRaisesRegex(ValueError, "overlap"):
                frozen.run(config, baseline, best, events, 100, choice, folder / "overlap")

    def test_extra_actual_cut_or_modified_frozen_cut_rejected(self):
        with tempfile.TemporaryDirectory(prefix="lambda-frozen-test-") as tmp:
            folder = Path(tmp)
            config, baseline, best, events, choice = fixture(folder)
            selected = json.loads(choice.read_text())
            selected["choices"][0]["combinedCutOverrides"]["minNHitsFit"] = 25
            frozen.base.write_json(choice, selected)
            with self.assertRaisesRegex(ValueError, "exactly match frozen overrides"):
                frozen.run(config, baseline, best, events, 100, choice, folder / "mismatch")

    def test_discovery_default_configuration_is_anchored(self):
        with tempfile.TemporaryDirectory(prefix="lambda-frozen-test-") as tmp:
            folder = Path(tmp)
            config, baseline, best, events, choice = fixture(folder)
            source = choice.parent / "sources.json"
            saved = json.loads(source.read_text())
            saved["datasets"][0]["effectiveCuts"]["minNHitsFit"] = 12
            frozen.base.write_json(source, saved)
            with self.assertRaisesRegex(ValueError, "default cuts differ"):
                frozen.run(config, baseline, best, events, 100, choice, folder / "wrong_default")

    def test_reconstructed_count_mismatch_fails(self):
        with tempfile.TemporaryDirectory(prefix="lambda-frozen-test-") as tmp:
            folder = Path(tmp)
            config, baseline, best, events, choice = fixture(folder)
            header, cuts, policy = frozen.base.read_metadata(best / "metadata.txt")
            header["reconstructedEvents"] = 99
            FIXTURE.write_metadata(best / "metadata.txt", header, cuts, policy)
            with self.assertRaisesRegex(ValueError, "reconstructed counts differ"):
                frozen.run(config, baseline, best, events, 100, choice, folder / "wrong_reco")


if __name__ == "__main__":
    unittest.main(verbosity=2)
