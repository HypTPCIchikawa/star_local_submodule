#!/usr/bin/env python3
"""Synthetic tests for the separate full-data exploratory refinement.

env -u PYTHONPATH /usr/bin/python3 tests/test_lambda_purity_refinement.py
"""
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
import refineLambdaKfPurity as refine

SPEC = importlib.util.spec_from_file_location("purity_fixture", REPO / "tests/test_lambda_purity_scan.py")
FIXTURE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(FIXTURE)
CONFIG = REPO / "config/qa/lambda_kf_purity70_refinement_20260927.yaml"


class PurityRefinementTest(unittest.TestCase):
    def test_declared_grid_and_previous_reference(self):
        policy, c, base_path = refine.load_config(CONFIG)
        points = refine.base.enumerate_points(c)
        self.assertEqual(len(points), 3756)
        self.assertEqual(len(set(points)), len(points))
        self.assertEqual(points[0], ())
        reference = tuple(sorted((key, float(value)) for key, value in policy["fixedReference"]["cuts"].items()))
        self.assertIn(reference, points)
        self.assertEqual(c["window"], [1.110, 1.122])
        self.assertEqual(policy["minimumSidebandCount"], 0)

    def test_exact_zero_and_one_count_upper_limits(self):
        self.assertAlmostEqual(refine.poisson_upper(0, .9), math.log(10), places=13)
        self.assertAlmostEqual(refine.poisson_upper(1, .9), 3.889720169867429, places=12)
        self.assertGreater(refine.poisson_upper(0, .9), 0)
        self.assertAlmostEqual(refine.poisson_cdf(1, refine.poisson_upper(1, .9)), .1, places=13)

    def test_poisson_discrete_coverage_and_monotonicity(self):
        confidence = .9
        upper = [refine.poisson_upper(k, confidence) for k in range(121)]
        self.assertTrue(all(a < b for a, b in zip(upper, upper[1:])))
        for mean in (.05, .5, 1., 2., 3., 5., 10., 20., 50.):
            rejection = sum(math.exp(-mean + k*math.log(mean)-math.lgamma(k+1))
                            for k, limit in enumerate(upper) if limit < mean)
            self.assertLessEqual(rejection, 1.-confidence+1e-12)
        for k in (0, 1, 10, 50, 120):
            self.assertAlmostEqual(refine.poisson_cdf(k, upper[k]), 1.-confidence, places=12)

    def test_zero_sideband_is_finite_not_background_free_score(self):
        score = 1.-.75*refine.poisson_upper(0, .9)/1000
        self.assertTrue(math.isfinite(score))
        self.assertLess(score, 1.)
        rows = [{"point_id": 1, "selectionEligible": 1, "minSignalRetentionAcrossSidebands": .72,
                 "robustPurityScore": score, "purity": 1., "changedParameters": 1, "S": 1000.}]
        self.assertEqual(refine.select_row(rows, .70)["point_id"], 1)

    def test_worst_sideband_retention_and_distinct_nominal_rule(self):
        common = {"selectionEligible": 1, "changedParameters": 1, "S": 700.}
        rows = [dict(common, point_id=1, minSignalRetentionAcrossSidebands=.71, robustPurityScore=.98, purity=.99),
                dict(common, point_id=2, minSignalRetentionAcrossSidebands=.72, robustPurityScore=.97, purity=.995),
                dict(common, point_id=3, minSignalRetentionAcrossSidebands=.69, robustPurityScore=.999, purity=1.)]
        self.assertEqual(refine.select_row(rows, .70)["point_id"], 1)
        self.assertEqual(refine.select_row(rows, .70, nominal=True)["point_id"], 2)
        with self.assertRaisesRegex(ValueError, "no eligible"):
            refine.select_row(rows, .80)

    def test_old_diagnostic_split_cannot_be_eligible(self):
        policy, c, _ = refine.load_config(CONFIG)
        result = refine.base.measure(100, 0, 0, .75)
        result.update(signalRetention=1., NwinRetention=1., baselineS=100., baselineNwin=100.,
                      poissonUpperMean=math.log(10), backgroundUpper=.75*math.log(10),
                      conservativePurity=1.-.75*math.log(10)/100)
        aggregate = {"robustPurityScore": result["conservativePurity"], "minSignalRetentionAcrossSidebands": 1.}
        d = {"id": "baseline"}
        for split in ("train", "validation"):
            row = refine.make_row(d, 0, (), split, "default", result, aggregate, policy)
            self.assertEqual(row["selectionEligible"], 0)
            self.assertEqual(row["splitRole"], "previously_inspected_stability_diagnostics_only")
        row = refine.make_row(d, 0, (), "full", "default", result, aggregate, policy)
        self.assertEqual(row["selectionEligible"], 1)

    def test_full_synthetic_pipeline_provenance_and_no_overwrite(self):
        untouched = [REPO / "common/macro/studyLambdaKfPurity.py", REPO / "config/qa/lambda_kf_purity_study_20260927.yaml"]
        hashes = [refine.base.sha256(path) for path in untouched]
        with tempfile.TemporaryDirectory(prefix="lambda-purity-refine-tests-") as temp:
            folder = Path(temp)
            c, base_cfg, export, datasets = FIXTURE.fixture(folder)
            policy, original, path = refine.load_config(CONFIG)
            policy["baseConfig"] = str(base_cfg)
            policy["grid"] = {"nSigmaProton": [2.5, 2.0]}
            policy["tripleAnchor"] = "nSigmaProton"
            policy["fixedReference"]["cuts"] = {"nSigmaProton": 2.5}
            policy["minimumWindowCount"] = 10
            config = folder / "refine.yaml"
            config.write_text(yaml.safe_dump(policy, sort_keys=False))
            output = folder / "results"
            summary = refine.run(config, datasets, output)
            self.assertEqual(summary["pointsPerDataset"], 3)
            self.assertEqual(summary["selectionSplit"], "full")
            self.assertIn("NOT a fresh holdout", summary["scope"])
            self.assertTrue(json.loads((output / "status.json").read_text())["complete"])
            sources = json.loads((output / "sources.json").read_text())
            self.assertEqual(sources["policy"]["sha256"], refine.base.sha256(config))
            self.assertEqual(sources["helpers"]["sha256"], hashes[0])
            event_counts = np.load(output / "chosen_event_counts.npz")
            self.assertEqual(event_counts["event_ids"].shape, (100, 2))
            self.assertTrue((event_counts["counts"][:, 90:, :] == 0).all())
            for line in (output / "SHA256SUMS").read_text().splitlines():
                expected, name = line.split("  ", 1)
                self.assertEqual(expected, refine.base.sha256(output/name))
            with self.assertRaisesRegex(ValueError, "refusing overwrite"):
                refine.run(config, datasets, output)
        self.assertEqual(hashes, [refine.base.sha256(path) for path in untouched])


if __name__ == "__main__":
    unittest.main(verbosity=2)

