#!/usr/bin/env python3
"""Synthetic file-only checks of frozen job generation; no ROOT or events."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import yaml

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "common/macro"))
import freezeLambdaKfUpstreamStudy as freeze


def fixture(root):
    master = dict(energies={"auau13p5": {}, "auau3p85": {}}, studyTag="synthetic",
                  scanEvents=10000, validationEvents=100000)
    for energy in master["energies"]:
        out = root / "rootfile" / (energy + "_anaLambda_KFParticle") / "synthetic"
        (out / "configs").mkdir(parents=True)
        (out / "scan").mkdir()
        original = dict(schemaVersion=2, minNHitsFit=15)
        (out / "configs/default_kf.yaml").write_text(yaml.safe_dump(original))
        (out / "configs/default_analysis.yaml").write_text(yaml.safe_dump(dict(dataset={}, analysis={})))
        (out / "configs/default_main_resolved.yaml").write_text(yaml.safe_dump(dict(kf="old", analysis="old")))
        (out / "scan/status.json").write_text(json.dumps(dict(complete=True)))
        best = dict(choice_id="best50", dataset_id="hits20", combinedCutOverrides={"minNHitsFit": 20, "minDecayLengthSignificance": 10.})
        (out / "scan/chosen_points.json").write_text(json.dumps(dict(choices=[best])))
        effective = dict(original, minDecayLengthSignificance=-1.)
        (out / "scan/sources.json").write_text(json.dumps(dict(datasets=[dict(effectiveCuts=effective)])))
    return master


class FreezeJobsTest(unittest.TestCase):
    def test_optional_effective_key_and_job_sample_separation(self):
        with tempfile.TemporaryDirectory(prefix="lambda-freeze-jobs-") as tmp:
            root = Path(tmp)
            master = fixture(root)
            with patch.object(freeze, "REPO", root): freeze.freeze(master)
            for energy in master["energies"]:
                out = root / "rootfile" / (energy + "_anaLambda_KFParticle") / "synthetic"
                folder = out / "configs/frozen"
                best = yaml.safe_load((folder / "best50_kf.yaml").read_text())
                self.assertEqual(best["minDecayLengthSignificance"], 10.)
                self.assertNotIn("minDecayLengthSignificance", yaml.safe_load((out / "configs/default_kf.yaml").read_text()))
                jobs = json.loads((folder / "jobs.json").read_text())
                self.assertEqual([j["events"] for j in jobs], [10000, 100000, 100000])
                self.assertTrue(jobs[0]["input"].endswith("/scan.list"))
                self.assertTrue(all(j["input"].endswith("/validation.list") for j in jobs[1:]))
                self.assertEqual(jobs[1]["overrides"], {})
                before = (folder / "best50_kf.yaml").read_bytes()
                with patch.object(freeze, "REPO", root), self.assertRaisesRegex(RuntimeError, "Refusing existing"):
                    freeze.freeze(master)
                self.assertEqual((folder / "best50_kf.yaml").read_bytes(), before)

    def test_both_discovery_sets_must_be_complete_before_writes(self):
        with tempfile.TemporaryDirectory(prefix="lambda-freeze-jobs-") as tmp:
            root = Path(tmp)
            master = fixture(root)
            out = root / "rootfile/auau3p85_anaLambda_KFParticle/synthetic"
            (out / "scan/status.json").write_text(json.dumps(dict(complete=False)))
            with patch.object(freeze, "REPO", root), self.assertRaisesRegex(RuntimeError, "Incomplete discovery"):
                freeze.freeze(master)
            self.assertFalse(any(root.rglob("frozen")))

    def test_unknown_override_rejected_without_output(self):
        with tempfile.TemporaryDirectory(prefix="lambda-freeze-jobs-") as tmp:
            root = Path(tmp)
            master = fixture(root)
            path = root / "rootfile/auau13p5_anaLambda_KFParticle/synthetic/scan/chosen_points.json"
            chosen = json.loads(path.read_text())
            chosen["choices"][0]["combinedCutOverrides"]["typoCut"] = 10
            path.write_text(json.dumps(chosen))
            with patch.object(freeze, "REPO", root), self.assertRaisesRegex(RuntimeError, "Unknown cut"):
                freeze.freeze(master)
            self.assertFalse(any(root.rglob("frozen")))


if __name__ == "__main__":
    unittest.main(verbosity=2)
