"""Synthetic contracts for the scoped restore audit; not physics fit tests."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import yaml


MODULE_PATH = Path(__file__).with_name("audit_restore.py")
SPEC = importlib.util.spec_from_file_location("restart_audit_restore", MODULE_PATH)
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


HEADER = """namespace NuclearIdDeDxVsMom {
double kMean[4][4] = {{10,1,0,0},{20,1,0,0},{30,1,0,0},{40,1,0,0}};
double kSigma[4][4] = {{1,0,0,0},{2,0,0,0},{3,0,0,0},{4,0,0,0}};
double kPRange[4][2] = {{0.2,5},{0.2,5},{0.2,5},{0.2,5}};
} // Synthetic data only.
"""


class RestoreContracts(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.project, self.before = self.root / "project", self.root / "before"
        self.old = dict(m2_selection=False, maxPOverQ=5,
                        speciesMinPOverQ={"d": 0.4},
                        highRigidityTof={"enabled": True, "thresholdPOverQ": 2},
                        dedxCalibration={"enabled": True, "calibrationFile": "calibration/old.yaml"})
        self.new = dict(self.old, dedxCalibration=dict(self.old["dedxCalibration"], enabled=False))
        self.before_main = {"analysis": "analysis/info.yaml", "nuclearid": "cuts/old.yaml"}
        self.after_main = dict(self.before_main, nuclearid="cuts/restored.yaml")
        for tree in (self.project, self.before):
            self.write(tree, "config/mainconf/main.yaml",
                       self.after_main if tree == self.project else self.before_main)
            self.write(tree, "config/analysis/info.yaml", {"starTag": {"year": "year2021"}})
            self.write(tree, "config/cuts/old.yaml", self.old)
        self.write(self.project, "config/cuts/restored.yaml", self.new)
        self.header = self.project / "StMaker/common/NuclearIdDeDxVsMom.h"
        self.header.parent.mkdir(parents=True)
        self.header.write_text(HEADER)
        self.settings = self.root / "settings.yaml"
        self.settings.write_text(yaml.safe_dump(dict(
            schemaVersion=1, gitRevision="explicit-commit",
            legacyHeader="StMaker/common/NuclearIdDeDxVsMom.h",
            mainconfs=[dict(mainconf="config/mainconf/main.yaml",
                           oldNuclearid="cuts/old.yaml", newNuclearid="cuts/restored.yaml")],
            rigiditySamples=[0.1, 0.5, 1, 3], dedxProbeMultipliers=[-2, 0, 2],
            comparisonTolerance={"relative": 1e-12, "absolute": 1e-12})))

    def tearDown(self):
        self.temp.cleanup()

    def write(self, tree, relative, value):
        path = tree / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(yaml.safe_dump(value))

    def run_audit(self):
        with mock.patch.object(audit.subprocess, "check_output",
                               side_effect=["a" * 40 + "\n", HEADER.encode()]):
            return audit.audit(self.project, self.settings, self.before)

    def test_exact_scoped_restore_passes(self):
        result = self.run_audit()
        self.assertEqual(result["status"], "passed")
        self.assertFalse(result["runtimeProviderVerified"])
        self.assertFalse(result["productionFitsExecuted"])
        self.assertTrue(result["mainconfs"][0]["nonPidSelectionUnchanged"])
        json.dumps(result, allow_nan=False)

    def test_active_calibration_is_rejected(self):
        self.write(self.project, "config/cuts/restored.yaml", self.old)
        with self.assertRaisesRegex(ValueError, "true -> false"):
            self.run_audit()

    def test_disabled_calibration_knots_are_rejected(self):
        changed = dict(self.new, dedxCalibration=dict(self.new["dedxCalibration"], knots=[]))
        self.write(self.project, "config/cuts/restored.yaml", changed)
        with self.assertRaisesRegex(ValueError, "inject knots"):
            self.run_audit()

    def test_coefficient_injection_is_rejected(self):
        changed = dict(self.new, pidCoefficients={"he4": [1, 2, 3]})
        self.write(self.project, "config/cuts/restored.yaml", changed)
        with self.assertRaisesRegex(ValueError, "beyond dedxCalibration"):
            self.run_audit()

    def test_cutoff_changes_are_rejected(self):
        self.write(self.project, "config/cuts/restored.yaml", dict(self.new, maxPOverQ=6))
        with self.assertRaisesRegex(ValueError, "beyond dedxCalibration"):
            self.run_audit()

    def test_tof_changes_are_rejected(self):
        changed = dict(self.new, highRigidityTof={"enabled": False, "thresholdPOverQ": 2})
        self.write(self.project, "config/cuts/restored.yaml", changed)
        with self.assertRaisesRegex(ValueError, "beyond dedxCalibration"):
            self.run_audit()

    def test_other_mainconf_change_is_rejected(self):
        self.write(self.project, "config/mainconf/main.yaml",
                   dict(self.after_main, maker="maker/new.yaml"))
        with self.assertRaisesRegex(ValueError, "Mainconf changed beyond"):
            self.run_audit()

    def test_other_dependency_change_is_rejected(self):
        self.write(self.project, "config/analysis/info.yaml", {"changed": True})
        with self.assertRaisesRegex(ValueError, "dependency changed"):
            self.run_audit()

    def test_original_dated_yaml_is_immutable(self):
        self.write(self.project, "config/cuts/old.yaml", self.new)
        with self.assertRaisesRegex(ValueError, "dated original"):
            self.run_audit()

    def test_current_header_must_match_git_bytes(self):
        self.header.write_text(HEADER + "\n// locally changed\n")
        with self.assertRaisesRegex(ValueError, "byte-identical"):
            self.run_audit()

    def test_missing_required_file_is_fatal(self):
        (self.project / "config/analysis/info.yaml").unlink()
        with self.assertRaises(FileNotFoundError):
            self.run_audit()

    def test_missing_settings_have_no_defaults(self):
        settings = yaml.safe_load(self.settings.read_text())
        del settings["rigiditySamples"]
        self.settings.write_text(yaml.safe_dump(settings))
        with self.assertRaisesRegex(ValueError, "missing restore settings"):
            self.run_audit()

    def test_unsupported_domain_is_not_constant_fallback(self):
        rows = self.run_audit()["samples"]
        unsupported = [row for row in rows if row["rigidity"] == 0.1]
        self.assertEqual(len(unsupported), 4)
        self.assertTrue(all(row["mean"] is None and row["sigma"] is None and
                            row["pulls"] is None and not row["supportedByGitDomain"]
                            for row in unsupported))

    def test_positive_mean_and_width_are_explicit(self):
        rows = self.run_audit()["samples"]
        he4 = next(row for row in rows if row["species"] == "he4" and row["rigidity"] == 1)
        self.assertEqual(he4["mean"], 40)
        self.assertEqual(he4["sigma"], 4)
        self.assertEqual(he4["pulls"][0]["dedx"], 32)
        self.assertEqual(he4["pulls"][0]["pull"], -2)

    def test_negative_width_does_not_fall_back(self):
        tables = audit.tables_from_header(HEADER)
        tables["sigma"][3] = [-1, 0, 0, 0]
        with self.assertRaisesRegex(ValueError, "invalid mean or width"):
            audit.sample_curves(tables, [1], [25])

    def test_relative_paths_cannot_escape_tree(self):
        with self.assertRaisesRegex(ValueError, "repository-relative"):
            audit.checked_relative("../outside.yaml")


if __name__ == "__main__":
    unittest.main()
