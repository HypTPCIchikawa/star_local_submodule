"""Host-only metadata tests. Fixtures are synthetic, not physical PID values."""

from __future__ import annotations

import copy
import importlib.util
from pathlib import Path
import tempfile
import unittest


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("restart_render_unit_subject", HERE / "render_tag_study.py")
subject = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(subject)


class RendererMetadataTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="restart-render-test-")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.nodes = self.directory / "nodes.tsv"
        self.identity = "synthetic_he4_tagd"
        self.write_nodes()
        self.plan = {
            "schemaVersion": 1,
            "productionEnabled": False,
            "coefficientsUpdated": False,
            "purityAvailable": False,
            "totalFits": 1,
            "inputs": {str(self.nodes): subject.digest(self.nodes)},
            "slices": [{
                "id": self.identity, "target": "he4", "mode": "TAG_D",
                "rBin": 88, "rLo": 0.87, "rHi": 0.88, "r": 0.875,
                "fitLo": 35.0, "fitHi": 45.0,
            }],
        }
        components = []
        for species, mean, sigma in (("d", 10.0, 1.0), ("t", 20.0, 1.5),
                                     ("he3", 30.0, 2.0), ("he4", 40.0, 2.0)):
            target = species == "he4"
            components.append({
                "species": species, "active": target, "freeMean": target, "freeSigma": target,
                "initialMean": mean, "initialSigma": sigma,
                "mean": mean + (0.3 if target else 0.0),
                "sigma": sigma + (0.1 if target else 0.0),
            })
        self.fits = {
            "productionEnabled": False, "purityAvailable": False,
            "rawReadback": True, "predictionReadback": True, "totalSlices": 1,
            "results": [{
                "id": self.identity, "target": "he4", "mode": "TAG_D", "rBin": 88,
                "rLo": 0.87, "rHi": 0.88, "fitLo": 35.0, "fitHi": 45.0,
                "fitCompleted": True, "curveAdoptionAllowed": False,
                "purityAvailable": False, "purity": None,
                "components": components, "pearsonPerNdf": 1.2,
                "status": 0, "finite": True, "pointAccepted": True,
                "amplitudeConverged": True, "atBoundary": False,
                "shapeCovarianceValid": True, "yieldCovarianceValid": True,
                "foreignPeakInsideFit": False,
            }],
        }
        self.completed = {
            "status": "completed", "inputsUnchanged": True,
            "productionEnabled": False, "coefficientsUpdated": False, "purityAvailable": False,
        }

    def write_nodes(self, omitted=None, weight=1.0):
        rows = ["id\tspecies\tweight\tmean\tsigma"]
        for species, mean, sigma in (("d", 10.0, 1.0), ("t", 20.0, 1.5), ("he3", 30.0, 2.0)):
            if species != omitted:
                rows.append("{}\t{}\t{}\t{}\t{}".format(self.identity, species, weight, mean, sigma))
        if omitted != "he4":
            rows.extend(["{}\the4\t{}\t39\t2".format(self.identity, weight),
                         "{}\the4\t{}\t41\t2".format(self.identity, weight)])
        self.nodes.write_text("\n".join(rows) + "\n", encoding="utf-8")

    def verify(self, plan=None, fits=None, completed=None):
        return subject.verify_diagnostic(plan or self.plan, fits or self.fits, completed or self.completed)

    def metric(self, fits=None):
        rows, summary = subject.metric_rows(self.plan, fits or self.fits)
        return dict(zip(subject.METRICS_HEADER, rows[0])), summary

    def target(self, fits=None):
        return next(row for row in (fits or self.fits)["results"][0]["components"] if row["species"] == "he4")

    def test_matching_diagnostic_metadata_is_accepted(self):
        self.assertIsNone(self.verify())

    def test_frozen_pre_recenter_backup_is_not_an_executed_node_table(self):
        backup = self.directory / 'pre_recenter_nodes.tsv'
        backup.write_bytes(self.nodes.read_bytes())
        self.plan['inputs'][str(backup)] = subject.digest(backup)
        row, _ = self.metric()
        self.assertEqual(row['mean'], 40.3)

    def test_row_identity_and_order_must_match(self):
        for key, value in (("id", "different"), ("target", "d"), ("mode", "TAG_T"), ("rBin", 89)):
            with self.subTest(key=key):
                fits = copy.deepcopy(self.fits)
                fits["results"][0][key] = value
                with self.assertRaises(ValueError):
                    self.verify(fits=fits)

    def test_reversed_rows_are_rejected(self):
        plan, fits = copy.deepcopy(self.plan), copy.deepcopy(self.fits)
        planned = copy.deepcopy(plan["slices"][0])
        planned["id"] = "second_synthetic_slice"
        result = copy.deepcopy(fits["results"][0])
        result["id"] = planned["id"]
        plan["slices"].append(planned); plan["totalFits"] = 2
        fits["results"].append(result); fits["totalSlices"] = 2
        self.assertIsNone(self.verify(plan=plan, fits=fits))
        fits["results"].reverse()
        with self.assertRaises(ValueError):
            self.verify(plan=plan, fits=fits)

    def test_he4_must_share_deuteron_tof_tag(self):
        for wrong_mode in ("TAG_HE3", "TAG_T", "ALL"):
            with self.subTest(mode=wrong_mode):
                plan, fits = copy.deepcopy(self.plan), copy.deepcopy(self.fits)
                plan["slices"][0]["mode"] = wrong_mode
                fits["results"][0]["mode"] = wrong_mode
                with self.assertRaises(ValueError):
                    self.verify(plan=plan, fits=fits)

    def test_only_target_shapes_may_be_active_or_free(self):
        for species, key, value in (("d", "active", True), ("t", "freeMean", True),
                                    ("he3", "freeSigma", True), ("he4", "active", False),
                                    ("he4", "freeMean", False), ("he4", "freeSigma", False)):
            with self.subTest(species=species, key=key):
                fits = copy.deepcopy(self.fits)
                next(row for row in fits["results"][0]["components"] if row["species"] == species)[key] = value
                with self.assertRaises(ValueError):
                    self.verify(fits=fits)

    def test_missing_or_duplicate_species_are_rejected(self):
        fits = copy.deepcopy(self.fits)
        fits["results"][0]["components"].pop()
        with self.assertRaises(ValueError):
            self.verify(fits=fits)
        fits = copy.deepcopy(self.fits)
        fits["results"][0]["components"][2]["species"] = "d"
        with self.assertRaises(ValueError):
            self.verify(fits=fits)

    def test_completed_fit_geometry_must_match_plan(self):
        for key in ("rLo", "rHi", "fitLo", "fitHi"):
            with self.subTest(key=key):
                fits = copy.deepcopy(self.fits)
                fits["results"][0][key] += 0.1
                with self.assertRaises(ValueError):
                    self.verify(fits=fits)

    def test_fit_failure_accepts_plan_geometry_without_fitted_components(self):
        fits = copy.deepcopy(self.fits)
        fits["results"][0] = {
            "id": self.identity, "target": "he4", "mode": "TAG_D", "rBin": 88,
            "fitCompleted": False, "curveAdoptionAllowed": False,
            "purityAvailable": False, "purity": None, "error": "synthetic nonconvergence",
        }
        self.assertIsNone(self.verify(fits=fits))
        row, summary = self.metric(fits=fits)
        self.assertEqual(row["mean"], 40.0)
        self.assertEqual(row["sigma"], 2.0)
        self.assertEqual(row["fitLo"], 35.0)
        self.assertEqual(row["fitHi"], 45.0)
        self.assertEqual(row["numerical"], 0)
        self.assertEqual(row["chi2Ndf"], "null")
        self.assertIn("original Git PID", row["reason"])
        self.assertEqual(summary["fitExecutionFailures"], 1)

    def test_fit_completion_flag_must_be_explicit_boolean(self):
        for value in (None, 0, 1, "false"):
            with self.subTest(value=value):
                fits = copy.deepcopy(self.fits)
                fits["results"][0]["fitCompleted"] = value
                with self.assertRaises(ValueError):
                    self.verify(fits=fits)

    def test_purity_or_production_adoption_is_rejected(self):
        cases = [("plan", "productionEnabled", True), ("plan", "coefficientsUpdated", True),
                 ("fits", "purityAvailable", True), ("completed", "productionEnabled", True),
                 ("row", "purityAvailable", True), ("row", "curveAdoptionAllowed", True),
                 ("row", "purity", 0.9)]
        for scope, key, value in cases:
            with self.subTest(scope=scope, key=key):
                plan, fits, completed = copy.deepcopy(self.plan), copy.deepcopy(self.fits), copy.deepcopy(self.completed)
                mapping = {"plan": plan, "fits": fits, "completed": completed, "row": fits["results"][0]}[scope]
                mapping[key] = value
                with self.assertRaises(ValueError):
                    self.verify(plan=plan, fits=fits, completed=completed)

    def test_readback_and_immutable_completion_are_required(self):
        for key in ("rawReadback", "predictionReadback"):
            with self.subTest(key=key):
                fits = copy.deepcopy(self.fits); fits[key] = False
                with self.assertRaises(ValueError):
                    self.verify(fits=fits)
        for key, value in (("status", "failed"), ("inputsUnchanged", False)):
            with self.subTest(key=key):
                completed = copy.deepcopy(self.completed); completed[key] = value
                with self.assertRaises(ValueError):
                    self.verify(completed=completed)

    def test_counts_and_duplicate_ids_must_match(self):
        fits = copy.deepcopy(self.fits); fits["totalSlices"] = 2
        with self.assertRaises(ValueError):
            self.verify(fits=fits)
        plan, fits = copy.deepcopy(self.plan), copy.deepcopy(self.fits)
        plan["slices"].append(copy.deepcopy(plan["slices"][0])); plan["totalFits"] = 2
        fits["results"].append(copy.deepcopy(fits["results"][0])); fits["totalSlices"] = 2
        with self.assertRaises(ValueError):
            self.verify(plan=plan, fits=fits)

    def test_valid_metric_uses_fitted_shape_but_does_not_adopt_it(self):
        row, summary = self.metric()
        self.assertEqual(row["mean"], 40.3)
        self.assertEqual(row["sigma"], 2.1)
        self.assertEqual(row["numerical"], 1)
        self.assertEqual(summary["numericallySuccessfulRows"], 1)
        for flag in ("scientificFitsAccepted", "purityAvailable", "productionEnabled", "coefficientsUpdated"):
            self.assertIs(summary[flag], False)

    def test_small_chi2_is_not_scientific_fit_acceptance(self):
        self.fits["results"][0]["pearsonPerNdf"] = 0.001
        row, summary = self.metric()
        self.assertEqual(row["chi2Ndf"], 0.001)
        self.assertIs(summary["scientificFitsAccepted"], False)
        self.assertIs(summary["purityAvailable"], False)
        self.assertIn("Scientific fit acceptance", row["reason"])

    def test_large_chi2_is_preserved_without_inventing_a_quality_cut(self):
        self.fits["results"][0]["pearsonPerNdf"] = 10000.0
        row, summary = self.metric()
        self.assertEqual(row["chi2Ndf"], 10000.0)
        self.assertEqual(row["numerical"], 1)
        self.assertIs(summary["scientificFitsAccepted"], False)

    def test_numerical_failure_flags_are_reported(self):
        changes = [(key, False) for key in ("finite", "pointAccepted", "amplitudeConverged",
                                           "shapeCovarianceValid", "yieldCovarianceValid")]
        changes.extend([("status", 4), ("atBoundary", True)])
        for key, value in changes:
            with self.subTest(key=key):
                fits = copy.deepcopy(self.fits); fits["results"][0][key] = value
                row, summary = self.metric(fits=fits)
                self.assertEqual(row["numerical"], 0)
                self.assertIn(key, row["reason"])
                self.assertEqual(summary["numericallyUnsuccessfulRows"], 1)

    def test_foreign_peak_warning_is_not_hidden(self):
        self.fits["results"][0]["foreignPeakInsideFit"] = True
        row, summary = self.metric()
        self.assertIn("Another original-PID peak", row["reason"])
        self.assertEqual(summary["foreignPeakInsideFitRows"], 1)
        self.assertIs(summary["scientificFitsAccepted"], False)

    def test_nonfinite_shape_uses_original_git_padding_and_is_not_numerical_success(self):
        for mean, sigma in ((None, None), (float("nan"), 2.0), (40.0, -1.0)):
            with self.subTest(mean=mean, sigma=sigma):
                fits = copy.deepcopy(self.fits)
                target = self.target(fits=fits); target["mean"] = mean; target["sigma"] = sigma
                row, summary = self.metric(fits=fits)
                self.assertEqual((row["mean"], row["sigma"]), (40.0, 2.0))
                self.assertEqual(row["numerical"], 0)
                self.assertIn("Invalid fitted shape", row["reason"])
                self.assertIs(summary["coefficientsUpdated"], False)

    def test_initial_shape_must_equal_immutable_git_nodes(self):
        self.target()["initialMean"] += 0.5
        with self.assertRaises(ValueError):
            self.metric()

    def test_missing_and_ambiguous_node_tables_are_rejected(self):
        plan = copy.deepcopy(self.plan); plan["inputs"] = {}
        with self.assertRaises(ValueError):
            subject.seed_shapes(plan)
        extra = self.directory / "another_nodes.tsv"
        extra.write_text(self.nodes.read_text(encoding="utf-8"), encoding="utf-8")
        plan = copy.deepcopy(self.plan); plan["inputs"][str(extra)] = subject.digest(extra)
        with self.assertRaises(ValueError):
            subject.seed_shapes(plan)

    def test_node_weights_and_species_are_validated(self):
        self.write_nodes(omitted="he3")
        with self.assertRaises(ValueError):
            subject.seed_shapes(self.plan)
        self.write_nodes(weight=0.0)
        with self.assertRaises(ValueError):
            subject.seed_shapes(self.plan)

    def test_nonfinite_chi2_is_displayed_as_unavailable(self):
        for value in (None, float("inf"), float("nan")):
            with self.subTest(value=value):
                fits = copy.deepcopy(self.fits); fits["results"][0]["pearsonPerNdf"] = value
                row, summary = self.metric(fits=fits)
                self.assertEqual(row["chi2Ndf"], "null")
                self.assertIs(summary["scientificFitsAccepted"], False)

    def test_reason_cannot_inject_additional_tsv_rows(self):
        fits = copy.deepcopy(self.fits)
        fits["results"][0] = {
            "id": self.identity, "target": "he4", "mode": "TAG_D", "rBin": 88,
            "fitCompleted": False, "curveAdoptionAllowed": False,
            "purityAvailable": False, "purity": None, "error": "line1\tfield\nline2\rline3",
        }
        row, _ = self.metric(fits=fits)
        self.assertFalse(any(character in row["reason"] for character in "\t\n\r"))

    def test_python38_path_containment_and_finite_scalar_helpers(self):
        self.assertTrue(subject.within(self.nodes, self.directory))
        self.assertFalse(subject.within(self.directory.parent, self.directory))
        self.assertTrue(subject.finite(1.5))
        self.assertFalse(subject.finite(True))
        self.assertFalse(subject.finite(float("nan")))


if __name__ == "__main__":
    unittest.main()
