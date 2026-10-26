"""Synthetic saved-fit contracts; the fake ROOT bytes are never ROOT-read."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


PATH = Path(__file__).with_name("summarize_tag_study.py")
SPEC = importlib.util.spec_from_file_location("restart_tag_summary", PATH)
summary = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(summary)


class SummaryContracts(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.fit_dir = self.root / "fits"
        self.fit_dir.mkdir()
        self.cache = self.root / "immutable_cache.root"
        self.cache.write_bytes(b"synthetic immutable cache bytes; not a ROOT file")
        self.project = self.root / "project"
        self.project.mkdir()
        self.source = self.project / "source.C"
        self.source.write_text("// synthetic frozen source\n")
        self.source_snapshot = self.fit_dir / "source/source.C"
        self.source_snapshot.parent.mkdir()
        self.source_snapshot.write_bytes(self.source.read_bytes())
        (self.fit_dir / "native_spectra.root").write_bytes(b"synthetic native spectra; not ROOT-read")
        self.plan = dict(schemaVersion=1, totalFits=4,
                         config={"baselinePlots": {"rLowMax": 2, "rHighMin": 2}},
                         inputs={str(self.cache): summary.digest(self.cache)},
                         sourceSha256={str(self.source): summary.digest(self.source)},
                         cache=str(self.cache), inputFingerprint="synthetic input fingerprint",
                         slices=[])
        self.fits = dict(totalSlices=4, rawReadback=True, predictionReadback=True,
                         productionEnabled=False, purityAvailable=False, results=[])
        for index, (r, chi2) in enumerate(((0.875, 1), (1.775, 9), (2.0, 2), (2.5, 4))):
            planned = dict(id="slice_" + str(index), target="he4", mode="gaussian_linear",
                           rBin=index, rLo=r - 0.02, rHi=r + 0.02, r=r,
                           fitLo=20, fitHi=60, fitWindowAxisClipped=False)
            self.plan["slices"].append(planned)
            components = []
            for species in summary.SPECIES:
                active = species == "he4"
                components.append(dict(species=species, active=active, freeMean=active,
                                       freeSigma=active, initialMean=40, initialSigma=3,
                                       mean=41, sigma=3.1, meanError=0.1, sigmaError=0.2,
                                       histogram="component_" + species))
            row = {key: planned[key] for key in ("id", "mode", "target", "rBin", "rLo", "rHi",
                                                "fitLo", "fitHi")}
            row.update(fitCompleted=True, fullNativeCounts=1000 * (index + 1),
                       fitNativeCounts=500, nFitBins=20, status=0, covarianceStatus=3,
                       finite=True, pointAccepted=True, amplitudeConverged=True, atBoundary=False,
                       yieldCovarianceValid=True, shapeCovarianceValid=True,
                       foreignPeakInsideFit=index == 1, curveAdoptionAllowed=False,
                       purityAvailable=False, purity=None, engineConditionalNdf=10,
                       fittedShapeParameterCount=2, pearson=10 * chi2, pearsonPerNdf=chi2,
                       components=components)
            self.fits["results"].append(row)
        self.completion = dict(status="completed", inputsUnchanged=True, productionEnabled=False,
                               coefficientsUpdated=False, purityAvailable=False)
        self.write()

    def tearDown(self):
        self.temp.cleanup()

    def write(self):
        for name, data in (("plan.json", self.plan), ("fits.json", self.fits),
                           ("fit.log.completed.json", self.completion)):
            (self.fit_dir / name).write_text(json.dumps(data, allow_nan=False))

    def run_summary(self):
        self.write()
        return summary.summarize(self.fit_dir, project=self.project)

    def target(self, row=0):
        return next(component for component in self.fits["results"][row]["components"]
                    if component["species"] == "he4")

    def test_complete_is_not_scientific_acceptance(self):
        result = self.run_summary()
        self.assertTrue(result["technicalExecutionCompleted"])
        self.assertTrue(result["allSlicesFitCompleted"])
        self.assertEqual(result["totals"]["numericalOK"], 4)
        self.assertFalse(result["curveAdoptionAllowed"])
        self.assertFalse(result["scientificAcceptanceSelected"])
        self.assertFalse(result["purityAvailable"])
        self.assertTrue(result["noRootReadbackPerformedBySummary"])
        self.assertIn("nativeRoot", result["hashesAtSummaryCreation"])

    def test_bad_status_remains_unassessed(self):
        self.fits["results"][0]["status"] = 4
        result = self.run_summary()
        self.assertEqual(result["totals"]["numericalOK"], 3)
        self.assertEqual(result["rows"][0]["assessment"], "unassessed")
        self.assertEqual(result["rows"][0]["chi2PerNdf"], 1)

    def test_null_chi2_is_not_zero(self):
        self.fits["results"][0]["pearsonPerNdf"] = None
        result = self.run_summary()
        self.assertIsNone(result["rows"][0]["chi2PerNdf"])
        self.assertFalse(result["rows"][0]["numericalOK"])
        self.assertIn("chi2_unassessed", result["rows"][0]["reasons"])

    def test_negative_chi2_is_unassessed(self):
        self.fits["results"][0]["pearsonPerNdf"] = -1
        result = self.run_summary()
        self.assertIsNone(result["rows"][0]["chi2PerNdf"])
        self.assertFalse(result["rows"][0]["numericalOK"])

    def test_row_order_mismatch_is_fatal(self):
        self.fits["results"][0], self.fits["results"][1] = \
            self.fits["results"][1], self.fits["results"][0]
        with self.assertRaisesRegex(ValueError, "row order/identity"):
            self.run_summary()

    def test_partial_mean_width_never_becomes_complete_curve(self):
        self.target()["sigma"] = None
        result = self.run_summary()
        row = result["rows"][0]
        self.assertEqual(row["fittedMean"], 41)
        self.assertIsNone(row["fittedSigma"])
        self.assertFalse(row["numericalOK"])
        self.assertIn("invalid_or_partial_target_shape", row["reasons"])

    def test_global_and_partition_medians_are_computed_from_rows(self):
        result = self.run_summary()["species"]["he4"]
        key = "chi2PerNdfNumericalOK"
        self.assertEqual(result["allR"][key]["statistics"]["median"], 3)
        self.assertEqual(result["lowR"][key]["statistics"]["median"], 5)
        self.assertEqual(result["highR"][key]["statistics"]["median"], 3)
        self.assertEqual(result["highR"]["total"], 2)  # exact R=2 belongs to TOF/highR.
        self.assertEqual(result["allR"][key]["sensitivity"]["countsAtOrBelow"],
                         {"2": 2, "3": 2, "5": 3, "10": 4})

    def test_at_boundary_is_not_numerical_ok(self):
        self.fits["results"][0]["atBoundary"] = True
        self.assertFalse(self.run_summary()["rows"][0]["numericalOK"])

    def test_changed_protected_input_is_fatal(self):
        self.cache.write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError, "Protected input changed"):
            self.run_summary()

    def test_unverified_prediction_readback_is_fatal(self):
        self.fits["predictionReadback"] = False
        with self.assertRaisesRegex(ValueError, "readback was not verified"):
            self.run_summary()

    def test_curve_adoption_permission_is_fatal(self):
        self.fits["results"][0]["curveAdoptionAllowed"] = True
        with self.assertRaisesRegex(ValueError, "permits curve/purity adoption"):
            self.run_summary()

    def test_nonnull_purity_is_fatal(self):
        self.fits["results"][0]["purity"] = 0.9
        with self.assertRaisesRegex(ValueError, "permits curve/purity adoption"):
            self.run_summary()

    def test_negative_native_counts_are_not_silently_clamped(self):
        self.fits["results"][0]["fullNativeCounts"] = -1
        row = self.run_summary()["rows"][0]
        self.assertIsNone(row["fullNativeCounts"])
        self.assertFalse(row["numericalOK"])

    def test_representative_cases_use_saved_original_seed(self):
        cases = self.run_summary()["representativeCases"]
        self.assertEqual([row["r"] for row in cases], [0.875, 1.775])
        self.assertEqual(cases[0]["originalMean"], 40)
        self.assertEqual(cases[0]["originalSigma"], 3)
        self.assertEqual(cases[0]["fittedMean"], 41)
        self.assertEqual(cases[0]["meanError"], 0.1)

    def test_failed_row_has_null_statistics_not_zero(self):
        planned = self.plan["slices"][0]
        failed = {key: planned[key] for key in ("id", "mode", "target", "rBin")}
        failed.update(fitCompleted=False, purityAvailable=False, purity=None,
                      curveAdoptionAllowed=False, error="synthetic empty slice")
        self.fits["results"][0] = failed
        result = self.run_summary()
        self.assertTrue(result["technicalExecutionCompleted"])
        self.assertFalse(result["allSlicesFitCompleted"])
        self.assertIsNone(result["rows"][0]["chi2PerNdf"])
        self.assertIsNone(result["rows"][0]["fullNativeCounts"])
        self.assertIsNone(result["rows"][0]["fittedMean"])
        self.assertEqual(result["totals"]["fitCompleted"], 3)

    def test_incomplete_execution_is_fatal(self):
        self.completion["status"] = "failed"
        with self.assertRaisesRegex(ValueError, "completion"):
            self.run_summary()

    def test_missing_species_population_is_null_not_zero_median(self):
        result = self.run_summary()["species"]["d"]["allR"]
        self.assertEqual(result["total"], 0)
        self.assertIsNone(result["chi2PerNdfNumericalOK"]["statistics"]["median"])

    def test_existing_summary_output_is_preserved(self):
        output = self.root / "summary.json"
        output.write_text("preserve this existing report")
        with self.assertRaises(FileExistsError):
            summary.save_report(self.run_summary(), output)
        self.assertEqual(output.read_text(), "preserve this existing report")

    def test_changed_live_source_does_not_change_historical_fit(self):
        self.source.write_text("// modified after fit")
        result = self.run_summary()
        verified = result["sourceFilesVerified"][str(self.source)]
        self.assertEqual(verified["authoritative"], "executed_snapshot")
        self.assertFalse(verified["liveSourceCompared"])

    def test_changed_source_snapshot_is_fatal(self):
        self.source_snapshot.write_text("// modified frozen fit source")
        with self.assertRaisesRegex(ValueError, "Protected source snapshot changed"):
            self.run_summary()


if __name__ == "__main__":
    unittest.main()
