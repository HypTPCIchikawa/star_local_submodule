"""Post-fit identifiability checks use synthetic nodes; no ROOT or fit runs."""
import copy
import csv
from pathlib import Path
import tempfile
import unittest

import compare_response_studies as compare


def make_run(family="gaussian_left_exp", mean=20., sigma=2., alpha=2.,
             lo=15., hi=25., deviance=10., variant="synthetic"):
    identity, target = "t_b0088", "t"
    planned = dict(id=identity, target=target, mode="TAG_T", r=.875, rBin=88,
                   rLo=.87, rHi=.88, fitLo=lo, fitHi=hi)
    shapes = []
    for species in compare.view.SPECIES:
        shapes.append(dict(species=species, initialMean=20., initialSigma=2.,
                           mean=mean if species==target else 20.,
                           sigma=sigma if species==target else 2., meanError=.1, sigmaError=.1))
    row = dict(fitCompleted=True, components=shapes, dataDeviance=deviance,
               pearson=15., pearsonPerNdf=1., engineConditionalNdf=15, nFitBins=20,
               responseAlpha=alpha, responseAlphaError=.2, status=0, atBoundary=False,
               foreignPeakInsideFit=False, curveAdoptionAllowed=False,
               responseOptimizerDiagnosticsAvailable=True, responseMinimizeReturned=True,
               responseMinimizeConverged=True, responseMigradStatus=0, responseMigradEdm=0.,
               responseHesseAttempted=True, responseHesseSucceeded=True, responseHesseStatus=0)
    node = dict(weight=1., mean=20., sigma=2.)
    return dict(directory=Path("/synthetic/" + variant), variant=variant, family=family,
                alphaStart=2., plan={"slices": [planned]}, fits={"results": [row]},
                metrics=[[0, identity, target, "TAG_T", .875, lo, hi, mean, sigma, 1., 1, "compact"]],
                originalShapes={identity: {species: (20., 2.) for species in compare.view.SPECIES}},
                originalNodes={identity: {species: [copy.deepcopy(node)] for species in compare.view.SPECIES}})


class AlphaVisibilityContracts(unittest.TestCase):
    def test_left_join_inside_has_both_regimes(self):
        run = make_run()
        result = compare.alpha_visibility(run, 0)
        self.assertTrue(result["alphaShapeVisibleInFit"])
        self.assertTrue(result["tailVisibleInFit"])
        self.assertTrue(result["coreVisibleInFit"])
        self.assertEqual(result["joinMin"], 16.)
        self.assertTrue(compare.numerical_ok(run, 0))

    def test_left_tail_below_window_is_not_identifiable(self):
        run = make_run(alpha=3.)
        before = copy.deepcopy(run["fits"])
        result = compare.candidate_record(run, 0)
        self.assertFalse(result["alphaShapeVisibleInFit"])
        self.assertFalse(result["numericalOK"])
        self.assertFalse(result["alphaVisibilityDiagnostics"]["tailVisibleInFit"])
        self.assertTrue(result["alphaVisibilityDiagnostics"]["coreVisibleInFit"])
        self.assertEqual(run["fits"], before)

    def test_left_only_exponential_without_core_is_not_identifiable(self):
        result = compare.alpha_visibility(make_run(mean=30., alpha=1.), 0)
        self.assertTrue(result["tailVisibleInFit"])
        self.assertFalse(result["coreVisibleInFit"])
        self.assertFalse(result["alphaShapeVisibleInFit"])

    def test_right_join_inside_has_both_regimes(self):
        result = compare.alpha_visibility(make_run(family="gaussian_right_exp"), 0)
        self.assertTrue(result["alphaShapeVisibleInFit"])
        self.assertEqual(result["joinMin"], 24.)

    def test_right_tail_above_window_is_not_identifiable(self):
        result = compare.alpha_visibility(make_run(family="gaussian_right_exp", alpha=3.), 0)
        self.assertFalse(result["tailVisibleInFit"])
        self.assertTrue(result["coreVisibleInFit"])
        self.assertFalse(result["alphaShapeVisibleInFit"])

    def test_right_only_exponential_without_core_is_not_identifiable(self):
        result = compare.alpha_visibility(make_run(family="gaussian_right_exp", mean=10., alpha=1.), 0)
        self.assertTrue(result["tailVisibleInFit"])
        self.assertFalse(result["coreVisibleInFit"])
        self.assertFalse(result["alphaShapeVisibleInFit"])

    def test_left_join_equal_lower_edge_has_no_tail_interior(self):
        result = compare.alpha_visibility(make_run(alpha=2.5), 0)
        self.assertEqual(result["joinMin"], 15.)
        self.assertFalse(result["tailVisibleInFit"])

    def test_right_join_equal_upper_edge_has_no_tail_interior(self):
        result = compare.alpha_visibility(make_run(family="gaussian_right_exp", alpha=2.5), 0)
        self.assertEqual(result["joinMax"], 25.)
        self.assertFalse(result["tailVisibleInFit"])

    def test_left_join_equal_upper_edge_has_no_core_interior(self):
        result = compare.alpha_visibility(make_run(mean=29.), 0)
        self.assertEqual(result["joinMin"], 25.)
        self.assertFalse(result["coreVisibleInFit"])

    def test_right_join_equal_lower_edge_has_no_core_interior(self):
        result = compare.alpha_visibility(make_run(family="gaussian_right_exp", mean=11.), 0)
        self.assertEqual(result["joinMax"], 15.)
        self.assertFalse(result["coreVisibleInFit"])

    def shifted_gl3(self, family):
        run = make_run(family=family, mean=21., sigma=4., alpha=1.)
        run["originalNodes"]["t_b0088"]["t"] = [
            dict(weight=5/18, mean=18., sigma=1.),
            dict(weight=8/18, mean=20., sigma=2.),
            dict(weight=5/18, mean=22., sigma=3.)]
        return run

    def test_gl3_left_uses_node_shift_and_relative_width_scale(self):
        result = compare.alpha_visibility(self.shifted_gl3("gaussian_left_exp"), 0)
        self.assertEqual([node["mean"] for node in result["fittedNodes"]], [19., 21., 23.])
        self.assertEqual([node["sigma"] for node in result["fittedNodes"]], [2., 4., 6.])
        self.assertEqual([node["join"] for node in result["fittedNodes"]], [17., 17., 17.])
        self.assertTrue(result["alphaShapeVisibleInFit"])

    def test_gl3_right_visibility_uses_extreme_native_joins(self):
        result = compare.alpha_visibility(self.shifted_gl3("gaussian_right_exp"), 0)
        self.assertEqual([node["join"] for node in result["fittedNodes"]], [21., 25., 29.])
        self.assertEqual(result["joinMin"], 21.)
        self.assertEqual(result["joinMax"], 29.)
        self.assertTrue(result["alphaShapeVisibleInFit"])

    def test_gaussian_alpha_visibility_is_not_applicable(self):
        run = make_run(family="gaussian", alpha=None)
        result = compare.alpha_visibility(run, 0)
        self.assertFalse(result["applicable"])
        self.assertIsNone(result["alphaShapeVisibleInFit"])
        self.assertTrue(compare.numerical_ok(run, 0))

    def test_missing_alpha_is_unassessed_not_a_fake_zero(self):
        run = make_run(alpha=None)
        result = compare.candidate_record(run, 0)
        self.assertIsNone(result["alphaShapeVisibleInFit"])
        self.assertIsNone(result["alpha"])
        self.assertFalse(result["numericalOK"])

    def test_reported_problematic_geometry_is_rejected_despite_good_hesse_flags(self):
        run = make_run(mean=22.806, sigma=1.673, alpha=2.5006, lo=19.4, hi=28.)
        run["fits"]["results"][0]["responseAlphaError"] = 458.
        result = compare.candidate_record(run, 0)
        self.assertLess(result["alphaVisibilityDiagnostics"]["joinMax"], 19.4)
        self.assertFalse(result["numericalOK"])
        self.assertTrue(result["optimizerFlags"]["responseHesseSucceeded"])
        self.assertEqual(result["optimizerFlags"]["responseMigradStatus"], 0)

    def test_visibility_is_not_an_implicit_alpha_precision_threshold(self):
        run = make_run(alpha=1.)
        run["fits"]["results"][0]["responseAlphaError"] = 458.
        self.assertTrue(compare.numerical_ok(run, 0))

    def test_display_prefers_identifiable_solution_even_at_higher_deviance(self):
        invisible = make_run(alpha=3., deviance=1., variant="invisible")
        visible = make_run(alpha=1., deviance=20., variant="visible")
        run, selected, metrics = compare.select_display([invisible, visible], 0)
        self.assertIs(run, visible)
        self.assertTrue(selected["numericalOK"])
        self.assertFalse(selected["scientificAdoptionAllowed"])
        self.assertEqual(metrics[10], 1)
        self.assertFalse(selected["candidates"][0]["alphaShapeVisibleInFit"])

    def test_no_identifiable_seed_uses_explicit_failed_display_not_fake_adoption(self):
        first = make_run(alpha=3., deviance=20., variant="first")
        second = make_run(alpha=4., deviance=10., variant="second")
        saved = copy.deepcopy([first["fits"], second["fits"]])
        run, selected, metrics = compare.select_display([first, second], 0)
        self.assertIs(run, second)
        self.assertFalse(selected["numericalOK"])
        self.assertFalse(selected["scientificAdoptionAllowed"])
        self.assertEqual(selected["displaySelection"], "minimum_finite_deviance_failed_solution_for_display")
        self.assertIn("tail/core shape not visible; alpha unidentifiable", metrics[11])
        self.assertEqual(metrics[10], 0)
        self.assertEqual([first["fits"], second["fits"]], saved)

    def test_raw_only_fallback_keeps_missing_values_unassessed(self):
        run = make_run()
        run["fits"]["results"][0] = dict(fitCompleted=False, dataDeviance=None, error="synthetic empty sample")
        run["metrics"][0][10] = 0
        _, selected, metrics = compare.select_display([run], 0)
        self.assertEqual(selected["displaySelection"], "raw_only_fallback")
        self.assertIsNone(selected["alphaShapeVisibleInFit"])
        self.assertIsNone(selected["dataDeviance"])
        self.assertFalse(selected["scientificAdoptionAllowed"])
        self.assertIn("RAW_ONLY_FALLBACK=1", metrics[11])


class ImmutableNodeContracts(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.path = Path(self.temp.name) / "nodes.tsv"
        self.run = make_run()

    def tearDown(self):
        self.temp.cleanup()

    def write(self, sigma=2., identity="t_b0088"):
        with self.path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=compare.view.NODE_HEADER, delimiter="\t")
            writer.writeheader()
            for species in compare.view.SPECIES:
                writer.writerow(dict(id=identity, species=species, weight=1., mean=20., sigma=sigma))

    def test_all_four_known_curves_are_read_without_altering_nodes(self):
        self.write()
        nodes = compare.read_original_nodes(self.path, self.run["plan"], self.run["originalShapes"])
        self.assertEqual(set(nodes["t_b0088"]), compare.view.SPECIES)
        self.assertEqual(nodes["t_b0088"]["t"], [dict(weight=1., mean=20., sigma=2.)])

    def test_nonpositive_original_sigma_is_fatal(self):
        self.write(sigma=0.)
        with self.assertRaisesRegex(ValueError, "Invalid immutable"):
            compare.read_original_nodes(self.path, self.run["plan"], self.run["originalShapes"])

    def test_wrong_case_id_is_fatal(self):
        self.write(identity="wrong")
        with self.assertRaisesRegex(ValueError, "ids differ"):
            compare.read_original_nodes(self.path, self.run["plan"], self.run["originalShapes"])

    def test_seed_weighted_width_mismatch_is_fatal(self):
        self.write(sigma=3.)
        with self.assertRaisesRegex(ValueError, "averages"):
            compare.read_original_nodes(self.path, self.run["plan"], self.run["originalShapes"])


if __name__ == "__main__":
    unittest.main()
