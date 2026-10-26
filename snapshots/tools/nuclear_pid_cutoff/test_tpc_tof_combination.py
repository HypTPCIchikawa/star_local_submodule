import math
import unittest

from combine_tpc_tof import check_coverage, checked_rows, nullable_number, quality


def row(lo, hi, **overrides):
    value = dict(species="d", rLo=str(lo), rHi=str(hi), purityNominal="0.99", purityLow="0.98", purityHigh="1", rawSelected="12345", quality="1")
    value.update(overrides)
    return value


class TpcTofCombinationTests(unittest.TestCase):
    def setUp(self):
        self.targets = {"d": {"displayMax": 2.5}}

    def test_frozen_values_are_preserved_not_recomputed(self):
        source = row(.4, .45, purityNominal="0.98341019779536831")
        output = list(checked_rows([source], "TPC", .4, 2, self.targets))[0]
        self.assertEqual(source["purityNominal"], output["purityNominal"])
        self.assertEqual(output["regime"], "TPC")

    def test_unassessed_missing_estimate_stays_missing(self):
        output = list(checked_rows([row(2, 2.05, quality="0", purityNominal="null", purityLow="", purityHigh="nan")], "TOF", 2, 2.5, self.targets))[0]
        self.assertEqual(output["purityNominal"], "null")
        self.assertEqual(output["quality"], "0")

    def test_supported_missing_estimate_is_rejected(self):
        with self.assertRaises(ValueError):
            list(checked_rows([row(2, 2.05, purityLow="null")], "TOF", 2, 2.5, self.targets))

    def test_duplicate_interval_is_rejected(self):
        with self.assertRaises(ValueError):
            list(checked_rows([row(2, 2.05), row(2, 2.05)], "TOF", 2, 2.5, self.targets))

    def test_regime_boundary_cannot_be_crossed(self):
        with self.assertRaises(ValueError):
            list(checked_rows([row(1.99, 2.04)], "TOF", 2, 2.5, self.targets))

    def test_target_upper_domain_is_enforced(self):
        with self.assertRaises(ValueError):
            list(checked_rows([row(2.5, 2.55)], "TOF", 2, 5, self.targets))

    def test_model_envelope_cannot_be_inverted(self):
        with self.assertRaises(ValueError):
            list(checked_rows([row(2, 2.05, purityLow="1", purityHigh="0.9")], "TOF", 2, 2.5, self.targets))

    def test_raw_count_is_required_even_if_fit_failed(self):
        with self.assertRaises(ValueError):
            list(checked_rows([row(2, 2.05, quality="0", rawSelected="nan")], "TOF", 2, 2.5, self.targets))

    def test_coverage_requires_every_interval(self):
        rows = [row(2 + i * .05, 2 + (i + 1) * .05) for i in range(10)]
        check_coverage(rows, "TOF", 2, 2.5, self.targets, .05)
        with self.assertRaises(ValueError):
            check_coverage(rows[:-1], "TOF", 2, 2.5, self.targets, .05)

    def test_quality_has_no_silent_default(self):
        for invalid in ("", "maybe", "nan"):
            with self.assertRaises(ValueError):
                quality(invalid)

    def test_null_is_not_zero_purity(self):
        self.assertTrue(math.isnan(nullable_number("None")))
        self.assertEqual(nullable_number("0"), 0)


if __name__ == "__main__":
    unittest.main()
