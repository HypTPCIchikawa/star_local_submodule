#!/usr/bin/env python3
"""Count-only mT quantile selection; no ROOT dependency or CF optimization."""

import copy
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "mt_binning", ROOT / "tools/kfparticle_summary/mt_binning.py")
MT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MT)


def fixture(weights, normalization=None):
    """Saved k* = 50 MeV toy bins, rebinned by two; mT = 50 MeV."""
    nx, ny = 20, len(weights)
    norm = list(weights if normalization is None else normalization)
    data = {"channels": {}}
    for mode in MT.MODES:
        data["channels"][mode] = {}
        for region in MT.REGIONS:
            scale = (1.0 if mode == "SE" else 5.0) * (1.0 if region == "signal" else 0.2)
            counts = [[0.0] * ny for _ in range(nx)]
            for y in range(ny):
                for x in range(6):
                    counts[x][y] = scale * weights[y] / 6.0
                for x in (6, 7):
                    counts[x][y] = scale * norm[y] / 2.0
            error2 = [[value * 1.25 for value in row] for row in counts]
            channel = {"kstarEdges": [i / 20.0 for i in range(nx + 1)],
                       "mtEdges": [i / 20.0 for i in range(ny + 1)],
                       "counts": counts, "error2": error2,
                       "flowCounts": [[0.0] * (ny + 2) for _ in range(nx + 2)],
                       "flowError2": [[0.0] * (ny + 2) for _ in range(nx + 2)]}
            for x in range(nx):
                for y in range(ny):
                    channel["flowCounts"][x + 1][y + 1] = counts[x][y]
                    channel["flowError2"][x + 1][y + 1] = error2[x][y]
            data["channels"][mode][region] = channel
    settings = {"mass_region": "signal", "n_slices": 4,
                "split_kstar": [0.0, 0.3],
                "split_strategy": "nearest_se_signal_quantiles_at_saved_mt_edges",
                "mt_extent": "occupied_regular_bins_all_channels"}
    base = {"rebin_group": 2, "normalization": [0.3, 0.4], "plot_x_range": [0.0, 0.3]}
    return data, settings, base


def set_cell(data, mode, region, x, y, count, error2=None):
    channel = data["channels"][mode][region]
    for regular, flow, value in (("counts", "flowCounts", count),
                                 ("error2", "flowError2", count if error2 is None else error2)):
        channel[regular][x][y] = value
        channel[flow][x + 1][y + 1] = value


def fixed_settings(settings, edges):
    result = dict(settings)
    result.update(split_strategy="fixed_saved_mt_edges", mt_extent="configured_range",
                  fixed_mt_edges=list(edges))
    return result


class MtBinning(unittest.TestCase):
    def test_quantile_full_output_regression(self):
        # Captured before introducing fixed edges, including every nested audit
        # field: the established quantile-mode output must remain unchanged.
        result = MT.choose_slices(*fixture([0.0] + [6.0] * 8 + [0.0]))
        payload = json.dumps(result, sort_keys=True, separators=(",", ":")).encode()
        self.assertEqual(hashlib.sha256(payload).hexdigest(),
                         "4af390486cf9f3001aaf955331258169afa678d94b2b5b7a23a603ab630ec979")

    def test_quantile_edges_disjoint_and_exact_closure(self):
        data, settings, base = fixture([0.0] + [12.0] * 8 + [0.0])
        result = MT.choose_slices(data, settings, base)
        self.assertEqual(result["boundary_bins"], [1, 3, 5, 7, 9])
        self.assertEqual(result["mt_edges"], [0.05, 0.15, 0.25, 0.35, 0.45])
        self.assertEqual(result["quantile_cost"], 0.0)
        self.assertEqual([s["split_fraction"] for s in result["slices"]], [0.25] * 4)
        covered = [y for s in result["slices"] for y in range(s["first"], s["stop"])]
        self.assertEqual(covered, list(range(1, 9)))
        for mode in MT.MODES:
            for region in MT.REGIONS:
                closure = result["closure"][mode][region]
                self.assertLess(closure["max_abs_count_difference"], 1e-12)
                self.assertLess(closure["max_abs_error2_difference"], 1e-12)
                self.assertEqual(closure["outside_extent"]["below"]["count"], 0.0)
                self.assertEqual(closure["outside_extent"]["above"]["count"], 0.0)
        self.assertEqual(result["slices"][0]["signal_cf_count_audit"],
                         {"display_bins": 3, "zero_se_display": 0,
                          "zero_me_display": 0, "positive_me_display": 3})

    def test_all_channel_tails_are_retained_even_above_display(self):
        data, settings, base = fixture([0.0] + [6.0] * 8 + [0.0, 0.0, 0.0])
        set_cell(data, "SE", "leftSB", 19, 0, 2.0)
        set_cell(data, "ME", "rightSB", 19, 11, 3.0)
        result = MT.choose_slices(data, settings, base)
        self.assertEqual(result["boundary_bins"], [0, 3, 5, 7, 12])
        self.assertEqual(result["slices"][0]["channels"]["SE"]["leftSB"]["rebinned"]["counts"][-1], 2.0)
        self.assertEqual(result["slices"][-1]["channels"]["ME"]["rightSB"]["rebinned"]["counts"][-1], 3.0)

    def test_quantile_ties_choose_lexicographically_lower_edges(self):
        data, settings, base = fixture([6.0] * 6)
        result = MT.choose_slices(data, settings, base)
        self.assertEqual(result["boundary_bins"], [0, 1, 3, 4, 6])
        self.assertAlmostEqual(result["quantile_cost"], 2.0 * (1.0 / 12.0) ** 2)

    def test_dynamic_program_matches_brute_force_with_positive_norm_constraint(self):
        weights, norm = [3.0, 1.0, 5.0, 2.0, 8.0, 4.0, 6.0, 1.0], [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]
        data, settings, base = fixture(weights, norm)
        result = MT.choose_slices(data, settings, base)
        valid = []
        for cuts in itertools.combinations(range(1, len(weights)), 3):
            boundaries = (0,) + cuts + (len(weights),)
            if all(sum(weights[a:b]) > 0.0 and sum(norm[a:b]) > 0.0
                   for a, b in zip(boundaries, boundaries[1:])):
                cost = sum((sum(weights[:edge]) / sum(weights) - (i + 1) / 4.0) ** 2
                           for i, edge in enumerate(cuts))
                valid.append((cost, list(boundaries)))
        expected = min(valid)
        self.assertEqual(result["boundary_bins"], expected[1])
        self.assertAlmostEqual(result["quantile_cost"], expected[0])
        for entry in result["slices"]:
            for mode in MT.MODES:
                self.assertGreater(entry["normalization_signal"][mode]["count"], 0.0)

    def test_zero_total_norm_is_rejected(self):
        data, settings, base = fixture([6.0] * 8, [0.0] * 8)
        with self.assertRaisesRegex(ValueError, "normalization integral"):
            MT.choose_slices(data, settings, base)

    def test_nonzero_total_norm_but_impossible_slice_norm_is_rejected(self):
        data, settings, base = fixture([6.0] * 8, [1.0] + [0.0] * 7)
        with self.assertRaisesRegex(ValueError, "no 4 nonoverlapping"):
            MT.choose_slices(data, settings, base)

    def test_relevant_mt_flow_rejected_with_preserved_audit(self):
        data, settings, base = fixture([6.0] * 8)
        data["channels"]["ME"]["leftSB"]["flowCounts"][1][0] = 2.0
        data["channels"]["ME"]["leftSB"]["flowError2"][1][0] = 3.0
        with self.assertRaises(MT.MtFlowError) as caught:
            MT.choose_slices(data, settings, base)
        audit = caught.exception.audit
        self.assertTrue(audit["relevant_mt_flow_nonzero"])
        self.assertEqual(audit["channels"]["ME"]["leftSB"]["split"]["mt_underflow"],
                         {"count": 2.0, "error2": 3.0})

    def test_norm_only_mt_flow_and_variance_only_flow_rejected(self):
        data, settings, base = fixture([6.0] * 8)
        data["channels"]["SE"]["signal"]["flowError2"][7][-1] = 2.0
        with self.assertRaises(MT.MtFlowError):
            MT.choose_slices(data, settings, base)

    def test_outside_relevant_flow_reported_without_reassignment(self):
        data, settings, base = fixture([6.0] * 8)
        data["channels"]["SE"]["signal"]["flowCounts"][20][0] = 4.0
        data["channels"]["SE"]["signal"]["flowError2"][20][0] = 4.0
        data["channels"]["ME"]["signal"]["flowCounts"][0][2] = 7.0
        result = MT.choose_slices(data, settings, base)
        self.assertEqual(result["boundary_bins"], [0, 2, 4, 6, 8])
        self.assertEqual(result["split_signal_total"], 48.0)
        self.assertTrue(result["flow"]["any_flow_nonzero"])
        self.assertFalse(result["flow"]["relevant_mt_flow_nonzero"])
        self.assertEqual(len(result["warnings"]), 1)

    def test_empty_individual_cf_bin_is_audit_not_new_cut(self):
        data, settings, base = fixture([6.0] * 8)
        for x in (0, 1):
            for y in range(2):
                set_cell(data, "ME", "signal", x, y, 0.0)
        result = MT.choose_slices(data, settings, base)
        self.assertEqual(result["boundary_bins"], [0, 2, 4, 6, 8])
        self.assertEqual(result["slices"][0]["channels"]["ME"]["signal"]["zero_stats"]["display"], [0])
        self.assertEqual(result["slices"][0]["signal_cf_count_audit"]["positive_me_display"], 2)

    def test_counts_and_sumw2_rebin_independently(self):
        data, settings, base = fixture([6.0] * 8)
        result = MT.choose_slices(data, settings, base)
        series = result["slices"][0]["channels"]["SE"]["signal"]["rebinned"]
        self.assertEqual(series["counts"][:3], [4.0, 4.0, 4.0])
        self.assertEqual(series["error2"][:3], [5.0, 5.0, 5.0])
        self.assertEqual(result["slices"][0]["split_signal"]["SE"], {"count": 12.0, "error2": 15.0})

    def test_input_is_not_mutated(self):
        args = fixture([6.0] * 8)
        original = copy.deepcopy(args)
        MT.choose_slices(*args)
        self.assertEqual(args, original)

    def test_invalid_schema_axis_flow_and_nonfinite_fail_loudly(self):
        for case in ("negative", "nonfinite", "flow_mismatch", "axis", "missing_flow", "range"):
            data, settings, base = fixture([6.0] * 8)
            channel = data["channels"]["SE"]["signal"]
            if case == "negative":
                set_cell(data, "SE", "signal", 0, 0, -1.0)
            elif case == "nonfinite":
                set_cell(data, "SE", "signal", 0, 0, math.inf)
            elif case == "flow_mismatch":
                channel["flowCounts"][1][1] += 1.0
            elif case == "axis":
                data["channels"]["ME"]["leftSB"]["mtEdges"][-1] += 0.01
            elif case == "missing_flow":
                del channel["flowCounts"]
            elif case == "range":
                base["normalization"] = [0.31, 0.4]
            with self.subTest(case=case), self.assertRaises((ValueError, KeyError)):
                MT.choose_slices(data, settings, base)

    def test_fixed_saved_edges_and_positive_counts(self):
        data, settings, base = fixture([6.0, 6.0, 12.0, 12.0, 18.0, 18.0, 24.0, 24.0])
        settings = fixed_settings(settings, [0.0, 0.1, 0.2, 0.3, 0.4])
        result = MT.choose_slices(data, settings, base)
        self.assertEqual(result["boundary_bins"], [0, 2, 4, 6, 8])
        self.assertEqual(result["mt_edges"], settings["fixed_mt_edges"])
        self.assertEqual([s["split_signal"]["SE"]["count"] for s in result["slices"]], [12, 24, 36, 48])
        self.assertIsNone(result["quantile_cost"])
        self.assertEqual(result["target_cumulative_fractions"], [])
        self.assertEqual(result["closure_scope"], "configured_range")
        self.assertEqual(result["selection"]["positivity_constraints"], [])
        for item in result["slices"]:
            status = item["count_status"]
            self.assertTrue(status["signal_cf_normalization_valid"])
            self.assertTrue(status["sideband_cf_normalization_valid"])
            self.assertTrue(status["split_signal_se_positive"])
            self.assertFalse(status["empty_all_regular_channels"])

    def test_fixed_edges_do_not_move_for_empty_intervals(self):
        data, settings, base = fixture([6.0, 6.0, 0.0, 0.0, 12.0, 12.0, 18.0, 18.0])
        settings = fixed_settings(settings, [0.0, 0.1, 0.2, 0.3, 0.4])
        result = MT.choose_slices(data, settings, base)
        self.assertEqual(result["boundary_bins"], [0, 2, 4, 6, 8])
        empty = result["slices"][1]
        self.assertEqual(empty["split_fraction"], 0.0)
        self.assertFalse(empty["count_status"]["signal_cf_normalization_valid"])
        self.assertFalse(empty["count_status"]["sideband_cf_normalization_valid"])
        self.assertTrue(empty["count_status"]["empty_all_regular_channels"])
        self.assertTrue(any("Fixed slice 1" in item for item in result["warnings"]))

    def test_fixed_all_empty_remains_auditable_without_division_by_zero(self):
        data, settings, base = fixture([0.0] * 8)
        settings = fixed_settings(settings, [0.0, 0.1, 0.2, 0.3, 0.4])
        result = MT.choose_slices(data, settings, base)
        self.assertEqual(result["split_signal_total"], 0.0)
        self.assertIsNone(result["occupied_extent"])
        self.assertEqual(len(result["slices"]), 4)
        for entry in result["slices"]:
            self.assertIsNone(entry["split_fraction"])
            self.assertIsNone(entry["cumulative_split_fraction"])
            self.assertTrue(entry["count_status"]["empty_all_regular_channels"])
        json.dumps(result, allow_nan=False)

    def test_fixed_zero_normalization_is_audited_not_edge_optimized(self):
        data, settings, base = fixture([6.0] * 8, [0.0] * 8)
        settings = fixed_settings(settings, [0.0, 0.1, 0.2, 0.3, 0.4])
        result = MT.choose_slices(data, settings, base)
        self.assertEqual(result["boundary_bins"], [0, 2, 4, 6, 8])
        for entry in result["slices"]:
            self.assertTrue(entry["count_status"]["split_signal_se_positive"])
            self.assertFalse(entry["count_status"]["signal_cf_normalization_valid"])
            self.assertEqual(entry["count_status"]["normalization_signal_positive"], {"SE": False, "ME": False})

    def test_fixed_excluded_tails_and_configured_range_closure(self):
        data, settings, base = fixture([6.0] * 12)
        settings = fixed_settings(settings, [0.1, 0.2, 0.3, 0.4, 0.5])
        set_cell(data, "ME", "rightSB", 19, 11, 17.0, 23.0)
        result = MT.choose_slices(data, settings, base)
        self.assertEqual(result["boundary_bins"], [2, 4, 6, 8, 10])
        self.assertEqual(result["split_signal_total"], 48.0)
        self.assertEqual(result["occupied_extent"]["first"], 0)
        self.assertEqual(result["occupied_extent"]["stop"], 12)
        tails = result["excluded_tails"]
        self.assertEqual(tails["below"]["SE"]["signal"]["split"]["count"], 12.0)
        self.assertEqual(tails["above"]["SE"]["signal"]["split"]["count"], 12.0)
        self.assertEqual(tails["above"]["ME"]["rightSB"]["rebinned"]["counts"][-1], 17.0)
        self.assertEqual(tails["above"]["ME"]["rightSB"]["rebinned"]["error2"][-1], 23.0)
        for mode in MT.MODES:
            for region in MT.REGIONS:
                self.assertLess(result["closure"][mode][region]["max_abs_count_difference"], 1e-12)
                for name, key in (("count", "counts"), ("error2", "error2")):
                    inside = result["totals"][mode][region]["all_regular_kstar"][name]
                    outside = sum(tails[side][mode][region]["all_regular_kstar"][name]
                                  for side in ("below", "above"))
                    original = sum(sum(row) for row in data["channels"][mode][region][key])
                    self.assertAlmostEqual(inside + outside, original)
        self.assertTrue(any("Populated regular mT tails" in item for item in result["warnings"]))

    def test_fixed_unaligned_out_of_range_or_duplicate_edges_are_rejected(self):
        cases = ([0.0, 0.11, 0.2, 0.3, 0.4],
                 [-0.1, 0.1, 0.2, 0.3, 0.4],
                 [0.0, 0.1, 0.2, 0.3, 0.5],
                 [0.0, 0.1, 0.1, 0.3, 0.4],
                 [0.0, 0.2, 0.1, 0.3, 0.4],
                 [0.0, 0.1, 0.2, 0.4])
        for edges in cases:
            data, settings, base = fixture([6.0] * 8)
            with self.subTest(edges=edges), self.assertRaises(ValueError):
                MT.choose_slices(data, fixed_settings(settings, edges), base)

    def test_fixed_mismatched_extent_and_missing_edges_fail(self):
        data, settings, base = fixture([6.0] * 8)
        settings["split_strategy"] = "fixed_saved_mt_edges"
        with self.assertRaisesRegex(ValueError, "configured_range"):
            MT.choose_slices(data, settings, base)
        settings["mt_extent"] = "configured_range"
        with self.assertRaises(KeyError):
            MT.choose_slices(data, settings, base)

    def test_fixed_relevant_flow_still_requires_explicit_resolution(self):
        data, settings, base = fixture([6.0] * 8)
        settings = fixed_settings(settings, [0.0, 0.1, 0.2, 0.3, 0.4])
        data["channels"]["SE"]["signal"]["flowCounts"][1][0] = 1.0
        with self.assertRaises(MT.MtFlowError):
            MT.choose_slices(data, settings, base)


if __name__ == "__main__":
    unittest.main()
