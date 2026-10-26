import math
import unittest

from render_quality_diagnostics import display_number, flag, focus_extent, quoted, verify_panel_policy


class QualityViewContracts(unittest.TestCase):
    def setUp(self):
        self.cfg = dict(focusSpecies={'he4': ['he3', 'he4']}, focusSigma=3.5,
                        focusMinimumWidth=4., plots=dict(dedxMin=0., dedxMax=120.))
        self.item = dict(numericalValid=True, components={
            'he3': dict(initialMean=20., initialSigma=1., mean=21., sigma=1.),
            'he4': dict(initialMean=30., initialSigma=2., mean=32., sigma=2.)})

    def test_local_range_contains_both_nuclei_and_saved_initial_fit_shapes(self):
        self.assertEqual(focus_extent(self.item, 'he4', self.cfg), (16.5, 39.))

    def test_failed_fit_is_not_used_to_set_focus_range(self):
        self.item['numericalValid'] = False
        self.item['components']['he4']['mean'] = math.nan
        self.assertEqual(focus_extent(self.item, 'he4', self.cfg), (16.5, 37.))

    def test_invalid_supported_fit_shape_rejected(self):
        self.item['components']['he4']['sigma'] = -1.
        with self.assertRaises(ValueError):
            focus_extent(self.item, 'he4', self.cfg)

    def test_view_clamped_to_observed_axis(self):
        self.cfg['plots'] = dict(dedxMin=18., dedxMax=35.)
        self.assertEqual(focus_extent(self.item, 'he4', self.cfg), (18., 35.))

    def test_missing_chi2_is_not_shown_as_zero(self):
        for value in (None, '', 'null', 'None', 'NA', float('nan')):
            self.assertEqual(display_number(value), 'NA')
        self.assertEqual(display_number('2.25'), '2.25')

    def test_flags_are_explicit(self):
        self.assertTrue(flag('True'))
        self.assertFalse(flag('0'))
        with self.assertRaises(ValueError):
            flag('unknown')

    def test_unsafe_root_argument_rejected(self):
        self.assertEqual(quoted('/tmp/input.root'), '"/tmp/input.root"')
        with self.assertRaises(ValueError):
            quoted('/tmp/"bad.root')

    def panel_log(self, linear_focus=1):
        return '\n'.join('SPECTRUM_VIEW_POLICY serial={} log={} nuclearFocus={} linearFocusOnly=1'.format(
            i, i % 2, linear_focus if i % 2 == 0 else int(i >= 4)) for i in range(8))

    def test_linear_panels_are_all_local_but_wide_log_is_kept(self):
        self.assertEqual(len(verify_panel_policy(self.panel_log(), 2, {'linearFocusOnly': 1})), 8)

    def test_full_range_linear_regression_rejected(self):
        with self.assertRaises(ValueError):
            verify_panel_policy(self.panel_log(linear_focus=0), 2, {'linearFocusOnly': 1})

    def test_missing_focus_audit_is_rejected(self):
        with self.assertRaises(ValueError):
            verify_panel_policy(self.panel_log(), 3, {'linearFocusOnly': 1})

    def test_no_silent_linear_focus_default(self):
        with self.assertRaises(KeyError):
            verify_panel_policy(self.panel_log(), 2, {})


if __name__ == '__main__':
    unittest.main()
