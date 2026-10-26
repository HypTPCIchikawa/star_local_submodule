"""No-ROOT presentation contracts; no statistical or PID certification."""
import copy
import importlib.util
from pathlib import Path
import unittest


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location('peak_anchored_view_subject', HERE/'render_peak_anchored.py')
VIEW = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VIEW)


def fixture():
    job = dict(id='gaussian_linear', response='gaussian')
    components = {name: dict(mean=center, sigma=1., initialMean=center, initialSigma=1.,
                             freeMean=False, freeSigma=False)
                  for name, center in (('d', 8.), ('t', 10.), ('he3', 25.), ('he4', 40.))}
    regions = {name: dict(lo=low, hi=high, binCount=50, pearson=8., maxAbsPearsonPull=1.5)
               for name, low, high in (('joint', 5., 50.), ('helium', 20., 50.), ('tControl', 5., 15.))}
    fit = dict(shapeSource='legacy_log10_polynomial', backgroundPeakShapesFree=False)
    row = dict(index=3, numericalValid=True, fitLo=5., fitHi=50., components=components,
               regions=regions, pearson=12., conditionalNdf=4, fitFailure='')
    plots = dict(dedxMin=0., dedxMax=120., displayDedxRebinFactor=1, logMinimum=.1)
    return job, fit, row, plots


def panels():
    output, serial = [], 0
    for role in VIEW.REGIONS:
        for pair in range(2):
            for model in range(2):
                for logarithmic in (0, 1):
                    output.append('PEAK_ANCHORED_PANEL serial=%d index=3 label=%s role=%s log=%d low=5 high=50'
                                  % (serial, 'abcd'[2*pair+model], role, logarithmic))
                    output.append('SPECTRUM_VISIBLE_AXES serial=%d' % serial)
                    serial += 1
    return '\n'.join(output)


class PeakAnchoredViewContracts(unittest.TestCase):
    def setUp(self):
        self.job, self.fit, self.row, self.plots = fixture()

    def metric(self):
        return dict(zip(VIEW.METRIC_HEADER, VIEW.metric_row(self.job, self.fit, self.row, self.plots)))

    def test_fixed_original_pid_shapes_and_actual_region_geometry(self):
        result = self.metric()
        self.assertEqual((result['muT'], result['muHe3'], result['muHe4']), (10., 25., 40.))
        self.assertEqual((result['heLo'], result['heHi'], result['tLo'], result['tHi']), (20., 50., 5., 15.))

    def test_neither_t_nor_he3_may_float_as_background(self):
        for species in ('t', 'he3'):
            for key in ('freeMean', 'freeSigma'):
                with self.subTest(species=species, key=key):
                    self.job, self.fit, self.row, self.plots = fixture()
                    self.row['components'][species][key] = True
                    with self.assertRaises(ValueError):
                        self.metric()

    def test_drifted_center_or_width_rejected_even_with_false_free_flags(self):
        for species in ('d', 't', 'he3', 'he4'):
            for key in ('mean', 'sigma'):
                with self.subTest(species=species, key=key):
                    self.job, self.fit, self.row, self.plots = fixture()
                    self.row['components'][species][key] += .01
                    with self.assertRaises(ValueError):
                        self.metric()

    def test_unknown_pid_source_or_free_background_peak_rejected(self):
        for key, value in (('shapeSource', 'independent_fit_points'), ('backgroundPeakShapesFree', True)):
            with self.subTest(key=key):
                self.job, self.fit, self.row, self.plots = fixture()
                self.fit[key] = value
                with self.assertRaises(ValueError):
                    self.metric()

    def test_failed_fit_is_kept_not_certified_or_omitted(self):
        self.row.update(numericalValid=False, pearson=None, conditionalNdf=0,
                        fitFailure='numerical failure\tretained\nfor diagnosis')
        result = self.metric()
        self.assertEqual((result['numerical'], result['jointChi2'], result['jointRatio']), (0, 'NA', 'NA'))
        self.assertNotIn('\t', result['reason'])
        self.assertNotIn('\n', result['reason'])

    def test_only_actual_joint_ndf_and_ratio_are_used(self):
        result = self.metric()
        self.assertEqual((result['jointChi2'], result['jointNdf'], result['jointRatio']), ('12', '4', '3'))
        self.assertEqual((result['heBins'], result['tBins']), (50, 50))
        self.assertFalse(any('ndf' in key.lower() for key in result if key.startswith(('he', 't'))))

    def test_zero_joint_ndf_has_no_invented_ratio(self):
        self.row['conditionalNdf'] = 0
        self.assertEqual(self.metric()['jointRatio'], 'NA')

    def test_full_axis_linear_nuclear_view_is_rejected(self):
        self.row['fitLo'] = self.row['regions']['joint']['lo'] = 0.
        self.row['fitHi'] = self.row['regions']['joint']['hi'] = 120.
        self.row['regions']['helium'].update(lo=0., hi=120.)
        with self.assertRaises(ValueError):
            self.metric()

    def test_region_outside_joint_fit_is_rejected(self):
        self.row['regions']['tControl']['lo'] = 4.
        with self.assertRaises(ValueError):
            self.metric()

    def test_region_binning_must_be_saved_and_positive(self):
        self.row['regions']['helium']['binCount'] = 0
        with self.assertRaises(ValueError):
            self.metric()

    def test_every_model_roi_and_scale_has_an_axis_audit(self):
        records = VIEW.verify_panels(panels(), [(3, .875)], [[name, name+'.root'] for name in 'abcd'], 4)
        self.assertEqual(len(records), 16)

    def test_missing_or_repeated_panel_rejected(self):
        for text in (panels().replace('PEAK_ANCHORED_PANEL serial=0', 'LOST_PANEL serial=0'),
                     panels().replace('serial=15 index=', 'serial=14 index=')):
            with self.subTest(text=text[:80]):
                with self.assertRaises(ValueError):
                    VIEW.verify_panels(text, [(3, .875)], [[name, name+'.root'] for name in 'abcd'], 4)

    def test_duplicating_one_model_cannot_hide_a_missing_model(self):
        text = panels().replace('label=d ', 'label=c ')
        with self.assertRaises(ValueError):
            VIEW.verify_panels(text, [(3, .875)], [[name, name+'.root'] for name in 'abcd'], 4)

    def test_requested_center_roundoff_keeps_exact_native_bin_identity(self):
        planned = dict(index=6, rBin=183, lo=1.82, hi=1.83, r=1.825)
        saved = dict(planned, r=.5*(planned['lo']+planned['hi']))
        self.assertNotEqual(saved['r'], planned['r'])
        self.assertEqual(VIEW.native_selection([saved], [planned], 1.e-9), [(6, saved['r'])])

    def test_wrong_native_bin_or_edge_cannot_hide_inside_rigidity_tolerance(self):
        planned = dict(index=6, rBin=183, lo=1.82, hi=1.83, r=1.825)
        original = dict(planned, r=.5*(planned['lo']+planned['hi']))
        for key, value in (('index', 5), ('rBin', 184), ('lo', 1.8200000000000003), ('hi', 1.84)):
            with self.subTest(field=key):
                saved = dict(original)
                saved[key] = value
                with self.assertRaises(ValueError):
                    VIEW.native_selection([saved], [planned], 1.e-9)

    def test_material_requested_rigidity_mismatch_is_rejected(self):
        planned = dict(index=6, rBin=183, lo=1.82, hi=1.83, r=1.826)
        saved = dict(planned, r=.5*(planned['lo']+planned['hi']))
        with self.assertRaises(ValueError):
            VIEW.native_selection([saved], [planned], 1.e-9)

    def test_saved_rigidity_must_be_the_actual_midpoint_not_a_rounded_label(self):
        planned = dict(index=6, rBin=183, lo=1.82, hi=1.83, r=1.825)
        with self.assertRaises(ValueError):
            VIEW.native_selection([planned], [planned], 1.e-9)


if __name__ == '__main__':
    unittest.main()
