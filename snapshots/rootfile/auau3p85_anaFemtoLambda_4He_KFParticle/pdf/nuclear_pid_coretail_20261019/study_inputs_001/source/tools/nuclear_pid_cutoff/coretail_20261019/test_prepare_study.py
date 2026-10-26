"""Pure configuration contracts; no fitting, freezer execution or production switch."""
import copy
import importlib.util
import os
from pathlib import Path
import unittest

import yaml


HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('coretail_prepare_test_subject', HERE/'prepare_study.py')
prepare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare)


class CoretailPreparation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.project = HERE.parents[2]
        cls.cuts = yaml.safe_load((cls.project/'config'/prepare.LOW_PROFILE).read_text())
        cls.base = yaml.safe_load((cls.project/'tools/nuclear_pid_cutoff/free_local_settings.yaml').read_text())

    def test_reference_is_explicit_and_freezer_compatible(self):
        prepare.validate_reference(self.cuts)
        selection, paths = prepare.inspect_reference(self.project, self.project/prepare.LOW_MAIN)
        self.assertIsNone(selection['calibration'])
        self.assertEqual(selection['cuts']['maxPOverQ'], 2.)
        self.assertTrue(paths)

    def test_relative_script_filename_has_absolute_snapshot_sources(self):
        relative = Path(os.path.relpath(HERE/'prepare_study.py', Path.cwd()))
        self.assertFalse(relative.is_absolute())
        sources = prepare.own_source_paths(relative)
        self.assertTrue(all(path.is_absolute() for path in sources))
        self.assertEqual(sources, ((HERE/'prepare_study.py').resolve(),
                                   (HERE/'test_prepare_study.py').resolve()))
        for path in sources:
            self.assertTrue(path.relative_to(self.project))

    def test_per_species_min_or_veto_keys_must_be_absent(self):
        for key in ('speciesMinPOverQ', 'speciesPOverQVeto'):
            cuts = dict(self.cuts, **{key: {}})
            with self.assertRaises(ValueError):
                prepare.validate_reference(cuts)

    def test_tof_key_must_be_absent_even_if_disabled(self):
        with self.assertRaises(ValueError):
            prepare.validate_reference(dict(self.cuts, highRigidityTof={'enabled': False}))

    def test_low_reference_cannot_inherit_interpolated_calibration(self):
        with self.assertRaises(ValueError):
            prepare.validate_reference(dict(self.cuts, dedxCalibration={'enabled': True}))

    def test_low_reference_cannot_extend_into_high_r(self):
        with self.assertRaises(ValueError):
            prepare.validate_reference(dict(self.cuts, maxPOverQ=5.))

    def test_settings_change_only_reference_and_output(self):
        result = prepare.derive_settings(self.base, prepare.LOW_MAIN, prepare.OUTPUT_BASE)
        changed = {key for key in result if result[key] != self.base[key]}
        self.assertEqual(changed, {'runtimeMainconf', 'baselineMainconf', 'outputBase'})
        self.assertEqual(result['inputStudy'], self.base['inputStudy'])
        self.assertEqual(result['tagStudies'], self.base['tagStudies'])
        self.assertEqual(result['support'], self.base['support'])
        self.assertNotIn('tailFraction', result.get('numeric', {}))

    def test_derivation_does_not_modify_parent_settings(self):
        before = copy.deepcopy(self.base)
        result = prepare.derive_settings(self.base, prepare.LOW_MAIN, prepare.OUTPUT_BASE)
        result['numeric']['shapeMaxCalls'] += 1
        self.assertEqual(self.base, before)

    def test_no_low_r_extrapolation(self):
        with self.assertRaises(ValueError):
            prepare.derive_settings(dict(self.base, rigidityMin=.2), prepare.LOW_MAIN, prepare.OUTPUT_BASE)

    def test_no_joining_high_r_to_all_cache(self):
        with self.assertRaises(ValueError):
            prepare.derive_settings(dict(self.base, rigidityMax=5.), prepare.LOW_MAIN, prepare.OUTPUT_BASE)

    def profile(self):
        return dict(schemaVersion=1, profiles=[dict(label='explicit_test_only',
                    tailFraction=.1, tailOffsetSeedSigma=-1., tailWidthRatio=2.)])

    def test_explicit_fixed_profile_is_preserved(self):
        doc = self.profile()
        self.assertEqual(prepare.validate_tail_profiles(doc), doc)

    def test_empty_grid_has_no_silent_defaults(self):
        with self.assertRaises(ValueError):
            prepare.validate_tail_profiles(dict(schemaVersion=1, profiles=[]))

    def test_gaussian_zero_fraction_control(self):
        doc = self.profile()
        doc['profiles'][0].update(tailFraction=0., tailWidthRatio=1.)
        self.assertEqual(prepare.validate_tail_profiles(doc), doc)

    def test_nonzero_tail_must_be_wider_than_core(self):
        doc = self.profile()
        doc['profiles'][0]['tailWidthRatio'] = 1.
        with self.assertRaises(ValueError):
            prepare.validate_tail_profiles(doc)

    def test_invalid_fraction_or_nonfinite_parameter(self):
        for key, bad in [('tailFraction', .5), ('tailFraction', -.1),
                         ('tailOffsetSeedSigma', float('nan')), ('tailWidthRatio', True)]:
            doc = self.profile()
            doc['profiles'][0][key] = bad
            with self.assertRaises(ValueError):
                prepare.validate_tail_profiles(doc)

    def test_duplicate_or_unsafe_labels(self):
        doc = self.profile()
        doc['profiles'].append(copy.deepcopy(doc['profiles'][0]))
        with self.assertRaises(ValueError):
            prepare.validate_tail_profiles(doc)
        doc = self.profile()
        doc['profiles'][0]['label'] = '../outside'
        with self.assertRaises(ValueError):
            prepare.validate_tail_profiles(doc)


if __name__ == '__main__':
    unittest.main()
