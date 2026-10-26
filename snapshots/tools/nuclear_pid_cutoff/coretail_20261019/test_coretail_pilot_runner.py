"""Small no-ROOT checks for study-only provenance/quoting/settings helpers."""
import importlib.util
import math
from pathlib import Path
import tempfile
import unittest

import yaml

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location('coretail_pilot_runner', HERE / 'run_coretail_pilot.py')
PILOT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PILOT)


class PilotHelpers(unittest.TestCase):
    def test_relative_close_and_nonfinite(self):
        self.assertTrue(PILOT.close(1e6, 1e6 + .001, 1e-8))
        self.assertFalse(PILOT.close(1., 1.1, 1e-8))
        self.assertFalse(PILOT.close(math.nan, math.nan, 1e-8))

    def test_root_path_quoting(self):
        self.assertEqual(PILOT.quote('/tmp/pilot file.root'), '"/tmp/pilot file.root"')
        for invalid in ('bad"path', 'bad\\path', 'bad\npath'):
            with self.assertRaises(ValueError):
                PILOT.quote(invalid)

    def test_identifier_is_not_a_path(self):
        self.assertEqual(PILOT.identifier('he4', 'target_only', 'gaussian_control'), 'he4_target_only_gaussian_control')
        for invalid in ('../escape', '/tmp/escape', 'with space', 'semicolon;exit', ''):
            with self.assertRaises(ValueError):
                PILOT.identifier('he4', invalid)

    def test_selection_preserves_old_windows_and_rejects_invented_bins(self):
        original = [dict(index=1, r=.895, lo=.89, hi=.90, fitLo=35., fitHi=80.),
                    dict(index=2, r=.905, lo=.90, hi=.91, fitLo=34., fitHi=79.)]
        selected = PILOT.select_slices(original, [.905], 1e-9)
        self.assertIs(selected[0], original[1])
        self.assertEqual((selected[0]['fitLo'], selected[0]['fitHi']), (34., 79.))
        with self.assertRaises(ValueError):
            PILOT.select_slices(original, [.90], 1e-9)
        with self.assertRaises(ValueError):
            PILOT.select_slices(original, [.895, .895], 1e-9)

    def test_no_output_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / 'record.json'
            PILOT.save(destination, {'productionEnabled': False})
            self.assertIs(PILOT.read(destination)['productionEnabled'], False)
            with self.assertRaises(FileExistsError):
                PILOT.save(destination, {})

    def test_user_superseded_pilot_is_not_executed(self):
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory)
            PILOT.save(destination / 'superseded.json', {'fitsExecuted': False})
            with self.assertRaisesRegex(ValueError, 'superseded by the user'):
                PILOT.run(destination, 'pending')
            self.assertFalse((destination / 'run.log').exists())

    def test_prepared_input_mutation_is_detected(self):
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / 'settings.tsv'
            PILOT.table(destination, (), [('tailFraction', .1)])
            plan = {'immutableInputs': {str(destination): PILOT.sha(destination)}}
            PILOT.unchanged(plan)
            with destination.open('a') as stream:
                stream.write('unexpected\t1\n')
            with self.assertRaises(ValueError):
                PILOT.unchanged(plan)

    def test_tsv_preserves_shapes_and_duplicate_settings_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            tsv = Path(directory) / 'data.tsv'
            PILOT.table(tsv, ('index', 'mean'), [(41, 1.2345678901234567)])
            self.assertEqual(float(PILOT.rows(tsv)[0]['mean']), 1.2345678901234567)
            duplicate = Path(directory) / 'settings.tsv'
            PILOT.table(duplicate, (), [('key', 1), ('key', 2)])
            with self.assertRaises(ValueError):
                PILOT.settings_tsv(duplicate)

    def test_explicit_small_scope_and_control(self):
        settings = yaml.safe_load((HERE / 'pilot_settings.yaml').read_text())
        self.assertEqual(sum(map(len, settings['representativeRigidity'].values())), 4)
        self.assertEqual(settings['profiles']['gaussian_control']['fraction'], 0)
        self.assertEqual(len(settings['profiles']), 4)
        self.assertEqual(set(settings['backgrounds']), {'linear', 'exponential'})
        self.assertEqual(set(settings['freePolicies']), {'target_only', 'local_neighbors'})
        self.assertEqual(settings['tailSpecies']['he4'], ['he3', 'he4'])

    def test_old_window_and_covariance_guard_source_contract(self):
        source = (HERE / 'FitCoreTailPilot.C').read_text()
        self.assertIn('CT::ComponentIntegral', source)
        self.assertNotIn('Integrate<GaussianEngine>', source)
        self.assertIn('fit.conditionalNdf==int(req.bins.size())-int(fit.shapeNames.size())-int(activeYields)', source)
        self.assertIn('if(fit.shapeCovarianceValid)NuclearFit::WriteMatrix', source)
        self.assertIn('Poisson native count/error/flow', source)


if __name__ == '__main__':
    unittest.main()
