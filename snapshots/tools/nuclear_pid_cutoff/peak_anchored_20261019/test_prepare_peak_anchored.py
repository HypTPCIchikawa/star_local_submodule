"""No-ROOT checks: fixed original PID, native regions, and immutable study helpers."""
import copy
import importlib.util
import json
import math
from pathlib import Path
import tempfile
import unittest

import yaml

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location('peak_anchored_prep', HERE / 'prepare_peak_anchored.py')
STUDY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(STUDY)


class PreparationTests(unittest.TestCase):
    def setUp(self):
        self.cfg = yaml.safe_load((HERE / 'study_settings.yaml').read_text())
        self.axis = dict(rRange=[0, 6], rBinWidth=.01, rBins=600,
                         dedxRange=[0, 120], dedxBinWidth=.1)

    def test_four_models_and_fixed_three_peak_contract(self):
        STUDY.validate_config(self.cfg)
        self.assertEqual(len(self.cfg['responseModels']) * len(self.cfg['backgroundModels']), 4)
        self.assertEqual(self.cfg['freeSpecies'], [])
        self.assertEqual(self.cfg['fitSpecies'], ['t', 'he3', 'he4'])
        self.assertEqual(len(self.cfg['representativeRigidity']), 7)

    def test_free_shape_or_wrong_basis_is_rejected(self):
        for key, value in (('freeSpecies', ['t']), ('shapePolicy', 'per_bin_fit'),
                           ('fitSpecies', ['d', 't', 'he3', 'he4'])):
            bad = copy.deepcopy(self.cfg)
            bad[key] = value
            with self.assertRaises(ValueError):
                STUDY.validate_config(bad)

    def test_native_r_centers_only(self):
        self.assertEqual(STUDY.native_bin(.825, self.axis, 1e-9)[0], 83)
        with self.assertRaises(ValueError):
            STUDY.native_bin(.820, self.axis, 1e-9)

    def test_roi_whole_native_bins_and_no_silent_axis_clamp(self):
        self.assertEqual(STUDY.snapped(20.17, 58.76, .5, self.axis, 1e-9), (20., 59.))
        with self.assertRaises(ValueError):
            STUDY.snapped(-.1, 58., .5, self.axis, 1e-9)
        with self.assertRaises(ValueError):
            STUDY.snapped(10.12, 20.34, .15, self.axis, 1e-9)

    def test_original_pid_gl3_and_three_peak_coverage(self):
        reader = STUDY.module(STUDY.PROJECT / self.cfg['runtimeReader'])
        tables = reader.parse_legacy_tables(STUDY.PROJECT / self.cfg['legacyHeader'])
        runtime = dict(tables=tables, calibration=None,
                       cuts=dict(maxPOverQ=2., maxNSigmaNuclear=2.))
        slices, shapes, gates = STUDY.make_slices(self.cfg, runtime, reader, self.axis)
        self.assertEqual((len(slices), len(shapes), len(gates)), (7, 84, 84))
        for item in slices:
            for node, (r, weight) in enumerate(reader.gl3_nodes(item['lo'], item['hi'])):
                row = [x for x in shapes if x[0] == item['index'] and x[4] == node]
                self.assertEqual(len(row), 4)
                for s, shape in enumerate(row):
                    self.assertEqual((shape[8], shape[9]), reader.legacy_parameters(tables, s, r))
                    self.assertGreater(weight, 0.)
                tmean, tsigma = reader.legacy_parameters(tables, 1, r)
                hmean, hsigma = reader.legacy_parameters(tables, 2, r)
                fmean, fsigma = reader.legacy_parameters(tables, 3, r)
                self.assertLessEqual(item['fitLo'], tmean - 2 * tsigma)
                self.assertLessEqual(item['regions']['helium']['lo'], hmean - 2 * hsigma)
                self.assertGreaterEqual(item['fitHi'], fmean + 4 * fsigma)

    def test_duplicate_or_invented_r_slice_is_rejected(self):
        reader = STUDY.module(STUDY.PROJECT / self.cfg['runtimeReader'])
        runtime = dict(tables=reader.parse_legacy_tables(STUDY.PROJECT / self.cfg['legacyHeader']),
                       calibration=None, cuts=dict(maxPOverQ=2.))
        bad = copy.deepcopy(self.cfg)
        bad['representativeRigidity'][1] = bad['representativeRigidity'][0]
        with self.assertRaises(ValueError):
            STUDY.make_slices(bad, runtime, reader, self.axis)

    def test_existing_record_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / 'plan.json'
            STUDY.write(destination, {'productionEnabled': False})
            with self.assertRaises(FileExistsError):
                STUDY.write(destination, {})

    def test_mutated_immutable_input_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / 'original.txt'
            destination.write_text('original')
            plan = dict(immutableInputs={str(destination): STUDY.sha(destination)})
            STUDY.unchanged(plan)
            destination.write_text('changed')
            with self.assertRaises(ValueError):
                STUDY.unchanged(plan)


if __name__ == '__main__':
    unittest.main()
