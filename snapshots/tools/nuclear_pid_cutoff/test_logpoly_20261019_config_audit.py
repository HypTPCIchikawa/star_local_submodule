"""Configuration/legacy-polynomial contracts; no ROOT event or production changes."""
import copy
from pathlib import Path
import unittest

import logpoly_20261019_config_audit as audit


class LogpolyStudyContracts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.project = Path(__file__).resolve().parents[2]
        cls.old = audit.load(cls.project/'config'/audit.OLD_PROFILE)
        cls.new = audit.load(cls.project/'config'/audit.NEW_PROFILE)
        cls.header = (cls.project/'StMaker/common/NuclearIdDeDxVsMom.h').read_text()

    def test_live_rollback_study_and_qa_refs(self):
        report, sources, profile = audit.audit_configs(self.project)
        self.assertEqual(len(report['mainconfs']), 5)
        self.assertFalse(report['productionSwitched'])
        self.assertFalse(profile['dedxCalibration']['enabled'])
        self.assertTrue(sources)

    def test_only_authorized_study_profile_changes(self):
        self.assertEqual(audit.audit_profile(self.old, self.new),
                         ['dedxCalibration', 'speciesPOverQVeto'])

    def test_accidental_extra_physics_change_rejected(self):
        doc = copy.deepcopy(self.new)
        doc['maxNSigmaNuclear'] += .1
        with self.assertRaises(ValueError):
            audit.audit_profile(self.old, doc)

    def test_not_silently_keep_old_veto_holes(self):
        doc = copy.deepcopy(self.new)
        doc['speciesPOverQVeto']['d'] = [{'min': 1.45, 'max': 1.5}]
        with self.assertRaises(ValueError):
            audit.audit_profile(self.old, doc)

    def test_calibration_table_cannot_remain_active(self):
        doc = copy.deepcopy(self.new)
        doc['dedxCalibration'] = copy.deepcopy(self.old['dedxCalibration'])
        with self.assertRaises(ValueError):
            audit.audit_profile(self.old, doc)

    def test_original_coefficients_preserved(self):
        tables = audit.parse_header(self.header)
        self.assertEqual(tables['kMean'], [
            (10.2875, -36.6639, 67.63, -45.2483),
            (19.2177, -70.3198, 113.876, -67.0591),
            (28.8747, -68.8637, 104.455, -56.9342),
            (39.974, -103.365, 139.182, -67.1446)])
        self.assertEqual(tables['kSigma'], [
            (.771189, -2.63609, 2.8044, 0.), (1.40996, -4.97147, 4.97608, 0.),
            (1.636, -4.35588, 4.40594, 0.), (2.24127, -7.08328, 8.07124, 0.)])

    def test_variable_cannot_silently_change_to_log(self):
        with self.assertRaises(ValueError):
            audit.parse_header(self.header.replace('TMath::Log10(p)', 'TMath::Log(p)'))

    def test_dense_curve_positive_and_independent_evaluation_agrees(self):
        report = audit.dense_check(audit.parse_header(self.header),
                                   self.new['speciesMinPOverQ'], self.new['maxPOverQ'], 1001)
        self.assertEqual(report['checks'], 4*1001*4)
        self.assertFalse(report['runtimeProviderVerified'])
        self.assertTrue(all(row['minimumSigma'] > 0 for row in report['species']))

    def test_no_low_r_extrapolation(self):
        minima = dict(self.new['speciesMinPOverQ'], he4=.1)
        with self.assertRaises(ValueError):
            audit.dense_check(audit.parse_header(self.header), minima, 5., 10)


if __name__ == '__main__':
    unittest.main()
