"""Pure-Python schema checks; do not certify the free response or isotope purity."""
import copy
import importlib.util
from pathlib import Path
import unittest

import yaml


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location('free_coretail_contract_runner', HERE/'run_coretail_pilot.py')
RUNNER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNNER)
NULL_SELECTION_FIELDS = ('purity', 'signalSelected', 'totalSelected', 'backgroundSelected',
                         'conditionalError', 'purityLowerIncludingOmittedBound',
                         'selectedColumnFractions', 'selectedColumnCounts')


def fixture():
    cfg = yaml.safe_load((HERE/'free_pilot_settings.yaml').read_text())
    job = dict(target='he4', background='linear', profile='free_symmetric_seed')
    selection = dict(available=False,
                     reason='group_purity_not_evaluated_in_response_quality_pilot',
                     **{key: None for key in NULL_SELECTION_FIELDS})
    result = dict(responseModel='gaussian_coretail_free_components', responseStudyMode='free_coretail',
                  sampleHistogram='ALL', targetSpecies='he4', backgroundModel='linear',
                  profileName='free_symmetric_seed', fullUncertaintyValidated=False,
                  purityPolicy='Not estimated in this response-quality pilot; null is not zero',
                  slices=[dict(components={name: dict(selection=copy.deepcopy(selection))
                                           for name in RUNNER.SPECIES})])
    return cfg, job, result


class FreeCoreTailContracts(unittest.TestCase):
    def setUp(self):
        self.cfg, self.job, self.result = fixture()

    def validate(self):
        RUNNER.validate_result_identity(self.result, self.cfg, self.job)

    def test_real_settings_declare_free_adapter_and_engine(self):
        self.assertEqual(self.cfg['responseStudyMode'], 'free_coretail')
        self.assertEqual(self.cfg['responseModel'], 'gaussian_coretail_free_components')
        self.assertEqual(Path(self.cfg['adapter']).stem, self.cfg['adapterFunction'])
        self.assertEqual(self.cfg['adapterFunction'], 'FitFreeCoreTailPilot')
        self.assertEqual(Path(self.cfg['engine']).name, 'FreeGaussianCoreTail.h')

    def test_seeds_are_distinct_from_fixed_profiles(self):
        profiles = self.cfg['profiles']
        self.assertEqual(set(profiles), {'gaussian_control', 'free_symmetric_seed'})
        self.assertEqual(profiles['gaussian_control']['fraction'], 0)
        self.assertGreater(profiles['free_symmetric_seed']['fraction'], 0)
        self.assertGreater(profiles['free_symmetric_seed']['widthRatio'], 1)

    def test_tail_species_match_explicit_free_core_policy(self):
        for policy, species in self.cfg['freePolicies'].items():
            self.assertEqual(species, self.cfg['tailPolicies'][policy])

    def test_valid_null_purity_schema_is_accepted(self):
        self.validate()

    def test_unavailable_purity_is_not_zero(self):
        self.result['slices'][0]['components']['he4']['selection']['purity'] = 0
        with self.assertRaises(ValueError):
            self.validate()

    def test_every_unavailable_estimator_field_must_stay_null(self):
        for key in NULL_SELECTION_FIELDS:
            with self.subTest(field=key):
                self.cfg, self.job, self.result = fixture()
                self.result['slices'][0]['components']['he4']['selection'][key] = 0
                with self.assertRaises(ValueError):
                    self.validate()

    def test_boolean_available_true_does_not_certify_a_pilot(self):
        self.result['slices'][0]['components']['he4']['selection']['available'] = True
        with self.assertRaises(ValueError):
            self.validate()

    def test_wrong_mode_and_response_are_rejected(self):
        for key, value in (('responseStudyMode', 'fixed'), ('responseModel', 'gaussian')):
            with self.subTest(field=key):
                self.cfg, self.job, self.result = fixture()
                self.result[key] = value
                with self.assertRaises(ValueError):
                    self.validate()

    def test_uncertainty_validation_is_not_claimed(self):
        self.result['fullUncertaintyValidated'] = True
        with self.assertRaises(ValueError):
            self.validate()

    def test_target_background_profile_and_sample_identity_checked(self):
        for key, value in (('targetSpecies', 'he3'), ('backgroundModel', 'exponential'),
                           ('profileName', 'positive_tail'), ('sampleHistogram', 'TAG_D')):
            with self.subTest(field=key):
                self.cfg, self.job, self.result = fixture()
                self.result[key] = value
                with self.assertRaises(ValueError):
                    self.validate()

    def test_folded_four_isotope_contract_is_complete(self):
        del self.result['slices'][0]['components']['he3']
        with self.assertRaises(ValueError):
            self.validate()

    def test_fixed_mode_keeps_existing_identity_contract(self):
        fixed = dict(responseModel='gaussian_coretail_fixed_profile', sampleHistogram='ALL',
                     targetSpecies='he4', backgroundModel='linear', profileName='symmetric_tail')
        RUNNER.validate_result_identity(fixed, {}, dict(target='he4', background='linear',
                                                       profile='symmetric_tail'))


if __name__ == '__main__':
    unittest.main()
