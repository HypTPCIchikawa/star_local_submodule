"""Input grouping and not-yet-certified quality labels, without ROOT or fits."""
import copy
import importlib.util
from pathlib import Path
import unittest


HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('coretail_render_test_subject', HERE/'render_pilot.py')
render = importlib.util.module_from_spec(spec)
spec.loader.exec_module(render)
helper = render.module('coretail_render_test_helper', render.QUALITY/'render_quality_diagnostics.py')


def fixture_plan():
    cfg = dict(representativeRigidity={'d': [1.475], 't': [1.375], 'he4': [.895, 1.775]},
               backgrounds=['linear', 'exponential'], freePolicies={'target_only': {}, 'local_neighbors': {}},
               profiles={name: {} for name in render.PROFILE_ORDER})
    jobs = []
    for target, rigidity in cfg['representativeRigidity'].items():
        slices = [dict(index=index, r=r, fitLo=10., fitHi=40.) for index, r in enumerate(rigidity)]
        for background in cfg['backgrounds']:
            for policy in cfg['freePolicies']:
                for profile in render.PROFILE_ORDER:
                    jobs.append(dict(id='_'.join((target, background, policy, profile)),
                                     target=target, background=background, freePolicy=policy,
                                     profile=profile, slices=copy.deepcopy(slices), rows=len(slices)))
    return dict(config=cfg, jobs=jobs)


class CoretailRenderInputs(unittest.TestCase):
    def setUp(self):
        self.cfg = dict(focusSpecies=render.FOCUS, focusSigma=3.5, focusMinimumWidth=4.,
                        plots=dict(dedxMin=0., dedxMax=120.))
        self.job = dict(target='he4', profile='symmetric_tail')
        self.item = dict(index=5, numericalValid=True, pearson=12., conditionalNdf=4,
                         fitLo=10., fitHi=40., components={
                             'he3': dict(initialMean=20., initialSigma=1., mean=21., sigma=1.),
                             'he4': dict(initialMean=30., initialSigma=2., mean=32., sigma=2.)})

    def quality(self):
        return dict(zip(render.QUALITY_HEADER,
                        render.quality_row(self.job, self.item, self.cfg, helper)))

    def test_all_48_jobs_become_12_books_64_pages(self):
        groups = render.grouped_jobs(fixture_plan())
        self.assertEqual(len(groups), 12)
        self.assertEqual(sum(group['pages'] for group in groups), 64)
        self.assertEqual(sorted(group['pages'] for group in groups), [4]*8+[8]*4)

    def test_each_book_preserves_four_profile_order(self):
        for group in render.grouped_jobs(fixture_plan()):
            self.assertEqual(tuple(job['profile'] for job in group['jobs']), render.PROFILE_ORDER)

    def test_missing_profile_cannot_silently_omit_failed_fit(self):
        plan = fixture_plan()
        plan['jobs'].pop()
        with self.assertRaises(ValueError):
            render.grouped_jobs(plan)

    def test_duplicate_job_rejected(self):
        plan = fixture_plan()
        plan['jobs'][-1] = copy.deepcopy(plan['jobs'][0])
        with self.assertRaises(ValueError):
            render.grouped_jobs(plan)

    def test_moved_old_window_rejected(self):
        plan = fixture_plan()
        plan['jobs'][0]['slices'][0]['fitLo'] += .1
        with self.assertRaises(ValueError):
            render.grouped_jobs(plan)

    def test_actual_conditional_ndf_used(self):
        row = self.quality()
        self.assertEqual((row['chi2'], row['ndf'], row['chi2Ndf']), ('12', '4', '3'))

    def test_support_is_explicitly_not_evaluated_not_fail_or_pass(self):
        row = self.quality()
        self.assertEqual((row['supported'], row['supportKnown']), (0, 0))
        self.assertEqual(row['reason'], render.PENDING)
        self.assertEqual(row['numerical'], 1)

    def test_failed_fit_remains_with_original_numerical_flag(self):
        self.item['numericalValid'] = False
        self.item['components']['he4']['mean'] = float('nan')
        row = self.quality()
        self.assertEqual(row['numerical'], 0)
        self.assertEqual((row['focusLo'], row['focusHi']), (16.5, 37.))

    def test_nonpositive_ndf_does_not_create_a_fake_ratio(self):
        self.item['conditionalNdf'] = 0
        self.assertEqual(self.quality()['chi2Ndf'], 'NA')

    def test_absent_chi2_stays_na_not_zero(self):
        self.item['pearson'] = None
        row = self.quality()
        self.assertEqual(row['chi2'], 'NA')
        self.assertEqual(row['chi2Ndf'], 'NA')

    def test_he4_focus_contains_both_helium_neighbors(self):
        row = self.quality()
        self.assertEqual((row['focusLo'], row['focusHi']), (16.5, 39.))

    def test_t_focus_does_not_expand_to_unrelated_d_peak(self):
        self.job['target'] = 't'
        self.item['components']['t'] = dict(initialMean=22., initialSigma=1., mean=22., sigma=1.)
        self.item['components']['d'] = dict(initialMean=3., initialSigma=1., mean=3., sigma=1.)
        row = self.quality()
        self.assertEqual((row['focusLo'], row['focusHi']), (18.5, 25.5))

    def test_nonboolean_saved_numerical_flag_rejected(self):
        self.item['numericalValid'] = 'true'
        with self.assertRaises(ValueError):
            self.quality()


if __name__ == '__main__':
    unittest.main()
