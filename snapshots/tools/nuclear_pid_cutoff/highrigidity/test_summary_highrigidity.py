import unittest
from unittest.mock import patch
import summarize_highrigidity as summary
import test_cutoff_summary as old_fixture


class HighSummaryTests(unittest.TestCase):
    def fixtures(self):
        fixture=old_fixture.CutoffSummaryTests();cfg=fixture.settings()
        cfg.update(nativeWidth=.01,aggregateWidth=.02)
        cfg['support']['minDropDeltaD']=25.
        rows=[fixture.row(0,900.,1000.),fixture.row(1,10.,100.)]
        for r in rows:
            r['components']['d']['conditionalDropDeltaD']=100.
            r['components']['d']['selection']['omittedLocalYieldUpperBound']=0.
        audit={i:dict(failures={s:[] for s in ('d','t','he3','he4')}) for i in (0,1)}
        return rows,audit,cfg

    def test_sum_signal_total_not_mean_of_purity(self):
        rows,audit,cfg=self.fixtures();value=summary.grouped_model(rows,'d',audit,cfg)
        self.assertAlmostEqual(value['purity'],910./1100.)
        self.assertTrue(value['supported'])

    def test_one_unassessed_native_cannot_be_hidden_by_aggregation(self):
        rows,audit,cfg=self.fixtures()
        rows[1]['components']['d']['selection']['gateAbsoluteResidualCount']=30.
        value=summary.grouped_model(rows,'d',audit,cfg)
        self.assertLess(value['localResidualFrac'],cfg['cutoff']['maxGateAbsoluteResidualFraction'])
        self.assertFalse(value['supported'])
        self.assertIn('one_or_more_native_fits_unassessed',value['reasons'])

    def test_missing_native_coverage_not_called_complete(self):
        rows,audit,cfg=self.fixtures();value=summary.grouped_model(rows[:1],'d',audit,cfg)
        self.assertFalse(value['supported'])
        self.assertIn('incomplete_native_coverage',value['reasons'])


if __name__=='__main__': unittest.main()
