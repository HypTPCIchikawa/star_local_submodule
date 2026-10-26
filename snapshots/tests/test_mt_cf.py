import copy
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools/kfparticle_summary'))
from run_mt_cf import calculate, report, select_purity_group, shared_purities, study_identity
from mass_half_statistics import ErrorModel


class MtCFTests(unittest.TestCase):
    def fixture(self):
        edges = [i*.005 for i in range(201)]
        channels = {}
        for sample,scale in [('SE',1.),('ME',20.)]:
            channels[sample] = {}
            for region,region_scale in [('signal',100.),('leftSB',3.),('rightSB',6.)]:
                counts = [[scale*region_scale*(y+1) for y in range(4)] for _ in range(200)]
                channels[sample][region] = dict(kstarEdges=edges,mtEdges=[2.4+.05*y for y in range(5)],
                                               counts=counts,error2=copy.deepcopy(counts))
        group = dict(id='cent0_10',centMin=7,centMax=8,originalKstarEdges=edges,method5Base={},methods={'method5':{}})
        for sample in ('SE','ME'):
            group['method5Base'][sample] = [dict(counts=[50.,50.],error2=[50.,50.]) for _ in range(200)]
            group['methods']['method5'][sample] = [dict(full={'valid':True},background=[500.,500.],
                backgroundCovariance=[[100.,10.],[10.,100.]],windowOriginalFirstBin=1,windowOriginalLastBin=200)
                for _ in range(20)]
        source = dict(groups=[group],settings={'purityRebin':10},
                      inclusive={r:dict(count=c,error2=c) for r,c in
                                 [('low',1000),('high',1000),('full',2000),('leftSB',100),('rightSB',100)]},
                      massPartition={'binned':{r:dict(lower=lo,upper=hi) for r,lo,hi in
                                             [('low',1,2),('high',2,3),('full',1,3),('leftSB',-2,1),('rightSB',3,6)]}})
        selection = {'slices':[dict(first=i,stop=i+1,low=2.4+.05*i,high=2.4+.05*(i+1)) for i in range(4)]}
        return dict(centMin=7,centMax=8,channels=channels),selection,source

    def test_uniform_cf_and_shared_cross_slice_error(self):
        data,selection,source = self.fixture()
        result,tsv = calculate(data,selection,dict(cent9_min=7,cent9_max=8,centrality_label='0-10%'),
            dict(rebin_group=4,normalization=[.3,.4],plot_x_range=[0,.3]),source)
        for part in result['slices']:
            for method in part['methods'].values():
                self.assertEqual(method['valid_plot_bins'],15)
                for point in method['points']:
                    self.assertAlmostEqual(point['value'],1.)
                    self.assertGreater(point['error'],0.)
        self.assertGreater(result['cross_slice_covariance']['method5_common_purity'][0][0][1],0.)
        self.assertEqual(result['cross_slice_covariance']['constant_purity'][0][0][1],0.)
        self.assertEqual(len(tsv.splitlines()),401)

    def test_purity_factors_use_one_shared_window_count(self):
        _,_,source = self.fixture()
        model = ErrorModel()
        purity = shared_purities(model,source['groups'][0],source['settings'])
        self.assertAlmostEqual(purity['SE'][0].value,.95)
        self.assertGreater(model.covariance(purity['SE'][0],purity['SE'][1]),0.)

    def test_zero_sideband_me_masks_only_constant_purity(self):
        data,selection,source = self.fixture()
        for region in ('leftSB','rightSB'):
            for ix in range(4):
                for key in ('counts','error2'):
                    data['channels']['ME'][region][key][ix][3] = 0.
        result,_ = calculate(data,selection,dict(cent9_min=7,cent9_max=8,centrality_label='0-10%'),
            dict(rebin_group=4,normalization=[.3,.4],plot_x_range=[0,.3]),source)
        last = result['slices'][3]['methods']
        self.assertIsNone(last['constant_purity']['points'][0])
        self.assertEqual(last['constant_purity']['valid_plot_bins'],14)
        self.assertIsNotNone(last['method5_common_purity']['points'][0])
        self.assertEqual(last['method5_common_purity']['valid_plot_bins'],15)

    def test_empty_slice_keeps_edges_and_masks_both_methods(self):
        data,selection,source = self.fixture()
        for sample in ('SE','ME'):
            for region in ('signal','leftSB','rightSB'):
                for key in ('counts','error2'):
                    for row in data['channels'][sample][region][key]:
                        row[3] = 0.
        result,_ = calculate(data,selection,dict(cent9_min=7,cent9_max=8,centrality_label='0-10%'),
            dict(rebin_group=4,normalization=[.3,.4],plot_x_range=[0,.3]),source)
        self.assertEqual(result['slices'][3]['low'],selection['slices'][3]['low'])
        for curve in result['slices'][3]['methods'].values():
            self.assertEqual(curve['valid_plot_bins'],0)
            self.assertTrue(all(p is None for p in curve['points']))

    def test_legacy_and_nonseven_centrality_names(self):
        legacy = study_identity(dict(cent9_min=7,cent9_max=8,centrality_label='0-10%'))
        self.assertEqual(legacy['stem'],'cf_mt_cent0_10_constant_vs_method5')
        self.assertEqual(legacy['output_directory'],'mt_cent0_10')
        wider = study_identity(dict(cent9_min=6,cent9_max=8,centrality_label='0–20%'))
        self.assertEqual(wider['label'],'0-20%')
        self.assertEqual(wider['stem'],'cf_mt_cent0_20_constant_vs_method5')
        self.assertEqual(wider['output_directory'],'mt_cent0_20')
        custom = study_identity(dict(cent9_min=6,cent9_max=8,centrality_label='0-20%',report_stem='custom_report'))
        self.assertEqual(custom['stem'],'custom_report')
        with self.assertRaises(ValueError):
            study_identity(dict(cent9_min=6,cent9_max=8,centrality_label='0-20%',report_stem='../outside'))

    def test_purity_group_must_match_centrality_and_label(self):
        data,selection,source = self.fixture()
        settings = dict(cent9_min=6,cent9_max=8,centrality_label='0-20%')
        with self.assertRaisesRegex(ValueError,'Exactly one'):
            select_purity_group(source,settings)
        source['groups'][0].update(centMin=6,id='cent0_20')
        self.assertEqual(select_purity_group(source,settings)['id'],'cent0_20')
        with self.assertRaisesRegex(ValueError,'label disagrees'):
            select_purity_group(source,dict(settings,centrality_label='0-10%'))
        with self.assertRaisesRegex(ValueError,'extraction'):
            calculate(data,selection,settings,dict(rebin_group=4,normalization=[.3,.4],plot_x_range=[0,.3]),source)
        data['centMin'] = 6
        result,_ = calculate(data,selection,settings,dict(rebin_group=4,normalization=[.3,.4],plot_x_range=[0,.3]),source)
        self.assertEqual(result['centrality']['purity_group_id'],'cent0_20')
        source['groups'].append(copy.deepcopy(source['groups'][0]))
        with self.assertRaisesRegex(ValueError,'Exactly one'):
            select_purity_group(source,settings)

    def test_report_labels_and_explicit_stem_follow_nonseven_settings(self):
        data,selection,source = self.fixture()
        data['centMin'] = 6
        source['groups'][0].update(centMin=6,id='cent0_20')
        settings = dict(cent9_min=6,cent9_max=8,centrality_label='0-20%',
                        report_stem='cf_mt_cent0_20_constant_vs_method5',
                        split_strategy='fixed_saved_mt_edges',split_kstar=[0,.3])
        base = dict(rebin_group=4,normalization=[.3,.4],plot_x_range=[0,.3])
        selection['warnings'] = []
        for i,part in enumerate(selection['slices']):
            part['channels'] = {'SE':{'signal':{'rebinned':{'counts':[400.*(i+1)]*50}}}}
        result,_ = calculate(data,selection,settings,base,source)
        with tempfile.TemporaryDirectory() as temporary:
            out = Path(temporary)
            (out/'logs').mkdir()
            with patch('run_mt_cf.pdf_pages',return_value=4), patch('subprocess.run') as latex:
                self.assertEqual(report(out,dict(settings=settings,base=base),data,result),4)
            tex = (out/(settings['report_stem']+'.tex')).read_text()
            readme = (out/'README.md').read_text()
            self.assertIn(r'Centrality 0-20\%',tex)
            self.assertIn(r'centrality-0-20\%',tex)
            self.assertIn('Centrality 0–20%',readme)
            self.assertIn('cent9=6・7・8',readme)
            self.assertIn(settings['report_stem']+'.pdf',readme)
            self.assertNotIn('0-10',tex)
            self.assertEqual(latex.call_args[0][0][-1],settings['report_stem']+'.tex')


if __name__=='__main__':
    unittest.main()
