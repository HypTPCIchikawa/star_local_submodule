#!/usr/bin/env python3
"""Downstream macro contracts; no ROOT runtime or production-file writes."""
import copy
import importlib.util
import json
from pathlib import Path
import unittest

import yaml

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    'kf_prepare', ROOT / 'tools/kfparticle_summary/prepare_macros.py')
PREPARE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PREPARE)


class KFSummaryReporting(unittest.TestCase):
    def test_relative_header_provenance_is_idempotent(self):
        old = dict(support=dict(source='tools/kfparticle_summary/KFSummarySupport.h', sha256='abc'), purity=.97)
        original = (json.dumps(old, indent=2) + '\n').encode()
        canonical = copy.deepcopy(old)
        canonical['support']['source'] = str(ROOT / old['support']['source'])
        self.assertEqual(PREPARE.stable_manifest_bytes(canonical, original), original)
        changed = copy.deepcopy(canonical)
        changed['purity'] = .96
        self.assertNotEqual(PREPARE.stable_manifest_bytes(changed, original), original)
        changed = copy.deepcopy(canonical)
        changed['support']['sha256'] = 'def'
        self.assertNotEqual(PREPARE.stable_manifest_bytes(changed, original), original)

    def test_both_energies_adapt_without_changing_sources(self):
        settings = yaml.safe_load((ROOT / 'tools/kfparticle_summary/summary_settings.yaml').read_text())
        maker = yaml.safe_load((ROOT / 'config/maker/maker_anaFemtoLambda_4He_KFParticle.yaml').read_text())
        window = dict(mean=maker['lambdaSignalMean'], sigma=maker['lambdaSignalSigma'],
                      nsigma=maker['lambdaSignalNSigma'], sideband_outer_factor=maker['lambdaSidebandOuterFactor'])
        window['sideband_scale'] = 1. / (window['sideband_outer_factor']-1.)
        for energy in ('13p5', '3p85'):
            base = ROOT / 'tools'
            if energy == '3p85':
                base /= 'auau3p85GeV'
            for name, pdf, pages in PREPARE.TASKS:
                with self.subTest(energy=energy, macro=name):
                    source = (base / name).read_bytes()
                    adapted, changes, audit = PREPARE.adapt(source.decode(), name,
                        Path('/tmp/input.root'), Path('/tmp/kf-summary-test/pdf'), .97, window, settings)
                    self.assertEqual((base / name).read_bytes(), source)
                    self.assertIn('TString inputFile = "/tmp/input.root"', adapted)
                    self.assertNotIn('fIn->Get("nuclearid")', adapted)
                    self.assertNotIn('TH2D', adapted)
                    self.assertNotIn('nSpecies = 4', adapted)
                    if name.startswith('Fit_'):
                        self.assertIn('KFReport::DrawLLWarnings(fit1,fit2)', adapted)
                        self.assertIn('KFReport::FitLL', adapted)
                    if 'method' in name:
                        self.assertIn('KFReport::Purity', adapted)
                        self.assertIn('nt > 0 && nm > 0', adapted)
                    self.assertGreater(pages, 0)
                    if name == 'Draw_LambdaInvMass.C':
                        self.assertEqual(pdf, 'Lambda_InvMass_legacy.pdf')

    def test_changed_normalization_reaches_every_cf_macro(self):
        settings = yaml.safe_load((ROOT / 'tools/kfparticle_summary/summary_settings.yaml').read_text())
        settings['normalization'] = [.3,.4]
        maker = yaml.safe_load((ROOT / 'config/maker/maker_anaFemtoLambda_4He_KFParticle.yaml').read_text())
        window = dict(mean=maker['lambdaSignalMean'],sigma=maker['lambdaSignalSigma'],
                      nsigma=maker['lambdaSignalNSigma'],sideband_outer_factor=maker['lambdaSidebandOuterFactor'])
        window['sideband_scale'] = 1./(window['sideband_outer_factor']-1.)
        for energy in ('13p5','3p85'):
            source_dir=ROOT/'tools'
            if energy=='3p85':
                source_dir/='auau3p85GeV'
            for name, _, _ in PREPARE.TASKS:
                if name=='Draw_LambdaInvMass.C':
                    continue
                with self.subTest(energy=energy,macro=name):
                    adapted, _, _=PREPARE.adapt((source_dir/name).read_text(),name,
                        Path('/tmp/input.root'),Path('/tmp/kf-summary-test/pdf'),.97,window,settings)
                    self.assertNotIn('0.40001',adapted)
                    self.assertNotIn('0.59999',adapted)
                    self.assertIn(PREPARE.num(.3+1e-5),adapted)
                    self.assertIn(PREPARE.num(.4-1e-5),adapted)


if __name__ == '__main__':
    unittest.main()
