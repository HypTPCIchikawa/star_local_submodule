#!/usr/bin/env python3
"""Check additive He3/mapped histogram contracts against immutable templates."""
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    return yaml.safe_load((ROOT / name).read_text())


def main():
    maker = load('config/maker/maker_auau3p85_anaFemtoLambda_4He_KFParticle_he3Control.yaml')
    assert maker['lambdaHe3ControlEnabled'] is True
    assert maker['speciesKeys'] == 'lambda,he4,he3' and maker['nChannels'] == 8
    for i in range(8):
        nucleus = 'he4' if i < 4 else 'he3'
        suffix = ['', '_signal', '_leftSB', '_rightSB'][i % 4]
        assert maker[f'channel_{i}_name'] == 'lambda_' + nucleus + suffix
        assert maker[f'channel_{i}_partA'] == 'lambda'
        assert maker[f'channel_{i}_partB'] == nucleus
        assert maker[f'channel_{i}_enabled'] and maker[f'channel_{i}_doMixing']
        if i >= 4:
            for key in ['signalMin', 'signalMax', 'normQMin', 'normQMax']:
                assert maker[f'channel_{i}_{key}'] == maker[f'channel_{i-4}_{key}']
    top = load('config/hist/hist_auau3p85_femtoLambda_he4_he3_control_kf.yaml')['histograms']
    nuclear = load('config/hist/hist_auau3p85_femtoLambda_he4_he3_control_nuclear.yaml')['histograms']
    old_top = load('config/hist/hist_femtoLambda_he4_kf.yaml')['histograms']
    old_nuclear = load('config/hist/hist_femtoLambda_he4_nuclear_tofExtended.yaml')['histograms']
    assert len(top) == 129 and len(nuclear) == 561 and len(top)+len(nuclear) == 690
    assert all(top[k] == v for k, v in old_top.items() if k != 'hMixSamplerQA')
    qa = top['hMixSamplerQA']
    assert {k: v for k, v in qa.items() if k != 'yAxis'} == {
        k: v for k, v in old_top['hMixSamplerQA'].items() if k != 'yAxis'}
    assert qa['yAxis'] == {'nBins': 8, 'min': -0.5, 'max': 7.5}
    assert old_top['hMixSamplerQA']['yAxis'] == {'nBins': 4, 'min': -0.5, 'max': 3.5}
    assert all(nuclear[k] == v for k, v in old_nuclear.items())
    he3_top = load('config/hist/hist_femtoLambda_he3_kf.yaml')['histograms']
    for key, spec in he3_top.items():
        if 'lambda_he3' not in key and key != 'hNucleus_PtVsYLab_he3':
            continue
        assert top[key] == spec
        mapped = key.replace('lambda_he3', 'lambda_he3_as_he4')
        if key == 'hNucleus_PtVsYLab_he3':
            mapped += '_as_he4'
        for field in spec:
            if field != 'title':
                assert top[mapped][field] == spec[field]
    he3_nuclear = load('config/hist/hist_femtoLambda_he3_nuclear.yaml')['histograms']
    pair_keys = [key for key in he3_nuclear if key.startswith(
        ('hKstar_', 'hKstarMass_', 'hKstarMt_', 'hQlab_', 'hDphiDeta_proton_'))]
    assert len(pair_keys) == 176
    for key in pair_keys:
        assert nuclear[key] == he3_nuclear[key]
        mapped = key.replace('3He', '3HeAs4He')
        for field, value in he3_nuclear[key].items():
            if field != 'title':
                assert nuclear[mapped][field] == value
    print('PASS helium control contract: 3 species, 8 channels, 690 histograms; immutable physics definitions retained; mixing QA covers all 8 channel indices')


if __name__ == '__main__':
    main()
