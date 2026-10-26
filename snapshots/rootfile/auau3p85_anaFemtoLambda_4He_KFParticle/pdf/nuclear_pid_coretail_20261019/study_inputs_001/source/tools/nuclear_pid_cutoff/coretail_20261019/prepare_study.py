#!/usr/bin/env python3
"""Prepare isolated core+tail study inputs; never fit, freeze, or switch production."""
import argparse
import copy
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import shutil

import yaml


PROJECT = Path(__file__).resolve().parents[3]
LOW_PROFILE = 'cuts/nuclearid/nuclearid_femtoLambda_auau3p85_logpoly_lowr_reference_20261019.yaml'
LOW_MAIN = 'config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_logpolyLowRReference20261019.yaml'
OUTPUT_BASE = 'rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_coretail_20261019'
TAIL_KEYS = ('tailFraction', 'tailOffsetSeedSigma', 'tailWidthRatio')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            result.update(block)
    return result.hexdigest()


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def owned(project, path, must_exist=True):
    path = Path(path)
    result = (path if path.is_absolute() else project/path).resolve(strict=must_exist)
    require(project.resolve() in result.parents, 'Require a project-owned path: ' + str(path))
    return result


def own_source_paths(script_path):
    """Python preserves a relative __file__ when invoked by a relative CLI path."""
    script_path = Path(script_path).resolve(strict=True)
    return script_path, script_path.with_name('test_prepare_study.py')


def validate_reference(cuts):
    forbidden = {'speciesMinPOverQ', 'speciesPOverQVeto', 'highRigidityTof',
                 'he4HighRigidityTof'} & set(cuts)
    require(not forbidden, 'Low-R reference must omit per-species/TOF selectors')
    require(cuts.get('dedxCalibration') == {'enabled': False},
            'Require explicit legacy-polynomial calibration off')
    require(cuts.get('maxPOverQ') == 2.0 and cuts.get('m2_selection') is False,
            'Require TPC-only exclusive R<2 reference')
    require(cuts.get('requireBestSpecies') is True, 'Keep original best-species gate')


def derive_settings(base, main_reference, output_base):
    require(base.get('schemaVersion') == 1, 'Unsupported base settings schema')
    require(base.get('fitWindowMode') == 'target_local', 'Use target-local control settings')
    require(base.get('allFit', {}).get('targetShapeMode') == 'free_target_mean_sigma',
            'Require the existing free-target control comparison')
    require(0.4 <= base['rigidityMin'] < base['rigidityMax'] <= 2.0,
            'Do not extrapolate below the shared legacy domain or into joint high-R TOF')
    require(set(base.get('responseModels', [])) == {'gaussian', 'lognormal'},
            'Preserve both existing Gaussian/log-normal controls')
    require(set(base.get('backgroundModels', [])) == {'linear', 'exponential'},
            'Preserve both existing background controls')
    result = copy.deepcopy(base)
    result['runtimeMainconf'] = main_reference
    result['baselineMainconf'] = main_reference
    result['outputBase'] = output_base
    require(all(result[key] == base[key] for key in base
                if key not in ('runtimeMainconf', 'baselineMainconf', 'outputBase')),
            'Unexpected fit/quality/input changes')
    return result


def validate_tail_profiles(doc):
    """Validate explicit fixed alternatives; provide no implicit grid or fitted tail."""
    require(isinstance(doc, dict) and set(doc) == {'schemaVersion', 'profiles'}
            and doc['schemaVersion'] == 1, 'Explicit tail profile schema must be version 1')
    profiles = doc['profiles']
    require(isinstance(profiles, list) and bool(profiles), 'Provide nonempty explicit profiles')
    labels = set()
    for profile in profiles:
        require(isinstance(profile, dict) and set(profile) == {'label', *TAIL_KEYS},
                'Fixed profile requires exactly label/fraction/offset/ratio')
        label = profile['label']
        require(isinstance(label, str) and re.fullmatch(r'[A-Za-z0-9_]+', label)
                and label not in labels, 'Profile labels must be safe and unique')
        labels.add(label)
        require(all(isinstance(profile[key], (int, float)) and not isinstance(profile[key], bool)
                    and math.isfinite(profile[key]) for key in TAIL_KEYS),
                'Tail controls must be finite numbers, not bools')
        fraction, ratio = profile['tailFraction'], profile['tailWidthRatio']
        require(0 <= fraction < .5, 'Use fixed core-dominant tail fraction in [0,0.5)')
        require(ratio >= 1 if fraction == 0 else ratio > 1,
                'Nonzero tail must be wider than the Gaussian core')
    return copy.deepcopy(doc)


def inspect_reference(project, mainconf):
    doc = yaml.safe_load(mainconf.read_text())
    require(doc['nuclearid'] == LOW_PROFILE, 'Use the explicit isolated low-R reference profile')
    cuts_path = owned(project, project/'config'/doc['nuclearid'])
    cuts = yaml.safe_load(cuts_path.read_text())
    validate_reference(cuts)
    tag = yaml.safe_load((project/'config'/doc['analysis']).read_text())['starTag']
    require(tuple(tag[key] for key in ('collision', 'year', 'productionTag', 'libraryTag'))
            == ('auau3p85', 'year2021', 'P24iy', 'SL24y'), 'Require 2021 P24iy / SL24y')
    runtime_path = project/'tools/nuclear_pid_tof/runtime_selection.py'
    freeze_path = project/'tools/nuclear_pid_cutoff/freeze_baseline.py'
    runtime = module('coretail_reference_runtime', runtime_path)
    selection = runtime.load_current_selection(project, mainconf)
    module('coretail_reference_freezer', freeze_path).validate_reference_cuts(selection['cuts'])
    require(selection['calibration'] is None, 'Reader must use only original polynomials')
    paths = {Path(path) for path in selection['protectedFiles']}
    paths.update((runtime_path, freeze_path, mainconf, cuts_path,
                  project/'StMaker/common/FemtoNuclearDedxCalibration.cxx',
                  project/'StMaker/common/FemtoLambdaLegacy.cxx'))
    paths.update(owned(project, project/'config'/ref) for ref in doc.values())
    return selection, paths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('base-settings', 'mainconf', 'output-dir'):
        parser.add_argument('--'+name, required=True, type=Path)
    parser.add_argument('--tail-profiles', type=Path,
                        help='Optional explicit fixed profiles; never invent a default grid')
    args = parser.parse_args()
    base_path, mainconf = (owned(PROJECT, path) for path in (args.base_settings, args.mainconf))
    output = owned(PROJECT, args.output_dir, must_exist=False)
    output_base = owned(PROJECT, OUTPUT_BASE, must_exist=False)
    require(output_base in output.parents, 'Write only below the new coretail output base')
    require(not output.exists(), 'Never overwrite a prepared study')
    selection, paths = inspect_reference(PROJECT, mainconf)
    base = yaml.safe_load(base_path.read_text())
    cfg = derive_settings(base, str(mainconf.relative_to(PROJECT)), OUTPUT_BASE)
    require(cfg['allFit']['maxNSigma'] == selection['cuts']['maxNSigmaNuclear'],
            'Control fit and actual fixed PID gate must agree')
    input_study = owned(PROJECT, cfg['inputStudy'])
    context_path, completed_path = input_study/'context.json', input_study/'completed.json'
    context, completed = json.loads(context_path.read_text()), json.loads(completed_path.read_text())
    require(completed.get('executionComplete') and completed.get('protectedPostRunVerified'),
            'Reuse only the completed native input cache')
    require('2021P24iy_partial_merge' in Path(context['inputFingerprint']['path']).name,
            'Original cache must come from the 2021 P24iy partial merge')
    cache = input_study/'cache.root'
    require(digest(cache) == completed['products']['cache.root'], 'Immutable cache SHA changed')
    paths.update((base_path, context_path, completed_path, input_study/'native_audit.json'))
    paths.update(own_source_paths(__file__))
    tail = None
    if args.tail_profiles:
        tail_path = owned(PROJECT, args.tail_profiles)
        tail = validate_tail_profiles(yaml.safe_load(tail_path.read_text()))
        paths.add(tail_path)
    protected = {str(path.relative_to(PROJECT)): digest(path) for path in paths}
    output.mkdir(parents=True)
    for path in paths:
        target = output/'source'/path.relative_to(PROJECT)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        require(digest(target) == protected[str(path.relative_to(PROJECT))], 'Source snapshot changed')
    (output/'free_local_settings.yaml').write_text(yaml.safe_dump(cfg, sort_keys=False))
    if tail is not None:
        (output/'fixed_tail_profiles.yaml').write_text(yaml.safe_dump(tail, sort_keys=False))
    report = dict(complete=True, productionSwitched=False, fitsExecuted=False,
                  baselineFreezeExecuted=False, purityCertified=False,
                  fixedTailProfilesSupplied=tail is not None,
                  tailGridPolicy='Only user-supplied fixed profiles; no automatic tail grid or free mixture',
                  tailOffsetConvention='Fixed dE/dx offset seeded by tailOffsetSeedSigma times '
                                       'the seed core sigma, not the fitted core sigma',
                  coretailRuntimeConfigured=False,
                  note='Prepared Gaussian/log-normal controls plus optional profile inputs. '
                       'The core+tail engine connection must be supplied separately before fits.',
                  mainconf=str(mainconf), inputCache=str(cache),
                  inputCacheSha256=completed['products']['cache.root'], protected=protected,
                  generatedSettingsChangedKeys=['runtimeMainconf', 'baselineMainconf', 'outputBase'])
    (output/'preparation.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    require(all(digest(PROJECT/path) == expected for path, expected in protected.items()),
            'Protected reference changed during preparation')
    require(digest(cache) == report['inputCacheSha256'], 'Input cache changed during preparation')
    products = {str(path.relative_to(output)): digest(path)
                for path in output.rglob('*') if path.is_file()}
    (output/'completed.json').write_text(json.dumps(dict(complete=True, products=products,
                                                       productionSwitched=False), indent=2)+'\n')
    print('CORETAIL_STUDY_INPUTS_PREPARED', output)


if __name__ == '__main__':
    main()
