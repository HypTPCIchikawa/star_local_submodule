#!/usr/bin/env python3
"""Audit isolated log-poly study configs; never switch production or run events."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import shutil

import yaml


PROJECT = Path(__file__).resolve().parents[2]
SPECIES = ('d', 't', 'he3', 'he4')
BASES = ('main_auau3p85_anaFemtoLambda_4He_KFParticle',
         'main_auau3p85_anaFemtoLambda_4He_KFParticle_highpurity',
         'main_auau3p85_anaFemtoLambda_d_KFParticle_highpurity',
         'main_auau3p85_anaFemtoLambda_t_KFParticle_highpurity',
         'main_auau3p85_anaFemtoLambda_3He_KFParticle_highpurity')
OLD_PROFILE = 'cuts/nuclearid/nuclearid_femtoLambda_auau3p85_purity50_production_20261017.yaml'
NEW_PROFILE = 'cuts/nuclearid/nuclearid_femtoLambda_auau3p85_logpoly_study_20261019.yaml'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            digest.update(block)
    return digest.hexdigest()


def load(path):
    doc = yaml.safe_load(Path(path).read_text())
    require(isinstance(doc, dict), 'YAML root must be a mapping: ' + str(path))
    return doc


def audit_profile(old, new):
    require(set(old) == set(new), 'Nuclear profile keys changed')
    changed = {key for key in old if old[key] != new[key]}
    require(changed == {'dedxCalibration', 'speciesPOverQVeto'},
            'Only calibration mode and diagnostic veto removal may change')
    require(new['dedxCalibration'] == {'enabled': False},
            'Study must explicitly disable the table, without a calibrationFile')
    require(new['speciesPOverQVeto'] == {name: [] for name in SPECIES},
            'Study must diagnose all historical holes, without new vetoes')
    tof = new['highRigidityTof']
    require(new['maxPOverQ'] == 5.0 and tof['enabled'] is True
            and tof['thresholdPOverQ'] == 2.0
            and tof['extendHe4BeyondMaxPOverQ'] is False,
            'Require common R<5 and all-species TOF at R>=2')
    require(set(tof['species']) == set(SPECIES), 'All four species need explicit TOF')
    require(tof['species']['he4'] == tof['species']['d'], 'He4 must use the D window')
    require(new['requireBestSpecies'] is True, 'Keep original winner selection')
    return sorted(changed)


def parse_header(text):
    """Read coefficients from the actual header, not a duplicate calibration table."""
    tables = {}
    for name, count in (('kMean', 4), ('kSigma', 4), ('kPRange', 2)):
        match = re.search(r'\b' + name + r'\[kNSpecies\]\s*=\s*\{(.*?)\};', text, re.S)
        require(match is not None, 'Missing original polynomial table: ' + name)
        rows = re.findall(r'\{([^{}]+)\}', match.group(1))
        require(len(rows) == 4, 'Expected four species in ' + name)
        parsed = [tuple(float(cell.strip()) for cell in row.split(',')) for row in rows]
        require(all(len(row) == count and all(math.isfinite(x) for x in row) for row in parsed),
                'Invalid finite polynomial coefficients: ' + name)
        tables[name] = parsed
    require(all(row[3] == 0 for row in tables['kSigma']), 'Original widths must be quadratic')
    require('TMath::Log10(p)' in text, 'Original variable must remain log10(R)')
    expression = 'returnpar.c0+par.c1*lp+par.c2*lp*lp+par.c3*lp*lp*lp;'
    require(expression in re.sub(r'\s+', '', text), 'Original polynomial evaluation changed')
    return tables


def polynomial_direct(coefficients, x):
    return coefficients[0] + coefficients[1]*x + coefficients[2]*x*x + coefficients[3]*x*x*x


def polynomial_horner(coefficients, x):
    return ((coefficients[3]*x + coefficients[2])*x + coefficients[1])*x + coefficients[0]


def dense_check(tables, minima, cap, points, stream=None):
    require(points >= 3, 'Use at least three explicit dense-R points')
    if stream:
        stream.write('species\trigidity\tmean\tsigma\n')
    summaries = []
    checks = 0
    for index, name in enumerate(SPECIES):
        low = minima[name]
        valid_low, valid_high = tables['kPRange'][index]
        require(valid_low <= low < cap <= valid_high, 'Study outside legacy domain: ' + name)
        minimum_mean = minimum_sigma = math.inf
        minimum_sigma_r = None
        maximum_disagreement = 0.0
        for point in range(points):
            r = low + (cap-low)*point/(points-1)
            x = math.log10(r)
            mean = polynomial_direct(tables['kMean'][index], x)
            sigma = polynomial_direct(tables['kSigma'][index], x)
            for key, direct in (('kMean', mean), ('kSigma', sigma)):
                horner = polynomial_horner(tables[key][index], x)
                require(math.isfinite(direct) and direct > 0, 'Nonpositive/nonfinite legacy curve')
                require(abs(direct-horner) <= 2.e-13*max(1., abs(direct)),
                        'Independent polynomial evaluation disagrees')
                maximum_disagreement = max(maximum_disagreement, abs(direct-horner))
                checks += 2
            minimum_mean = min(minimum_mean, mean)
            if sigma < minimum_sigma:
                minimum_sigma, minimum_sigma_r = sigma, r
            if stream:
                stream.write('{}\t{:.17g}\t{:.17g}\t{:.17g}\n'.format(name, r, mean, sigma))
        summaries.append(dict(species=name, low=low, capExclusive=cap,
                              endpointIncludedForMathematicsOnly=True, points=points,
                              minimumMean=minimum_mean, minimumSigma=minimum_sigma,
                              minimumSigmaRigidity=minimum_sigma_r,
                              maximumIndependentEvaluationDifference=maximum_disagreement))
    return dict(checks=checks, species=summaries, runtimeProviderVerified=False,
                note='Source-coefficient arithmetic and positive smooth polynomial contract only; '
                     'not a C++ runtime, fit-quality, selection-purity or efficiency certification.')


def audit_configs(project):
    config = project/'config'
    sources = {config/OLD_PROFILE, config/NEW_PROFILE,
               project/'StMaker/common/NuclearIdDeDxVsMom.h',
               project/'StMaker/common/FemtoNuclearDedxCalibration.cxx',
               project/'StMaker/common/FemtoLambdaLegacy.cxx',
               Path(__file__).resolve(),
               project/'tools/nuclear_pid_cutoff/test_logpoly_20261019_config_audit.py'}
    old, new = load(config/OLD_PROFILE), load(config/NEW_PROFILE)
    changes = audit_profile(old, new)
    sources.add(config/old['dedxCalibration']['calibrationFile'])
    mains = []
    for base in BASES:
        paths = [config/'mainconf'/(base+suffix+'.yaml') for suffix in
                 ('', '_logpolyStudy20261019', '_preLogpoly20261019')]
        live, study, rollback = [load(path) for path in paths]
        require(live == rollback, 'Rollback is not YAML-equivalent to live: ' + base)
        require(live['nuclearid'] == OLD_PROFILE, 'Live nuclear selection was changed: ' + base)
        expected = dict(live, nuclearid=NEW_PROFILE)
        require(study == expected, 'Study main must differ only in nuclearid: ' + base)
        for path, doc in zip(paths, (live, study, rollback)):
            sources.add(path)
            for ref in doc.values():
                require(isinstance(ref, str), 'Main references must be strings')
                resolved = (config/ref).resolve(strict=True)
                require(resolved == project.resolve() or project.resolve() in resolved.parents,
                        'Reference escapes repository')
                sources.add(resolved)
            tag = load(config/doc['analysis'])['starTag']
            require(tuple(tag[key] for key in ('collision', 'year', 'productionTag', 'libraryTag'))
                    == ('auau3p85', 'year2021', 'P24iy', 'SL24y'), 'Wrong dataset/toolchain')
        mains.append(dict(live=str(paths[0].relative_to(project)),
                          study=str(paths[1].relative_to(project)),
                          rollback=str(paths[2].relative_to(project)),
                          liveUnchanged=True, rollbackYamlEquivalent=True,
                          studyChangedKeys=['nuclearid']))
    qa_on = config/'mainconf'/(BASES[0]+'_logpolyStudy20261019.yaml')
    qa_off = config/'mainconf'/(BASES[0]+'_logpolyStudyQaOff20261019.yaml')
    doc = load(qa_on)
    require('nuclearPidPairQaHist' in doc, 'Study canonical QA must remain on')
    expected_off = {key: value for key, value in doc.items() if key != 'nuclearPidPairQaHist'}
    require(load(qa_off) == expected_off, 'QA-off changes physics selection')
    sources.add(qa_off)
    return dict(profileChangedKeys=changes, mainconfs=mains,
                qaOn=str(qa_on.relative_to(project)), qaOff=str(qa_off.relative_to(project)),
                qaOnlyDifference=True, productionSwitched=False, purityCertified=False), sources, new


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--points', required=True, type=int)
    args = parser.parse_args()
    require(args.points >= 3, 'Use at least three explicit dense-R points')
    output = args.output_dir.resolve()
    require(not output.exists(), 'Output directory must be new')
    audit, sources, profile = audit_configs(PROJECT)
    protected = {str(path.relative_to(PROJECT)): sha(path) for path in sources}
    tables = parse_header((PROJECT/'StMaker/common/NuclearIdDeDxVsMom.h').read_text())
    output.mkdir(parents=True)
    for path in sources:
        dest = output/'source'/path.relative_to(PROJECT)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, dest)
        require(sha(dest) == protected[str(path.relative_to(PROJECT))], 'Snapshot differs')
    with (output/'dense_curve.tsv').open('x') as stream:
        dense = dense_check(tables, profile['speciesMinPOverQ'], profile['maxPOverQ'],
                            args.points, stream)
    with (output/'accepted_ranges.tsv').open('x') as stream:
        stream.write('species\tmin\tmax\n')
        for name in SPECIES:
            stream.write('{}\t{}\t{}\n'.format(name, profile['speciesMinPOverQ'][name],
                                              profile['maxPOverQ']))
    with (output/'tof_windows.tsv').open('x') as stream:
        stream.write('species\tstudiedMin\tstudiedMax\n')
        for name in SPECIES:
            window = profile['highRigidityTof']['species'][name]
            stream.write('{}\t{}\t{}\n'.format(name, window['mass2Min'], window['mass2Max']))
    report = dict(complete=True, audit=audit, denseR=dense, protected=protected,
                  sourceCoefficients=tables, productionSwitched=False,
                  physicsPurityCertified=False, noEventLoop=True)
    (output/'audit.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    require(all(sha(PROJECT/path) == digest for path, digest in protected.items()),
            'Protected source/config changed during audit')
    products = {str(path.relative_to(output)): sha(path)
                for path in output.rglob('*') if path.is_file()}
    (output/'completed.json').write_text(json.dumps(dict(complete=True, products=products,
                                                       productionSwitched=False), indent=2)+'\n')
    print('LOGPOLY_STUDY_CONFIG_AUDIT_PASS', dense['checks'], output)


if __name__ == '__main__':
    main()
