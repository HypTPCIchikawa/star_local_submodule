#!/usr/bin/env python3
"""Run a frozen, opt-in ROOT5 joint helium coefficient fit; never publish PID."""
import argparse
from datetime import datetime, timezone, timedelta
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys

import yaml

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent))
from calibration_candidates import polynomial_fit, polynomial_value


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fingerprint(path):
    st = path.stat()
    return dict(path=str(path), size=st.st_size, mtime_ns=st.st_mtime_ns,
                device=st.st_dev, inode=st.st_ino)


def save(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def legacy_coefficients(header):
    text = header.read_text()
    result = []
    arrays = []
    for name in ('kMean', 'kSigma'):
        block = re.search(r'\b'+name+r'\[kNSpecies\]\s*=\s*\{(.*?)\};', text, re.S)
        if not block:
            raise ValueError('Cannot read the frozen legacy coefficient table '+name)
        rows = [[float(x.strip()) for x in row.split(',')]
                for row in re.findall(r'\{([^{}]+)\}', block.group(1))]
        if len(rows) != 4 or any(len(row) != 4 for row in rows):
            raise ValueError('Legacy table structure changed.')
        arrays.append(rows)
    for species in (2, 3):
        if arrays[1][species][3] != 0:
            raise ValueError('Existing sigma formula is no longer quadratic.')
        result += arrays[0][species]+arrays[1][species][:3]
    return result


def extrema(coefficients, lo, hi):
    c = list(coefficients)+[0.]*(4-len(coefficients))
    points = [lo, hi]
    a, b, d = 3*c[3], 2*c[2], c[1]
    if abs(a) < 1.e-14:
        if abs(b) > 1.e-14:
            points.append(-d/b)
    elif b*b-4*a*d >= 0:
        discriminant = math.sqrt(b*b-4*a*d)
        points += [(-b+discriminant)/(2*a), (-b-discriminant)/(2*a)]
    values = [polynomial_value(c, x) for x in points if lo <= x <= hi]
    return min(values), max(values)


def feasible(values, legacy, cfg):
    if len(values) != 14 or not all(math.isfinite(v) for v in values):
        return False
    lo, hi = math.log10(cfg['rigidityMin']), math.log10(cfg['rigidityMax'])
    for offset in (0, 7):
        mean, width = values[offset:offset+4], values[offset+4:offset+7]
        if extrema(mean, lo, hi)[0] <= 0 or extrema(width, lo, hi)[0] < cfg['minimumSigma']:
            return False
        if cfg['requireMonotonicMean'] and extrema([mean[1], 2*mean[2], 3*mean[3]], lo, hi)[1] > 0:
            return False
        for bound, endpoint in (('meanShiftMinSigma', 0), ('meanShiftMaxSigma', 1)):
            difference = [mean[j]-legacy[offset+j]-cfg[bound]*(legacy[offset+4+j] if j<3 else 0.)
                          for j in range(4)]
            extremum = extrema(difference, lo, hi)[endpoint]
            if (endpoint == 0 and extremum < 0) or (endpoint == 1 and extremum > 0):
                return False
        upper = [width[j]-cfg['maximumSigmaScale']*legacy[offset+4+j] for j in range(3)]
        if extrema(upper, lo, hi)[1] > 0:
            return False
        scales = [1., 3., 5., 10., .1, .3, .5]
        if any(abs(values[offset+j]-legacy[offset+j]) > cfg['coefficientBoundSteps']*scales[j]
               for j in range(7)):
            return False
    return extrema([values[7+j]-values[j] for j in range(4)], lo, hi)[0] >= cfg['minimumMeanGap']


def blend_candidate(candidate, legacy, cfg, generation):
    if feasible(candidate, legacy, cfg):
        return candidate, 1.
    low, high = 0., 1.
    for _ in range(generation['blendBisectionIterations']):
        fraction = (low+high)/2
        trial = [a+fraction*(b-a) for a,b in zip(legacy, candidate)]
        if feasible(trial, legacy, cfg):
            low = fraction
        else:
            high = fraction
    fraction = low*generation['blendSafetyFactor']
    trial = [a+fraction*(b-a) for a,b in zip(legacy, candidate)]
    if not feasible(trial, legacy, cfg):
        raise ValueError('Cannot create a feasible common-coefficient seed.')
    return trial, fraction


def make_seeds(legacy, cfg, generation, knot_path):
    seeds, audit = [], []
    for mode in generation['modes']:
        candidate = list(legacy)
        detail = dict(name=mode, interpretation='Initialization only; not truth, prior or fit constraint.')
        if mode == 'legacy_low_shift':
            # A polynomial in log(R), not independent per-bin mean/width parameters.
            pivot = math.log10(cfg['rigidityMax'])
            span = math.log10(cfg['rigidityMin'])-pivot
            delta = generation['lowShiftMeanSigma']*polynomial_value(legacy[11:14], math.log10(cfg['rigidityMin']))
            candidate[7] += -pivot*delta/span
            candidate[8] += delta/span
            candidate[11:14] = [x*generation['lowShiftWidthScale'] for x in legacy[11:14]]
        elif mode == 'old_knots_poly':
            doc = yaml.safe_load(knot_path.read_text())
            for index, name in enumerate(('he3', 'he4')):
                knots = doc['species'][name]['knots']
                xs = [math.log10(k['rigidity']) for k in knots]
                candidate[index*7:index*7+4] = polynomial_fit(xs, [k['mean'] for k in knots], 3)
                candidate[index*7+4:index*7+7] = polynomial_fit(xs, [k['sigma'] for k in knots], 2)
            detail['source'] = str(knot_path)
            detail['sourceSha256'] = digest(knot_path)
            detail['note'] = 'Regression of the old independent-fit knots ONLY to initialize the raw-data fit.'
        elif mode != 'legacy':
            raise ValueError('Unknown seed mode '+mode)
        values, fraction = blend_candidate(candidate, legacy, cfg, generation)
        detail.update(unblendedCoefficients=candidate, coefficients=values, feasibleBlendFraction=fraction)
        seeds.append((mode, values)); audit.append(detail)
    if not seeds or len({name for name,_ in seeds}) != len(seeds):
        raise ValueError('Seeds must be nonempty and uniquely named.')
    return seeds, audit


def validate(cfg):
    sample, engine, generation = cfg['sample'], cfg['engine'], cfg['seedGeneration']
    if cfg['schemaVersion'] != 1 or (sample['collision'], sample['year'], sample['production']) != ('auau3p85',2021,'P24iy'):
        raise ValueError('Only auau3p85 / 2021 P24iy is authorized for this study.')
    if engine['sampleMode'] not in ('all','tag_not_tag') or engine['backgroundModel'] not in ('linear','exponential'):
        raise ValueError('Unknown data partition or background model.')
    if not 0.4 <= engine['rigidityMin'] < engine['rigidityMax'] <= 2.5:
        raise ValueError('Unsupported calibration domain.')
    if engine['rigidityStep'] not in (0.01,0.02) or engine['rQuadratureOrder'] not in (1,3,5):
        raise ValueError('R bins must align with native input; quadrature must be 1,3 or 5.')
    if engine['dedxRebinFactor'] not in (5,10):
        raise ValueError('Use alias-safe dE/dx bins: 0.5 or 1 keV/cm.')
    n = round((engine['rigidityMax']-engine['rigidityMin'])/engine['rigidityStep'])
    indices = engine['sliceIndices']
    if not isinstance(indices,list) or len(set(indices)) != len(indices) or any(
            isinstance(i,bool) or not isinstance(i,int) or not 0<=i<n for i in indices):
        raise ValueError('Invalid/duplicate native-aligned slice indices.')
    if not 0 < generation['blendSafetyFactor'] < 1 or generation['blendBisectionIterations'] < 1:
        raise ValueError('Invalid seed blending policy.')
    if cfg['executionTimeoutSeconds'] <= 0:
        raise ValueError('Invalid ROOT execution timeout.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--settings', type=Path, default=HERE/'global_logpoly_settings.yaml')
    parser.add_argument('--mainconf', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--prepare-only', action='store_true')
    parser.add_argument('--seeds', nargs='+', choices=['legacy','legacy_low_shift','old_knots_poly'])
    args = parser.parse_args()
    settings_path = args.settings.resolve(strict=True)
    cfg = yaml.safe_load(settings_path.read_text())
    if args.seeds:
        cfg['seedGeneration']['modes'] = args.seeds
    validate(cfg)
    mainconf = args.mainconf.resolve(strict=True)
    sample, engine = cfg['sample'], cfg['engine']
    if mainconf != (PROJECT/sample['mainconf']).resolve(strict=True):
        raise ValueError('Explicit mainconf differs from the configured runtime.')
    runtime = yaml.safe_load(mainconf.read_text())
    analysis_path = (PROJECT/'config'/runtime['analysis']).resolve(strict=True)
    analysis = yaml.safe_load(analysis_path.read_text())
    tags = analysis.get('starTag', {})
    if (str(tags.get('libraryTag')), str(tags.get('year')), str(tags.get('productionTag')),
            str(tags.get('collision')), str(analysis.get('analysis',{}).get('mode'))) != ('SL24y','year2021','P24iy','auau3p85','fxtmult'):
        raise ValueError('Runtime must be SL24y / 2021 P24iy / auau3p85 FXT.')
    source = (PROJECT/sample['input']).resolve(strict=True)
    source_summary = (PROJECT/sample['sourceSummary']).resolve(strict=True)
    state = json.loads(source_summary.read_text())
    initial = fingerprint(source)
    if state['phase'] != 'complete' or state['provenance']['input'] != initial:
        raise ValueError('Original ROOT differs from the verified source summary.')
    output = args.output_dir.resolve()
    # Even an empty existing directory is rejected: no uncertain RECREATE target.
    if output.exists():
        raise ValueError('Choose a nonexistent output directory; existing products are never replaced.')
    allowed = (PROJECT/sample['outputBase']).resolve()
    if allowed not in output.parents or output == allowed:
        raise ValueError('Study outputs must be a separate child of '+str(allowed))
    header = PROJECT/'StMaker/common/NuclearIdDeDxVsMom.h'
    legacy = legacy_coefficients(header)
    if not feasible(legacy, legacy, engine):
        raise ValueError('Legacy seed is infeasible under this explicit study configuration.')
    knots = (PROJECT/sample['legacyKnots']).resolve(strict=True)
    seeds, seed_audit = make_seeds(legacy, engine, cfg['seedGeneration'], knots)
    output.mkdir(parents=True)
    (output/'logs').mkdir()
    copies = [HERE/'GlobalHeliumFit.C', Path(__file__).resolve(),
              HERE.parent/'calibration_candidates.py', header,
              PROJECT/'script/singularity_root_batch.sh', mainconf, analysis_path, settings_path, knots]
    source_hashes = {}
    for src in copies:
        dst = output/'source'/src.relative_to(PROJECT)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(src),str(dst))
        source_hashes[str(src)] = digest(src)
        if digest(dst) != source_hashes[str(src)]:
            raise ValueError('Source copy verification failed.')
    (output/'settings.yaml').write_text(yaml.safe_dump(cfg, sort_keys=False,allow_unicode=True))
    (output/'settings.tsv').write_text(''.join('%s\t%s\n'%(key, ','.join(map(str,value)) if isinstance(value,list)
        else int(value) if isinstance(value,bool) else value) for key,value in engine.items()))
    (output/'seeds.tsv').write_text(''.join(name+'\t'+'\t'.join(format(v,'.17g') for v in values)+'\n'
                                            for name,values in seeds))
    save(output/'seed_audit.json',seed_audit)
    context = dict(schemaVersion=1,input=initial,settings=cfg,mainconf=str(mainconf),
        mainconfSha256=digest(mainconf),sourceSummary=str(source_summary),sourceSummarySha256=digest(source_summary),
        sources=source_hashes,createdJst=datetime.now(timezone(timedelta(hours=9))).isoformat(),
        productionEnabled=False,note='Charge-inclusive, centrality-integrated, empirical component model; no CF subtraction.')
    save(output/'context.json',context)
    macro = output/'source/tools/he4_purity/global_logpoly/GlobalHeliumFit.C'
    call = str(macro)+'++('+','.join(json.dumps(str(p)) for p in
        (source,output,output/'settings.tsv',output/'seeds.tsv'))+')'
    command = ['bash',str(PROJECT/'script/singularity_root_batch.sh'),str(mainconf),
               'root4star','-l','-b','-q',call]
    save(output/'command.json',command)
    print('PREPARED' if args.prepare_only else 'RUN',str(output),flush=True)
    if args.prepare_only:
        return
    log = output/'logs/fit.log'
    with log.open('w') as stream:
        stream.write('COMMAND: '+json.dumps(command)+'\n');stream.flush()
        result = subprocess.run(command,cwd=str(PROJECT),stdout=stream,stderr=subprocess.STDOUT,
                                timeout=cfg['executionTimeoutSeconds'])
    contents = log.read_text(errors='replace')
    if result.returncode or 'GLOBAL_HELIUM_COMPLETE' not in contents or any(
            word in contents for word in ('Error in <ACLiC>','segmentation violation','fatal error:','Error: Function')):
        raise RuntimeError('Joint fit did not finish cleanly; see '+str(log))
    data = json.loads((output/'global_fit.json').read_text(),parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))
    if data['input'] != str(source) or data['globalParameterCount'] != 14 or len(data['runs']) != len(seeds):
        raise ValueError('Unexpected source or number of fitted shared coefficients/seeds.')
    for run in data['runs']:
        if not feasible(run['coefficients'],legacy,engine):
            raise ValueError('Joint fit exported a nonphysical/unsupported coefficient function.')
    if fingerprint(source) != initial or any(digest(Path(p))!=v for p,v in source_hashes.items()):
        raise ValueError('Input/source changed during execution; this study is not publishable.')
    products = {str(p.relative_to(output)):digest(p) for p in
                (output/'global_fit.json',output/'input_slices.root',output/'settings.yaml',output/'settings.tsv',output/'seeds.tsv')}
    save(output/'completed.json',dict(executionComplete=True,productionEnabled=False,
        physicsValidationComplete=False,products=products,
        completedJst=datetime.now(timezone(timedelta(hours=9))).isoformat(),
        note='Execution completion is NOT fit convergence, identifiable isotope assignment or validated PID.'))
    print('GLOBAL_RUNNER_COMPLETE',str(output),flush=True)


if __name__ == '__main__':
    main()
