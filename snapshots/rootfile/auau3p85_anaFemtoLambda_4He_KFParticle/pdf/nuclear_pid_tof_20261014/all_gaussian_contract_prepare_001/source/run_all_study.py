#!/usr/bin/env python3
"""Immutable, study-only fixed-curve nuclear fits and conditional TPC purity."""
import argparse
from datetime import datetime, timedelta, timezone
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

import yaml

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
SAMPLES = ('ALL', 'REMAINING', 'TAG_D', 'TAG_T', 'TAG_HE3')
SPECIES = ('d', 't', 'he3', 'he4')
TOP_KEYS = set('''schemaVersion runtimeMainconf inputStudy outputBase legacyHeader
executionTimeoutSeconds backgroundModels partitions tagComponents primaryComponents seed numeric
support curves allFit sourceSummary protectedRelativePaths'''.split())
READBACK_SOURCE = r'''#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TSystem.h"
#include "TString.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
void ValidateAllProduct(const char* cachePath,const char* productPath,const char* sample,
                        int firstR,int rows,double axisTolerance,double sumTolerance){
  try{
    TFile cache(cachePath,"READ"),product(productPath,"READ");
    if(cache.IsZombie()||product.IsZombie()||product.TestBit(TFile::kRecovered))
      throw std::runtime_error("unreadable/recovered ALL product");
    TH2* input=dynamic_cast<TH2*>(cache.Get(sample));if(!input)throw std::runtime_error("missing exact sample");
    const char* parts[]={"smooth","d","t","he3","he4"};
    for(int index=0;index<rows;++index){
      TH1* raw=dynamic_cast<TH1*>(product.Get(Form("raw_native_%03d",index)));
      TH1* total=dynamic_cast<TH1*>(product.Get(Form("prediction_total_%03d",index)));
      TH1* mask=dynamic_cast<TH1*>(product.Get(Form("prediction_fit_mask_%03d",index)));
      if(!raw||!total||!mask||raw->GetNbinsX()!=input->GetNbinsY())throw std::runtime_error("incomplete native/model product");
      TH1* component[5];for(int j=0;j<5;++j){component[j]=dynamic_cast<TH1*>(product.Get(Form("prediction_%s_%03d",parts[j],index)));if(!component[j])throw std::runtime_error("missing model component");}
      for(int b=0;b<=raw->GetNbinsX()+1;++b){
        double original=input->GetBinContent(firstR+index,b),saved=raw->GetBinContent(b);
        double error=input->GetBinError(firstR+index,b),savedError=raw->GetBinError(b);
        if(original!=saved||std::fabs(error-savedError)>sumTolerance*std::max(1.,std::fabs(error)))
          throw std::runtime_error("native raw count/error/flow changed");
        if(b>=1&&b<=raw->GetNbinsX()&&std::fabs(raw->GetBinLowEdge(b)-input->GetYaxis()->GetBinLowEdge(b))>axisTolerance)
          throw std::runtime_error("native dEdx axis changed");
        double sum=0;for(int j=0;j<5;++j)sum+=component[j]->GetBinContent(b);double expected=total->GetBinContent(b);
        if(!TMath::Finite(sum)||!TMath::Finite(expected)||std::fabs(sum-expected)>sumTolerance*std::max(1.,std::fabs(expected)))
          throw std::runtime_error("model component sum changed");
        if(mask->GetBinContent(b)!=0&&mask->GetBinContent(b)!=1)throw std::runtime_error("invalid fit mask");
      }
    }
    product.Close();cache.Close();std::cout<<"NUCLEAR_ALL_READBACK_PASS rows="<<rows<<std::endl;
  }catch(const std::exception& e){std::cerr<<e.what()<<std::endl;gSystem->Exit(2);}
}
'''


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def fingerprint(path):
    path = Path(path)
    stat = path.stat()
    return dict(path=str(path), size=stat.st_size, mtime_ns=stat.st_mtime_ns,
                device=stat.st_dev, inode=stat.st_ino)


def save(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def read_json(path):
    def invalid(value):
        raise ValueError('Nonfinite JSON constant: '+value)
    return json.loads(Path(path).read_text(), parse_constant=invalid)


def project_path(value):
    path = Path(value)
    path = (path if path.is_absolute() else PROJECT/path).resolve(strict=True)
    if PROJECT not in path.parents:
        raise ValueError('Expected a project-owned input/source: '+str(path))
    return path


def finite(value, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError('Expected an explicit finite number.')
    if positive and value <= 0:
        raise ValueError('Expected a positive mathematical setting.')
    return value


def verify_manifest(folder, manifest, field, prefix=None):
    signatures = {}
    for name, expected in manifest[field].items():
        path = folder/name if prefix is None else folder/prefix/name
        if digest(path) != expected:
            raise ValueError('Completed source signature differs: '+str(path))
        signatures[str(path.resolve())] = expected
    return signatures


def load_snapshot_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def validate_candidate(candidate):
    if candidate.get('schemaVersion') != 1 or candidate.get('productionEnabled') is not False:
        raise ValueError('Only schemaVersion=1, productionEnabled=false curves are accepted.')
    if set(candidate.get('species', {})) != set(SPECIES):
        raise ValueError('Candidate must declare d/t/he3/he4, including unsupported empty segments.')
    for species in SPECIES:
        previous = 0.
        for segment in candidate['species'][species].get('segments', []):
            lo, hi = segment['domain']
            finite(lo, True)
            finite(hi, True)
            if hi <= lo or lo < previous or segment['family'] not in ('logpoly', 'pchip_logR'):
                raise ValueError('Reversed, overlapping or unknown curve segment.')
            previous = hi


def raw_parameters(runtime, tables, candidate, curves, species_index, rigidity):
    value = curves.evaluate_curve(candidate, SPECIES[species_index], rigidity)
    for field in ('mean', 'sigma', 'meanError', 'sigmaError'):
        finite(value[field])
    if value['supported']:
        if value['mean'] <= 0 or value['sigma'] <= 0 or value['meanError'] < 0 or value['sigmaError'] < 0:
            raise ValueError('Supported curve has invalid shape/errors.')
        return value, 'tag_curve_supported'
    mean, sigma = runtime.legacy_parameters(tables, species_index, rigidity,
                                           ignore_domain=True)
    if not math.isfinite(mean) or not math.isfinite(sigma) or mean <= 0 or sigma <= 0:
        raise ValueError('Unvalidated nuisance initialization is not mathematically positive.')
    return dict(mean=mean, sigma=sigma, meanError=0., sigmaError=0., supported=False), \
        'legacy_polynomial_unvalidated_nuisance_not_PID_fallback'


def write_tsv(path, header, rows):
    with path.open('w') as stream:
        if header:
            stream.write('\t'.join(header)+'\n')
        for row in rows:
            stream.write('\t'.join(format(v, '.17g') if isinstance(v, float) else str(v)
                                   for v in row)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ('settings', 'calibration', 'mainconf', 'output-dir'):
        parser.add_argument('--'+option, type=Path, required=True)
    parser.add_argument('--window-calibration', type=Path, action='append', default=[],
                        help='Repeat the same full family comparison set for identical fit windows.')
    parser.add_argument('--background', choices=('linear', 'exponential'), required=True)
    parser.add_argument('--sample', choices=SAMPLES, required=True)
    parser.add_argument('--rigidity-min', type=float)
    parser.add_argument('--rigidity-max', type=float)
    parser.add_argument('--prepare-only', action='store_true')
    args = parser.parse_args()
    settings_path = project_path(args.settings)
    cfg = yaml.safe_load(settings_path.read_text())
    needed = {'schemaVersion', 'runtimeMainconf', 'inputStudy', 'outputBase', 'legacyHeader',
              'executionTimeoutSeconds', 'backgroundModels', 'partitions', 'seed', 'numeric', 'allFit'}
    if cfg.get('schemaVersion') != 1 or not needed <= set(cfg) or set(cfg)-TOP_KEYS:
        raise ValueError('Invalid explicit ALL study configuration schema.')
    mainconf = project_path(args.mainconf)
    if mainconf != project_path(cfg['runtimeMainconf']):
        raise ValueError('Explicit mainconf differs from configured runtime.')
    if args.background not in cfg['backgroundModels'] or not cfg['allFit']['requireBestSpecies']:
        raise ValueError('Background or best-species policy is not declared in settings.')
    timeout = cfg['executionTimeoutSeconds']
    if isinstance(timeout, bool) or not isinstance(timeout, int) or timeout <= 0:
        raise ValueError('Execution timeout must be a positive declared integer.')
    inputs = project_path(cfg['inputStudy'])
    completion_path, context_path = inputs/'completed.json', inputs/'context.json'
    completed, input_context = read_json(completion_path), read_json(context_path)
    if not completed.get('executionComplete') or completed.get('productionEnabled') is not False \
            or not completed.get('protectedPostRunVerified'):
        raise ValueError('INPUT source is not complete/protected study-only data.')
    protected = verify_manifest(inputs, completed, 'products')
    protected.update(verify_manifest(inputs, completed, 'sourceSha256', 'source'))
    protected[str(completion_path)] = digest(completion_path)
    for path, expected in input_context['protectedPreRun'].items():
        if digest(Path(path)) != expected:
            raise ValueError('Previously protected original changed: '+path)
        protected[path] = expected
    source = Path(input_context['inputFingerprint']['path']).resolve(strict=True)
    if 'auau3p85' not in source.name or '2021P24iy_partial_merge' not in source.name \
            or fingerprint(source) != input_context['inputFingerprint']:
        raise ValueError('Only the original unchanged 3p85 / 2021 P24iy partial merge is allowed.')
    input_cfg = yaml.safe_load((inputs/'source/input_settings.yaml').read_text())
    source_summary = project_path(input_cfg['sourceSummary'])
    summary = read_json(source_summary)
    if summary['phase'] != 'complete' or summary['provenance']['input'] != fingerprint(source):
        raise ValueError('Original ROOT no longer matches its verified summary.')
    if input_context['runtimeMainconf'] != str(mainconf):
        raise ValueError('Runtime differs from immutable input preparation.')
    main = yaml.safe_load(mainconf.read_text())
    analysis_path = project_path(Path('config')/main['analysis'])
    analysis = yaml.safe_load(analysis_path.read_text())
    tag = analysis['starTag']
    if (tag['collision'], tag['year'], tag['productionTag'], tag['libraryTag']) != \
            ('auau3p85', 'year2021', 'P24iy', 'SL24y'):
        raise ValueError('Wrong collision/year/production/runtime.')
    audit = read_json(inputs/'native_audit.json')
    if not all(audit.get(k) is True for k in ('axesExactlyEqual', 'tagSumSubsetOfAll',
                                              'cacheCountErrorFlowPreserved', 'writtenCacheReadbackExact')):
        raise ValueError('Native cache count/axis/subset audit failed.')
    cache = (inputs/'cache.root').resolve(strict=True)
    axis = audit['nativeAxis']
    partition = cfg['partitions']['native']
    minimum = partition['min'] if args.rigidity_min is None else args.rigidity_min
    maximum = partition['max'] if args.rigidity_max is None else args.rigidity_max
    finite(minimum, True)
    finite(maximum, True)
    tolerance = finite(cfg['numeric']['axisTolerance'], True)
    step = finite(axis['rBinWidth'], True)
    if abs(partition['step']-step) > tolerance or maximum <= minimum:
        raise ValueError('Only the declared non-overlapping native-R partition is allowed.')
    rmin, rmax = axis['rRange']
    first = round((minimum-rmin)/step)+1
    last = round((maximum-rmin)/step)
    if first < 1 or last < first or last > axis['rBins'] or \
            abs(rmin+(first-1)*step-minimum) > tolerance or abs(rmin+last*step-maximum) > tolerance:
        raise ValueError('Requested R edges are not exact native bins.')
    calibration_path = project_path(args.calibration)
    window_paths = sorted({calibration_path}|{project_path(p) for p in args.window_calibration})
    candidate_objects = {path:read_json(path) for path in window_paths}
    for candidate in candidate_objects.values():
        validate_candidate(candidate)
    output = args.output_dir.resolve()
    allowed = project_path(cfg['outputBase'])
    if output.exists() or output.parent != allowed or output == inputs:
        raise ValueError('Choose a nonexistent sibling study under the configured new output base.')
    output.mkdir()
    snapshot = output/'source'
    snapshot.mkdir()
    source_files = [Path(__file__).resolve(), HERE/'FitAllNuclearPurity.C', HERE/'NuclearFitCommon.h',
                    HERE/'runtime_selection.py', HERE/'nuclear_curves.py', settings_path,
                    mainconf, analysis_path, PROJECT/'script/singularity_root_batch.sh']
    if len({p.name for p in source_files}) != len(source_files):
        raise ValueError('Source snapshot basename collision.')
    snapshots = {}
    for path in source_files:
        expected = digest(path)
        protected[str(path)] = expected
        dest = snapshot/path.name
        shutil.copy2(path, dest)
        if digest(dest) != expected:
            raise ValueError('Source changed while copying snapshot.')
        snapshots[str(dest)] = expected
    runtime = load_snapshot_module('runtime_selection', snapshot/'runtime_selection.py')
    curves = load_snapshot_module('nuclear_curves', snapshot/'nuclear_curves.py')
    current = runtime.load_current_selection(PROJECT, mainconf)
    if project_path(cfg['legacyHeader']) != Path(current['tables']['source']):
        raise ValueError('Actual header differs from configured legacy source.')
    if cfg['allFit']['maxNSigma'] != current['cuts']['maxNSigmaNuclear']:
        raise ValueError('The study cut must match the actual current nSigma gate.')
    protected.update(current['protectedFiles'])
    for path, candidate in candidate_objects.items():
        protected[str(path)] = digest(path)
        dest = snapshot/('calibration_'+str(window_paths.index(path))+'.json')
        shutil.copy2(path, dest)
        if read_json(dest) != candidate:
            raise ValueError('Candidate changed after it was read.')
        snapshots[str(dest)] = digest(dest)
    ylo, yhi = axis['dedxRange']
    half_width = finite(cfg['seed']['fitWindowSigmaHalfWidth'], True)
    alignment = finite(cfg['seed']['fitWindowAlignmentKeV'], True)
    if abs(alignment-axis['dedxBinWidth']*cfg['numeric']['displayDedxRebinFactor']) > tolerance:
        raise ValueError('Fit window alignment must equal exact display native-bin grouping.')
    rows, windows, row_audit = [], [], []
    candidate = candidate_objects[calibration_path]
    for index, rbin in enumerate(range(first, last+1)):
        lo, hi = rmin+(rbin-1)*step, rmin+rbin*step
        intervals, detail = [], []
        for node, (rigidity, weight) in enumerate(runtime.gl3_nodes(lo, hi)):
            for species_index, species in enumerate(SPECIES):
                value, policy = raw_parameters(runtime, current['tables'], candidate, curves,
                                               species_index, rigidity)
                mean, sigma, current_policy = runtime.current_parameters(current, species_index, rigidity)
                active = runtime.current_active(current, species_index, rigidity)
                rows.append([index, rbin, lo, hi, node, weight, rigidity, species, value['mean'],
                             value['sigma'], value['meanError'], value['sigmaError'], int(value['supported']),
                             mean, sigma, int(active), int(value['supported'])])
                detail.append(dict(node=node, species=species, candidatePolicy=policy,
                                   candidateSupported=bool(value['supported']), currentPolicy=current_policy,
                                   currentActive=active))
                # A common window includes all declared comparison functions, current
                # active selection, and the same unvalidated nuisance seed where absent.
                parameters = [(mean, sigma)] if active else []
                for other in candidate_objects.values():
                    shape, _ = raw_parameters(runtime, current['tables'], other, curves,
                                               species_index, rigidity)
                    parameters.append((shape['mean'], shape['sigma']))
                for center, width in parameters:
                    a, b = max(ylo, center-half_width*width), min(yhi, center+half_width*width)
                    if b > a:
                        intervals.append((a, b))
        if not intervals:
            raise ValueError('No observable nuclear/reference shape interval; do not invent a full-axis fit.')
        fitlo = max(ylo, ylo+math.floor((min(a for a, b in intervals)-ylo)/alignment)*alignment)
        fithi = min(yhi, ylo+math.ceil((max(b for a, b in intervals)-ylo)/alignment)*alignment)
        if fithi <= fitlo:
            raise ValueError('Empty whole-bin nuclear comparison window.')
        windows.append([index, rbin, fitlo, fithi])
        row_audit.append(dict(index=index, rBin=rbin, lo=lo, hi=hi, fitLo=fitlo, fitHi=fithi,
                              components=detail))
    numeric = dict(cfg['numeric'])
    for key, value in numeric.items():
        finite(value)
    numeric.update(sampleHistogram=args.sample, backgroundModel=args.background,
                   backgroundStartSlope=finite(cfg['allFit']['backgroundStartSlope']),
                   maxNSigma=finite(cfg['allFit']['maxNSigma'], True),
                   minimumPurityDenominator=finite(cfg['allFit']['minimumPurityDenominator']),
                   maximumConditionalChi2Ndf=finite(cfg['allFit']['maximumConditionalChi2Ndf'], True))
    shapes_path, windows_path = output/'shapes.tsv', output/'fit_windows.tsv'
    write_tsv(shapes_path, 'index rBin rLo rHi node weight rValue species mean sigma mean_error sigma_error supported current_mean current_sigma current_active candidate_active'.split(), rows)
    write_tsv(windows_path, 'index rBin fitLo fitHi'.split(), windows)
    write_tsv(output/'settings.tsv', None, [[k, v] for k, v in numeric.items()])
    save(output/'effective_settings.json', dict(sourceSettings=str(settings_path), configured=cfg,
         background=args.background, sample=args.sample, nativeRigidityRange=[minimum, maximum],
         comparisonCalibrationFiles=list(map(str, window_paths)), macroSettings=numeric))
    save(output/'shape_audit.json', dict(rows=row_audit, supportedAndReferenceExplicit=True,
         nuisancePolicy='Unsupported components use unvalidated raw legacy polynomials only as nuisance seeds; not adopted PID outside supported domain.',
         windowPolicy='Same current + all comparison-candidate/reference +-declared sigma union, native-axis clipped, whole-display-bin outward snapping.',
         uncertaintyPolicy='Conditional shape knots/errors only; no full TAG/ALL correlated or transfer uncertainty claim.'))
    readback = snapshot/'ValidateAllProduct.C'
    readback.write_text(READBACK_SOURCE)
    for path in (shapes_path, windows_path, output/'settings.tsv', output/'effective_settings.json',
                 output/'shape_audit.json', readback):
        snapshots[str(path)] = digest(path)
    context = dict(schemaVersion=1, inputStudy=str(inputs), inputCompletionSha256=digest(completion_path),
         cache=str(cache), cacheSha256=digest(cache), originalInputFingerprint=fingerprint(source),
         runtimeMainconf=str(mainconf), calibration=str(calibration_path),
         comparisonCalibrations=list(map(str, window_paths)), sourceSettings=str(settings_path),
         sampleHistogram=args.sample, backgroundModel=args.background, firstNativeRBin=first,
         nativeRows=last-first+1, productionEnabled=False, physicsValidationComplete=False,
         protected=protected, snapshots=snapshots,
         createdJst=datetime.now(timezone(timedelta(hours=9))).isoformat())
    save(output/'context.json', context)

    def preserve():
        if fingerprint(source) != context['originalInputFingerprint']:
            raise ValueError('Original ROOT fingerprint changed.')
        for path, expected in {**protected, **snapshots}.items():
            if digest(Path(path)) != expected:
                raise ValueError('Protected original or consumed snapshot changed: '+path)

    call = str(snapshot/'FitAllNuclearPurity.C')+'++('+','.join(json.dumps(str(p)) for p in
            (cache, output, output/'settings.tsv', shapes_path, windows_path))+')'
    command = ['bash', str(PROJECT/'script/singularity_root_batch.sh'), str(mainconf),
               'root4star', '-l', '-b', '-q', call]
    save(output/'command.json', command)
    preserve()
    if args.prepare_only:
        save(output/'prepared.json', dict(prepared=True, executionComplete=False, sourceInputPreserved=True,
                                         productionEnabled=False, nativeRows=last-first+1))
        print('NUCLEAR_ALL_PREPARED', output, flush=True)
        return
    print('NUCLEAR_ALL_RUN', output, flush=True)

    def execute(command, logfile, marker):
        try:
            with logfile.open('w') as log:
                run = subprocess.run(command, cwd=PROJECT, stdout=log, stderr=subprocess.STDOUT,
                                     timeout=timeout)
        finally:
            preserve()
        if run.returncode or marker not in logfile.read_text():
            raise RuntimeError('ROOT stage failed: '+str(logfile))

    try:
        execute(command, output/'fit.log', 'NUCLEAR_ALL_PURITY_COMPLETE')
        result = read_json(output/'all_fits.json')
        if result['productionEnabled'] is not False or result['sampleHistogram'] != args.sample \
                or len(result['slices']) != last-first+1:
            raise ValueError('Incomplete or unexpected ALL result identity.')
        for index, row in enumerate(result['slices']):
            expected = windows[index]
            if (row['index'], row['rBin'], row['fitLo'], row['fitHi']) != tuple(expected) \
                    or set(row['components']) != set(SPECIES):
                raise ValueError('ALL result native/window/species coverage mismatch.')
        validation_call = str(readback)+'+('+','.join((json.dumps(str(cache)),
             json.dumps(str(output/'all_fits.root')), json.dumps(args.sample), str(first),
             str(last-first+1), format(tolerance, '.17g'),
             format(cfg['numeric']['countVarianceRelativeTolerance'], '.17g')))+')'
        validation_command = command[:-1]+[validation_call]
        save(output/'readback_command.json', validation_command)
        execute(validation_command, output/'readback.log', 'NUCLEAR_ALL_READBACK_PASS')
        preserve()
        products = {str(path.relative_to(output)):digest(path) for path in output.rglob('*')
                    if path.is_file() and snapshot not in path.parents and path.name != 'completed.json'}
        save(output/'completed.json', dict(executionComplete=True, sourceInputPreserved=True,
             cacheReadbackPreserved=True, productClosedAndReadbackVerified=True,
             productionEnabled=False, physicsValidationComplete=False,
             numericalConverged=result['numericalConverged'], products=products,
             completedJst=datetime.now(timezone(timedelta(hours=9))).isoformat(),
             note='Execution completion is not shape identification, TOF transfer validation or adopted PID.'))
        print('NUCLEAR_ALL_COMPLETE', output, flush=True)
    except Exception as error:
        save(output/'failed.json', dict(executionComplete=False, productionEnabled=False,
                                       error=str(error), partialProductsPreserved=True))
        raise


if __name__ == '__main__':
    main()
