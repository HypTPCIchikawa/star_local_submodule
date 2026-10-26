#!/usr/bin/env python3
"""Post-hoc calibration-support diagnostics; never adopt PID, purity or cutoffs.

Chi-square levels are required reporting arguments, not physics acceptance cuts.
Precision and one-sided ROI coverage are reported without arbitrary rejection.
"""
import argparse
import csv
import math
from pathlib import Path
import shutil

from prepare_tag_study import owned, read, require, save, sha, tsv


SPECIES = ('d', 't', 'he3', 'he4')


def finite(value):
    return type(value) in (int, float) and math.isfinite(value)


def positive(value):
    return finite(value) and value > 0


def analyze_row(row, planned):
    for key in ('id', 'target', 'mode', 'rBin'):
        require(row.get(key) == planned[key], 'Fit/plan identity mismatch: ' + key)
    output = dict(id=planned['id'], species=planned['target'], mode=planned['mode'],
                  rBin=planned['rBin'], r=planned['r'],
                  fitLo=planned['fitLo'], fitHi=planned['fitHi'],
                  axisClipped=planned['fitWindowAxisClipped'],
                  fitCompleted=row.get('fitCompleted') is True,
                  numericalFlagsOK=False, phaseDiagnosticsAvailable=False,
                  phaseFlagsOK=None, foreignPeakInsideFit=None,
                  necessaryConditions={}, supportEligible=False,
                  chi2PerNdf=None, targetWindowYield=None, targetEffectiveCount=None,
                  originalMean=None, originalSigma=None, mean=None, sigma=None,
                  meanError=None, sigmaError=None, meanErrorOverSigma=None,
                  sigmaErrorOverSigma=None, leftCoverageSigma=None, rightCoverageSigma=None,
                  zeroTargetYield=False, zeroTargetEffectiveCount=False,
                  targetModeInsideFit=False,
                  descriptiveOnly=True, curveAdoptionAllowed=False)
    if not output['fitCompleted']:
        output['necessaryConditions']['fitCompleted'] = False
        output['error'] = row.get('error', 'unrecorded')
        return output
    require(row.get('curveAdoptionAllowed') is False and row.get('purityAvailable') is False,
            'Support assessment may only read non-adopting diagnostic fits')
    require(isinstance(row.get('components'), list) and
            len([c for c in row['components'] if c.get('species') == row['target']]) == 1,
            'Require exactly one saved target component')
    component = next(c for c in row['components'] if c['species'] == row['target'])
    require(row['columnNames'].count(row['target']) == 1, 'Ambiguous target-yield column')
    column = row['columnNames'].index(row['target'])
    target_yield = row['columnWindowYields'][column]
    # Effective counts are component ordered, not background-column ordered.
    component_index = next(i for i, c in enumerate(row['components'])
                           if c['species'] == row['target'])
    effective = row['effectiveCount'][component_index]
    output['numericalFlagsOK'] = (type(row.get('status')) is int and row['status'] == 0 and
        all(row.get(key) is True for key in ('finite', 'pointAccepted', 'amplitudeConverged',
                                            'yieldCovarianceValid', 'shapeCovarianceValid')) and
        row.get('atBoundary') is False)
    available = row.get('responseOptimizerDiagnosticsAvailable')
    require(available in (None, True, False), 'Malformed phase availability flag')
    output['phaseDiagnosticsAvailable'] = available is True
    if available is True:
        output['phaseFlagsOK'] = (
            all(row.get(key) is True for key in ('responseMinimizeReturned', 'responseMinimizeConverged',
                    'responseHesseAttempted', 'responseHesseSucceeded')) and
            row.get('responseMigradStatus') == 0 and row.get('responseHesseStatus') == 0 and
            finite(row.get('responseMigradEdm')) and row['responseMigradEdm'] >= 0)
    output['foreignPeakInsideFit'] = row.get('foreignPeakInsideFit')
    output['targetModeInsideFit'] = (finite(component.get('mean')) and
                                    planned['fitLo'] < component['mean'] < planned['fitHi'])
    output['necessaryConditions'] = dict(
        fitCompleted=True, numericalFlags=output['numericalFlagsOK'],
        optimizerPhasesIfAvailable=available is not True or output['phaseFlagsOK'],
        noForeignPeakInsideFit=row.get('foreignPeakInsideFit') is False,
        positiveFiniteTargetYield=positive(target_yield),
        positiveFiniteTargetEffectiveCount=positive(effective),
        positiveFiniteTargetMean=positive(component.get('mean')),
        positiveFiniteTargetWidth=positive(component.get('sigma')),
        positiveFiniteMeanError=positive(component.get('meanError')),
        positiveFiniteWidthError=positive(component.get('sigmaError')),
        targetModeInsideFit=output['targetModeInsideFit'],
        finiteNonnegativeChi2=finite(row.get('pearsonPerNdf')) and row['pearsonPerNdf'] >= 0)
    output['supportEligible'] = all(output['necessaryConditions'].values())
    output['zeroTargetYield'] = finite(target_yield) and target_yield == 0
    output['zeroTargetEffectiveCount'] = finite(effective) and effective == 0
    for value, key in ((target_yield, 'targetWindowYield'), (effective, 'targetEffectiveCount'),
                       (row.get('pearsonPerNdf'), 'chi2PerNdf')):
        output[key] = value if finite(value) else None
    for source, target in (('initialMean', 'originalMean'), ('initialSigma', 'originalSigma'),
                           ('mean', 'mean'), ('sigma', 'sigma'),
                           ('meanError', 'meanError'), ('sigmaError', 'sigmaError')):
        value = component.get(source)
        output[target] = value if finite(value) else None
    sigma = component.get('sigma')
    if positive(sigma):
        if finite(output['meanError']):
            output['meanErrorOverSigma'] = output['meanError'] / sigma
        if finite(output['sigmaError']):
            output['sigmaErrorOverSigma'] = output['sigmaError'] / sigma
    left = row.get('responseSigmaLeft', sigma)
    right = row.get('responseSigmaRight', sigma)
    if finite(output['mean']):
        if positive(left):
            output['leftCoverageSigma'] = (output['mean'] - planned['fitLo']) / left
        if positive(right):
            output['rightCoverageSigma'] = (planned['fitHi'] - output['mean']) / right
    return output


def contiguous_runs(rows, eligible_ids):
    """Consecutive sampled-plan points, not proof every R in the span is validated."""
    result, current = [], []
    for row in sorted(rows, key=lambda r: r['rBin']):
        if row['id'] in eligible_ids:
            current.append(row)
        elif current:
            result.append(current)
            current = []
    if current:
        result.append(current)
    return [dict(firstR=run[0]['r'], lastR=run[-1]['r'], sampledPointCount=len(run),
                 ids=[r['id'] for r in run],
                 unsampledInteriorsAreValidated=False) for run in result]


def summarize(rows, levels):
    require(levels and len(set(levels)) == len(levels) and
            all(positive(v) for v in levels), 'Require unique positive descriptive chi2 levels')
    result = {}
    for species in SPECIES:
        sample = [r for r in rows if r['species'] == species]
        thresholds = []
        for level in levels:
            selected = [r for r in sample if r['supportEligible'] and r['chi2PerNdf'] <= level]
            thresholds.append(dict(chi2PerNdfAtMost=level, selectedCount=len(selected),
                lowRCount=sum(r['r'] < 2 for r in selected),
                highRCount=sum(r['r'] >= 2 for r in selected),
                ids=[r['id'] for r in selected],
                consecutiveSampledRuns=contiguous_runs(sample, {r['id'] for r in selected}),
                enoughPointsForFourMeanParameters=len(selected) >= 4,
                enoughPointsForThreeWidthParameters=len(selected) >= 3,
                pointCountAloneDoesNotValidatePolynomial=True,
                scientificAcceptanceCriterion=False))
        result[species] = dict(total=len(sample), numericalFlagsOK=sum(r['numericalFlagsOK'] for r in sample),
            supportEligible=sum(r['supportEligible'] for r in sample),
            zeroTargetYieldIds=[r['id'] for r in sample if r['zeroTargetYield']],
            zeroTargetEffectiveCountIds=[r['id'] for r in sample if r['zeroTargetEffectiveCount']],
            targetModeOutsideFitIds=[r['id'] for r in sample
                                     if r['fitCompleted'] and not r['targetModeInsideFit']],
            foreignPeakInsideFitCount=sum(r['foreignPeakInsideFit'] is True for r in sample),
            thresholds=thresholds)
    return result


def run(fit_dir, output, levels):
    require(not output.exists() and output.parent == fit_dir,
            'Require a NEW direct child quality output under the selected fit directory')
    sources = (fit_dir / 'plan.json', fit_dir / 'fits.json',
               fit_dir / 'fit.log.completed.json', Path(__file__).resolve())
    protected = {str(p): sha(p) for p in sources}
    plan, fits, completion = (read(p) for p in sources[:3])
    require('recenterSurvey' not in plan,
            'Unimplemented full-survey ROI assessment is not supported by this reporting version')
    require(completion.get('status') == 'completed' and completion.get('inputsUnchanged') is True and
            fits.get('status') == 'completed' and fits.get('rawReadback') is True and
            fits.get('predictionReadback') is True and fits.get('productionEnabled') is False and
            fits.get('purityAvailable') is False and fits.get('curveAdoptionAllowed') is False,
            'An incomplete or production-adopting fit cannot be assessed')
    require(fits['totalSlices'] == len(fits['results']) == len(plan['slices']) == plan['totalFits'],
            'Fit/plan row-count mismatch')
    rows = [analyze_row(row, planned) for row, planned in zip(fits['results'], plan['slices'])]
    summary = summarize(rows, levels)
    require(all(sha(Path(path)) == value for path, value in protected.items()),
            'Protected support inputs changed while reading')
    output.mkdir()
    shutil.copy2(str(sources[3]), str(output / sources[3].name))
    save(output / 'quality.json', dict(schemaVersion=1, sourceFit=str(fit_dir),
        model=fits.get('model', 'unrecorded'),
        responseMeanMeaning=fits.get('responseMeanMeaning', 'Gaussian location/mean'),
        responseSigmaMeaning=fits.get('responseSigmaMeaning', 'Gaussian width'),
        descriptiveChi2Levels=levels, thresholdsAreAcceptanceCriteria=False,
        necessaryConditionsAreNotSufficientForPhysics=True, rows=rows, bySpecies=summary,
        productionEnabled=False, coefficientsUpdated=False, purityAvailable=False,
        curveAdoptionAllowed=False, isolatedVetoesCreated=False,
        limitations=['Low chi2 with zero target yield does not identify a peak.',
            'Finite covariance and positive error flags alone do not identify a peak.',
            'Precision and asymmetric ROI coverage are diagnostics, not new implicit cuts.',
            'A target mode outside its fitted interval cannot directly calibrate that peak.',
            'Consecutive support refers to sampled R points only; gaps and extrapolation are not validated.',
            'No polynomial coefficients are fitted or adopted by this tool.']))
    columns = ['id', 'species', 'mode', 'rBin', 'r', 'numericalFlagsOK',
        'phaseDiagnosticsAvailable', 'phaseFlagsOK', 'foreignPeakInsideFit', 'supportEligible',
        'chi2PerNdf', 'targetWindowYield', 'targetEffectiveCount', 'zeroTargetYield',
        'zeroTargetEffectiveCount', 'originalMean', 'mean', 'meanError',
        'originalSigma', 'sigma', 'sigmaError', 'meanErrorOverSigma', 'sigmaErrorOverSigma',
        'fitLo', 'fitHi', 'axisClipped', 'leftCoverageSigma', 'rightCoverageSigma', 'targetModeInsideFit']
    tsv(output / 'quality.tsv', columns, ([row.get(key) for key in columns] for row in rows))
    require(all(sha(Path(path)) == value for path, value in protected.items()),
            'Protected support inputs changed during output')
    save(output / 'completed.json', dict(status='completed', protectedInputs=protected,
        inputsUnchanged=True, fitResultsModified=False, rows=len(rows),
        products={name: sha(output / name) for name in ('quality.json', 'quality.tsv', sources[3].name)},
        productionEnabled=False, coefficientsUpdated=False, purityAvailable=False))
    for species, values in summary.items():
        print(species, 'n=' + str(values['total']), 'identified=' + str(values['supportEligible']),
              'threshold_counts=' + str([q['selectedCount'] for q in values['thresholds']]), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fit-dir', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--chi2-levels', required=True, nargs='+', type=float)
    args = parser.parse_args()
    run(owned(args.fit_dir), owned(args.output_dir, False), args.chi2_levels)
