#!/usr/bin/env python3
"""Replay the read-only saved-histogram audit; never fit or overwrite products."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import shutil
import tarfile

import yaml

JOBID = '0BEFC64E62F19BBA534A129CABC9B106'
ANA = 'auau3p85_anaFemtoLambda_4He_KFParticle'
MEMBERS = (
    'StMaker/common/FemtoLambdaLegacy.cxx',
    'StMaker/common/FemtoLambdaLegacy.h',
    'config/cuts/nuclearid/nuclearid_femtoLambda_legacy.yaml',
    'config/hist/hist_femtoLambda_he4_nuclear.yaml',
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def read(path):
    return json.loads(path.read_text())


def create(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('xb') as stream:
        stream.write(data)


def save(path, value):
    create(path, (json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n').encode())


def snapshot(path, data):
    if path.exists():
        require(path.is_file() and path.read_bytes() == data,
                'Existing snapshot differs; refusing overwrite: '+str(path))
    else:
        create(path, data)


def fingerprint(path):
    stat = path.stat()
    return dict(path=str(path), size=stat.st_size, mtime_ns=stat.st_mtime_ns,
                device=stat.st_dev, inode=stat.st_ino)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    project, output = args.project.resolve(strict=True), args.output.resolve()
    if args.verify:
        done = read(output/'completed.json')
        require(all(sha(output/p) == expected for p, expected in done['products'].items()),
                'Saved audit product hash changed')
        require(all(sha(Path(p)) == expected for p, expected in done['originalInputs'].items()),
                'Small original input changed')
        audit = read(output/'audit.json')
        require(fingerprint(Path(audit['originalRoot']['fingerprint']['path'])) ==
                audit['originalRoot']['fingerprint'], 'Large ROOT fingerprint changed')
        print('SAVED_HIGH_R_INPUT_AUDIT_VERIFY_PASS')
        return
    require(not (output/'audit.json').exists() and not (output/'completed.json').exists(),
            'Existing audit cannot be replaced; choose a new output or use --verify')
    output.mkdir(parents=True, exist_ok=True)
    require(set(p.name for p in output.iterdir()) <=
            {'README.md', 'replay_audit.py', 'snapshot', 'attempt_001_failure.md'},
            'Output already contains products')
    base = project/'rootfile'/ANA
    run = project/'job/run/merge_notify'/(ANA+'_'+JOBID)
    runmeta = project/'job/run/runmeta'/('runmeta_'+ANA+'_'+JOBID+'.json')
    meta = read(runmeta)
    bundle = Path(meta['runtime']['runtimeBundle']['path'])
    require(sha(bundle) == meta['runtime']['runtimeBundle']['sha256'],
            'Submit-time runtime bundle hash mismatch')
    accepted = run/'manual_acceptance_20261003/result.json'
    acceptance = read(accepted)
    source = base/(ANA+'_2021P24iy_partial_merge.root')
    cache_dir = base/'pdf/nuclear_pid_tof_20261014/inputs_001'
    cache_context, cache_done = read(cache_dir/'context.json'), read(cache_dir/'completed.json')
    require(fingerprint(source) == cache_context['inputFingerprint'],
            'Original ROOT does not match existing native cache provenance')
    for name in ('context.json', 'native_audit.json', 'cache.root'):
        require(sha(cache_dir/name) == cache_done['products'][name], 'Cache evidence changed: '+name)
    snapshot_paths = []
    old_text = {}
    with tarfile.open(str(bundle), 'r:gz') as archive:
        for member in MEMBERS:
            item = archive.getmember(member)
            require(item.isfile(), 'Runtime member must be a regular file')
            data = archive.extractfile(item).read()
            path = output/'snapshot/submitted_runtime'/member
            snapshot(path, data)
            snapshot_paths.append(path)
            old_text[member] = data.decode()
    original_inputs = [runmeta, bundle, accepted, run/'merged_audit.tsv',
                       cache_dir/'context.json', cache_dir/'completed.json',
                       cache_dir/'native_audit.json', cache_dir/'effective_settings.json',
                       cache_dir/'source/ExtractTofInputs.C', cache_dir/'source/prepare_inputs.py',
                       project/'config/cuts/nuclearid/nuclearid_femtoLambda_auau3p85_purity50_20261015.yaml',
                       project/'StMaker/common/FemtoLambdaLegacy.cxx',
                       project/'StMaker/common/NuclearIdDeDxVsMom.h',
                       project/'config/calibration/dedx_helium_auau3p85_2021P24iy_from040.yaml']
    for original in original_inputs:
        require(original.is_file(), 'Missing audit input '+str(original))
        dest = output/'snapshot/project'/original.relative_to(project)
        snapshot(dest, original.read_bytes())
        snapshot_paths.append(dest)
    small_hashes = {str(p): sha(p) for p in original_inputs}
    small_hashes[str(cache_dir/'cache.root')] = sha(cache_dir/'cache.root')
    old = yaml.safe_load(old_text['config/cuts/nuclearid/nuclearid_femtoLambda_legacy.yaml'])
    final_path = project/'config/cuts/nuclearid/nuclearid_femtoLambda_auau3p85_purity50_20261015.yaml'
    final = yaml.safe_load(final_path.read_text())
    old_cxx = old_text['StMaker/common/FemtoLambdaLegacy.cxx']
    require('QA population is independent of final best-species / p-over-q selection.' in old_cxx,
            'Submitted-source QA independence comment not found; audit manually')
    expected_filling = 'if(m2>mTofMean[s]-mM2SigmaCut*mTofSigma[s] && m2<mTofMean[s]+mM2SigmaCut*mTofSigma[s])'
    require(expected_filling in old_cxx, 'Submitted-source tag fill expression changed')
    merged_lines = (run/'merged_audit.tsv').read_text().splitlines()
    # Original C++ audit writes real tabs; accept the escaped display form too.
    records = [line.replace('\\t', '\t').split('\t') for line in merged_lines]
    hist_rows = [row for row in records if row[0] == 'HIST']
    require(len(hist_rows) == acceptance['histogramCount'], 'Unexpected merged histogram count')
    names = [row[1].replace(';1', '') for row in hist_rows]
    require(not any(row[2].startswith('TH3') for row in hist_rows), 'Unexpected joint TH3: inspect manually')
    required = ['true/hDedxP', 'true/hDedxP_d_m2', 'true/hDedxP_t_m2',
                'true/hDedxP_3He_m2', 'true/hM2P']
    require(all(name in names for name in required), 'Required QA histogram missing')
    absent = ['true/hDedxP_4He_TOFControl', 'true/hNSigmaHe4VsP_TOFControl',
              'true/hM2P_4He_TPC', 'true/hDedxP_4He_selected', 'true/hM2P_4He_selected']
    require(not any(name in names for name in absent), 'He4 audit assumption no longer valid')
    windows = []
    for species, key, cache_name in (('d','Deuteron','TAG_D'), ('t','Triton','TAG_T'),
                                     ('he3','He3','TAG_HE3')):
        mean, sigma = old['nuclearTofMean'+key], old['nuclearTofSigma'+key]
        lo, hi = mean-old['m2SigmaCut']*sigma, mean+old['m2SigmaCut']*sigma
        cached_name = {'d':'D','t':'T','he3':'HE3'}[species]
        settings = cache_context['settings']
        configured = final['highRigidityTof']['species'][species]
        require(abs(lo-settings['tag'+cached_name+'Mass2Min']) < 1.e-12 and
                abs(hi-settings['tag'+cached_name+'Mass2Max']) < 1.e-12,
                'Submitted window does not match existing cache metadata')
        require(abs(lo-configured['mass2Min']) < 1.e-12 and abs(hi-configured['mass2Max']) < 1.e-12,
                'Configured high-R m2 window differs from saved sample')
        windows.append(dict(species=species, sourceHistogram={'d':'true/hDedxP_d_m2',
                            't':'true/hDedxP_t_m2','he3':'true/hDedxP_3He_m2'}[species],
                            cacheKey=cache_name, savedMass2Min=lo, savedMass2Max=hi,
                            strictBounds=True, finalConfigMatched=True,
                            finalRigidityRange=[final['highRigidityTof']['thresholdPOverQ'],
                                                final['maxPOverQ']], finalUpperBoundExclusive=True))
    save(output/'audit.json', dict(
        schemaVersion=1, createdAt=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        jobid=JOBID, dataset='production_3p85GeV_fixedTarget_2021 / P24iy',
        originalRoot=dict(fingerprint=fingerprint(source), wholeFileHashCalculated=False,
                          priorAcceptanceSha256=acceptance['mergedSha256'],
                          inputSubjobs=acceptance['inputCount'], partial=acceptance['partial'],
                          excludedIndices=acceptance['excludedIndices']),
        submitRuntime=dict(bundleSha256=small_hashes[str(bundle)],
                           rootVersion=meta['runtime']['rootVersion'],
                           libraryTag=meta['runtime']['libraryTag'],
                           starHostSys=meta['runtime']['starHostSys'],
                           dirtyWorktreeRecorded=meta['git']['isDirty']),
        cache=dict(path=str(cache_dir/'cache.root'), sha256=small_hashes[str(cache_dir/'cache.root')],
                   priorManifestHashesVerified=True, originalFingerprintMatched=True,
                   nativeAudit=read(cache_dir/'native_audit.json'), tagWindows=windows),
        originalHistograms=names, absentHe4JointControls=absent,
        noTH3InOriginalHistogramInventory=True,
        sampling=dict(trackQuality='nHitsDedx >= 15 and gPt >= 0.1; dedx > 0',
                      rigidity='StPicoTrack::pMom().Mag(), raw measured |p/q|; no He charge correction',
                      tof='isTofTrack, valid traits index, beta > 0',
                      mass2='R^2 * (1/beta^2 - 1), raw (m/q)^2',
                      qaBeforeFinalTPCBestSpeciesAndMaxR=True, tagIsPureIsotopeTruth=False,
                      tagsAreSubsetsOfAll=True, combineAsIndependentSamplesAllowed=False),
        he4=dict(finalMass2Window=final['highRigidityTof']['species']['he4'],
                 exactWindowAvailable=False, exactFinalPurityAssessableFromThisRoot=False,
                 narrowerTagDAvailable=True, tagDMayReplaceFinalWindow=False,
                 sourceOfLimitation='No preTPC He4 3<m2<4 dEdx histogram, TH3 or dEdx-m2 joint saved',
                 extensionConfigured=final['highRigidityTof']['extendHe4BeyondMaxPOverQ'],
                 nativeAxisUpper=6, overflowMustBeSeparate=True),
        interpretation=[
            'Existing prePID QA can be re-gated with current TPC response; old selected nuclear histogram is not the new selected sample.',
            'Raw TPC projection purity after m2 requires fitting isotope components; TOF matching does not imply purity=1.',
            'Finite native R and dEdx bins require explicit gate-integration approximation.',
            'QA is charge-inclusive and centrality-integrated, not CF-pair-, centrality- or mT-conditioned purity.',
            'The recorded partial merge does not prove full 23-billion-event coverage.',
            'No fit, PID calibration change, production update, batch submission, hadd or mail in this audit.']))
    create(output/'SHA256SUMS', ''.join(sha(p)+'  '+str(p.relative_to(output))+'\n'
                                      for p in sorted(snapshot_paths)).encode())
    products = {str(p.relative_to(output)): sha(p) for p in output.rglob('*')
                if p.is_file() and p.name != 'completed.json'}
    save(output/'completed.json', dict(executionComplete=True, noFit=True,
         productionEnabled=False, originalRootNotReadOrHashed=True,
         products=products, originalInputs=small_hashes))
    print('SAVED_HIGH_R_INPUT_AUDIT_COMPLETE histograms='+str(len(names)))


if __name__ == '__main__':
    main()
