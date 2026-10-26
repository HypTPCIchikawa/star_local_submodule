#!/usr/bin/env python3
"""Compare immutable Gaussian/log-normal TAG predictions; never refit or adopt PID."""
import argparse
from datetime import datetime, timedelta, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess

import yaml

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
PRIMARY = {'TAG_D': ('d', 'he4'), 'TAG_T': ('t',), 'TAG_HE3': ('he3',)}


def read_json(path):
    def invalid(value):
        raise ValueError('Nonfinite JSON constant '+value)
    return json.loads(Path(path).read_text(), parse_constant=invalid)


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            value.update(block)
    return value.hexdigest()


def fingerprint(path):
    value = Path(path).stat()
    return dict(path=str(path), size=value.st_size, mtime_ns=value.st_mtime_ns,
                device=value.st_dev, inode=value.st_ino)


def save(path, obj):
    if path.exists():
        raise ValueError('Existing output will not be overwritten: '+str(path))
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def tsv(path, rows):
    if path.exists():
        raise ValueError('Existing output will not be overwritten: '+str(path))
    def encode(value):
        if isinstance(value, float):
            if not math.isfinite(value):
                raise ValueError('Nonfinite TSV data.')
            return format(value, '.17g')
        return str(value)
    path.write_text(''.join('\t'.join(map(encode, row))+'\n' for row in rows))


def project_path(path):
    path = Path(path).resolve(strict=True)
    if PROJECT not in path.parents:
        raise ValueError('Expected project-owned input/config/source: '+str(path))
    return path


def ratio(numerator, denominator):
    return numerator/denominator if denominator > 0 else None


def check_study(study, response, background, mainconf, protected):
    context, completed = (read_json(study/name) for name in ('context.json', 'completed.json'))
    if not completed.get('executionComplete') or completed.get('productionEnabled') is not False \
            or not completed.get('inputFingerprintPreserved') or not completed.get('protectedSourcesPreserved'):
        raise ValueError('Incomplete or unprotected source study: '+str(study))
    if context['runtimeMainconf'] != str(mainconf) or context['background'] != background:
        raise ValueError('Different declared runtime/background study.')
    if response == 'gaussian':
        if context['selectedPartition'] != 'native':
            raise ValueError('Gaussian source must be the completed native study.')
    elif context.get('responseModel') != 'lognormal' \
            or context.get('responseParameterization') != 'mode_and_local_curvature' \
            or context.get('curveAdoptionAllowed') is not False:
        raise ValueError('Log-normal source lacks explicit response meaning/adoption guards.')
    for name, expected in completed['products'].items():
        path = study/name
        if digest(path) != expected:
            raise ValueError('Completed source artifact changed: '+str(path))
        protected[str(path)] = expected
    for name in ('context.json', 'completed.json'):
        protected[str(study/name)] = digest(study/name)
    for block in ('protected', 'snapshots'):
        for path, expected in context[block].items():
            if digest(path) != expected:
                raise ValueError('Source/cache/snapshot changed: '+path)
            protected[path] = expected
    original = Path(context['inputFingerprint']['path']).resolve(strict=True)
    if '2021P24iy_partial_merge' not in original.name or fingerprint(original) != context['inputFingerprint']:
        raise ValueError('Source 2021 P24iy ROOT fingerprint changed.')
    fit = read_json(study/'tag_fits.json')
    if len(fit['slices']) != completed['nSlices']:
        raise ValueError('Saved slice coverage disagrees with completion manifest.')
    lookup = {}
    for row in fit['slices']:
        if row['sample'] not in PRIMARY or row['firstNativeRBin'] != row['lastNativeRBin'] \
                or abs(row['hi']-row['lo']-.01) > 1.e-9:
            raise ValueError('Expected single non-overlapping native .01 rigidity cells.')
        key = (row['sample'], row['firstNativeRBin'])
        if key in lookup:
            raise ValueError('Repeated native cell within one TAG spectrum.')
        lookup[key] = row
    return context, completed, lookup


def table_pair(gaussian, lognormal, background, coverage, tolerance):
    same = ('firstNativeRBin', 'lastNativeRBin', 'firstFitDedxBin', 'lastFitDedxBin', 'fitCount')
    boundaries = ('lo', 'hi', 'r', 'fitLo', 'fitHi')
    keys = set(lognormal)
    if coverage == 'full' and keys != set(gaussian):
        raise ValueError('Full comparison requires identical native TAG coverage.')
    if not keys or not keys <= set(gaussian):
        raise ValueError('Log-normal source does not match Gaussian native cells.')
    rows = []
    for key in sorted(keys):
        g, l = gaussian[key], lognormal[key]
        for name in same:
            if g[name] != l[name]:
                raise ValueError('Different native bins/count for '+str((key, name)))
        for name in boundaries:
            if abs(g[name]-l[name]) > tolerance:
                raise ValueError('Different native/fit geometry: '+str((key, name)))
        for species in PRIMARY[key[0]]:
            row = dict(background=background, sample=key[0], species=species,
                       rBin=key[1], r=g['r'], rLo=g['lo'], rHi=g['hi'], fitLo=g['fitLo'],
                       fitHi=g['fitHi'], fitCount=g['fitCount'], gaussianIndex=g['index'],
                       lognormalIndex=l['index'], deltaD=l['dataDeviance']-g['dataDeviance'],
                       deltaDIsSignificance=0)
            for prefix, source in (('gaussian', g), ('lognormal', l)):
                component = source['components'][species]
                flags = []
                if not source['numericalConverged']:
                    flags.append('numerical_invalid')
                if component['zeroYield']:
                    flags.append('zero_yield')
                if component['centerOutsideDetector']:
                    flags.append('outside_detector')
                if source['atBoundary']:
                    flags.append('parameter_boundary')
                if not component['identityValidated']:
                    flags.append('identity_unvalidated')
                row.update({prefix+'Center': component['mean'], prefix+'Width': component['sigma'],
                    prefix+'CenterError': component['meanError'], prefix+'WidthError': component['sigmaError'],
                    prefix+'WindowYield': component['windowYield'], prefix+'EffectiveCount': component['effectiveCount'],
                    prefix+'DetectorAcceptance': component['detectorAcceptance'], prefix+'WindowAcceptance': component['windowAcceptance'],
                    prefix+'NumericalValid': int(source['numericalConverged']), prefix+'ZeroYield': int(component['zeroYield']),
                    prefix+'OutsideDetector': int(component['centerOutsideDetector']), prefix+'Boundary': int(source['atBoundary']),
                    prefix+'IdentityValidated': int(component['identityValidated']), prefix+'Status': source['status'],
                    prefix+'CovarianceStatus': source['covarianceStatus'], prefix+'EDM': source['edm'],
                    prefix+'AmpKkt': source['amplitudeKkt'], prefix+'D': source['dataDeviance'],
                    prefix+'Pearson': source['chi2'], prefix+'ConditionalNdf': source['conditionalNdf'],
                    prefix+'DPerCount': source['dataDeviance']/max(1.,source['fitCount']),
                    prefix+'BgSlope': source['backgroundSlope'], prefix+'Warnings': ','.join(flags) or 'none'})
            rows.append(row)
    return rows, keys


def summaries(rows):
    result = []
    for background in ('linear', 'exponential'):
        for sample in PRIMARY:
            # TAG_D d/He4 are one dataset. Count it only once here.
            selected = [r for r in rows if r['background'] == background and r['sample'] == sample
                        and r['species'] == PRIMARY[sample][0]]
            if not selected:
                continue
            record = dict(background=background, sample=sample, nativeSlices=len(selected),
                          fitWindowCount=sum(r['fitCount'] for r in selected))
            for prefix in ('gaussian', 'lognormal'):
                d = sum(r[prefix+'D'] for r in selected)
                chi2 = sum(r[prefix+'Pearson'] for r in selected)
                ndf = sum(r[prefix+'ConditionalNdf'] for r in selected)
                record.update({prefix+'D': d, prefix+'DPerCount': ratio(d,record['fitWindowCount']),
                               prefix+'PearsonPerConditionalNdf': ratio(chi2,ndf),
                               prefix+'NumericallyConvergedSlices': sum(r[prefix+'NumericalValid'] for r in selected),
                               prefix+'BoundarySlices': sum(r[prefix+'Boundary'] for r in selected)})
            record['lognormalOverGaussianDeviance'] = ratio(record['lognormalD'],record['gaussianD'])
            record['lognormalLowerDevianceSlices'] = sum(r['deltaD'] < 0 for r in selected)
            result.append(record)
    return result


def validate_pdf(path, expected):
    data = path.read_bytes()
    if not data.startswith(b'%PDF-') or b'%%EOF' not in data[-1024:]:
        raise ValueError('Incomplete PDF: '+str(path))
    pages = len(re.findall(rb'/Type\s*/Page\b',data))
    boxes = {}
    for box in (b'MediaBox',b'CropBox'):
        values = re.findall(rb'/'+box+rb'\s*\[([^]]+)\]',data)
        boxes[box.decode()] = [list(map(float,v.split())) for v in values]
    if pages != expected or not boxes['MediaBox']:
        raise ValueError('PDF page count/media box mismatch.')
    media = boxes['MediaBox'][0]
    for crop in boxes['CropBox']:
        if len(crop) != 4 or crop[0] < media[0] or crop[1] < media[1] \
                or crop[2] > media[2] or crop[3] > media[3]:
            raise ValueError('ROOT PDF CropBox extends beyond MediaBox: '+str(crop))
    return dict(pages=pages, mediaBoxes=boxes['MediaBox'], cropBoxes=boxes['CropBox'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for response in ('gaussian', 'lognormal'):
        for background in ('linear', 'exponential'):
            parser.add_argument('--'+response+'-'+background, type=Path, required=True)
    parser.add_argument('--coverage', choices=('full','explicit_subset'), required=True,
                        help='Subset comparison is only allowed when explicitly selected, never automatic fallback.')
    parser.add_argument('--mainconf', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--representative-rigidity', type=float, action='append', required=True)
    parser.add_argument('--display-rebin', type=int, required=True)
    parser.add_argument('--roi-half-width', type=float, required=True)
    parser.add_argument('--center-display-max', type=float, required=True)
    parser.add_argument('--width-display-max', type=float, required=True)
    parser.add_argument('--timeout-seconds', type=int, required=True)
    parser.add_argument('--prepare-only', action='store_true')
    args = parser.parse_args()
    if args.display_rebin <= 0 or args.timeout_seconds <= 0 or len(set(args.representative_rigidity)) != len(args.representative_rigidity) \
            or not all(math.isfinite(v) and v > 0 for v in args.representative_rigidity+
                       [args.roi_half_width,args.center_display_max,args.width_display_max]):
        raise ValueError('Invalid explicit display-only parameters.')
    mainconf = project_path(args.mainconf)
    config = yaml.safe_load(mainconf.read_text())
    analysis = project_path(PROJECT/'config'/config['analysis'])
    tag = yaml.safe_load(analysis.read_text())['starTag']
    if (tag['collision'],tag['year'],tag['productionTag'],tag['libraryTag']) != ('auau3p85','year2021','P24iy','SL24y'):
        raise ValueError('Only 3p85 2021 P24iy SL24y runtime is permitted.')
    protected = {str(p):digest(p) for p in (Path(__file__).resolve(),HERE/'CompareNuclearResponses.C',
                  mainconf,analysis,PROJECT/'script/singularity_root_batch.sh')}
    studies, sources = {}, {}
    for response in ('gaussian','lognormal'):
        for background in ('linear','exponential'):
            study = project_path(getattr(args,response+'_'+background))
            key = (response,background)
            studies[key] = study
            sources[key] = check_study(study,response,background,mainconf,protected)
    if len(set(studies.values())) != 4 or len({p.parent for p in studies.values()}) != 1:
        raise ValueError('Choose four distinct sibling completed studies.')
    inputs = [value[0]['inputFingerprint'] for value in sources.values()]
    if any(value != inputs[0] for value in inputs):
        raise ValueError('Response studies have different original ROOT inputs.')
    if len({value[0]['inputStudy'] for value in sources.values()}) != 1:
        raise ValueError('Response studies have different immutable input caches.')
    output = args.output_dir.resolve()
    if output.exists() or output.parent != next(iter(studies.values())).parent:
        raise ValueError('Output must be a new sibling comparison directory.')
    rows, model_keys = [], []
    for background in ('linear','exponential'):
        records, keys = table_pair(sources['gaussian',background][2],sources['lognormal',background][2],
                                  background,args.coverage,1.e-9)
        rows.extend(records);model_keys.append(keys)
    if model_keys[0] != model_keys[1]:
        raise ValueError('Background comparisons have different native TAG coverage.')
    rlo, rhi = min(r['rLo'] for r in rows),max(r['rHi'] for r in rows)
    for representative in args.representative_rigidity:
        if any(not any(r['species'] == sp and r['background'] == bg and abs(r['r']-representative) < 1.e-8
                       for r in rows) for sp in ('d','t','he3','he4') for bg in ('linear','exponential')):
            raise ValueError('Representative does not exist as a matched native center: '+str(representative))
    output.mkdir();snapshot=output/'source';snapshot.mkdir();snapshots={}
    for src in (Path(__file__).resolve(),HERE/'CompareNuclearResponses.C',mainconf,analysis,PROJECT/'script/singularity_root_batch.sh'):
        dest=snapshot/src.name;shutil.copy2(src,dest)
        if digest(src)!=digest(dest):
            raise ValueError('Source changed while copying snapshot.')
        snapshots[str(dest)]=digest(dest)
    settings=dict(coverage=args.coverage,representativeRigidities=args.representative_rigidity,
        displayRebin=args.display_rebin,roiHalfWidth=args.roi_half_width,centerDisplayMax=args.center_display_max,
        widthDisplayMax=args.width_display_max,axisRange=[0.,120.],paperCentimeters=[26.,20.],
        timeoutSeconds=args.timeout_seconds,productionEnabled=False,noFit=True,noGateChanges=True,
        gaussianCenterMeaning='arithmetic mean',lognormalCenterMeaning='mode',
        gaussianWidthMeaning='standard deviation of native-R projected Gaussian',
        lognormalWidthMeaning='local curvature scale of native-R projected log-normal, NOT standard deviation',
        gaussianLognormalAreNonNested=True,deltaDIsSignificance=False,
        tagDHe4AreNotIndependentDatasets=True,identitiesValidated=False,
        warningsPolicy='All unvalidated identities are red points; no joining unsupported rows. Display-clipped values remain exact in table.')
    save(output/'comparison_settings.json',settings)
    columns=list(rows[0]);tsv(output/'comparison.tsv',[columns]+[[r[k] for k in columns] for r in rows])
    summaries_data=summaries(rows);save(output/'comparison.json',dict(settings=settings,summaries=summaries_data,rows=rows))
    macro_settings=dict(rigidityMin=rlo,rigidityMax=rhi,axisMin=0,axisMax=120,
        displayRebin=args.display_rebin,roiHalfWidth=args.roi_half_width,centerDisplayMax=args.center_display_max,
        widthDisplayMax=args.width_display_max,representativeRigidities=','.join(format(v,'.17g') for v in args.representative_rigidity),
        representativeTolerance=1.e-8,auditRelativeTolerance=1.e-12)
    for response in ('gaussian','lognormal'):
        for background in ('linear','exponential'):
            macro_settings[response+'Root_'+background]=str(studies[response,background]/'tag_fits.root')
    tsv(output/'settings.tsv',macro_settings.items())
    for name in ('comparison_settings.json','comparison.tsv','comparison.json','settings.tsv'):
        snapshots[str(output/name)]=digest(output/name)
    save(output/'context.json',dict(schemaVersion=1,sourceStudies={r+'_'+b:str(p) for (r,b),p in studies.items()},
        originalInputFingerprint=inputs[0],protected=protected,snapshots=snapshots,productionEnabled=False,
        noFit=True,createdJst=datetime.now(timezone(timedelta(hours=9))).isoformat()))
    def preserve():
        for path,expected in {**protected,**snapshots}.items():
            if digest(path)!=expected:
                raise ValueError('Original/source/snapshot changed: '+path)
        if fingerprint(Path(inputs[0]['path']))!=inputs[0]:
            raise ValueError('Original ROOT fingerprint changed.')
    preserve()
    if args.prepare_only:
        save(output/'prepared.json',dict(prepared=True,productionEnabled=False,noFit=True,rows=len(rows)))
        print('NUCLEAR_RESPONSE_COMPARISON_PREPARED',output,flush=True)
        return
    macro=snapshot/'CompareNuclearResponses.C'
    call=str(macro)+'++('+','.join(json.dumps(str(p)) for p in (output,output/'settings.tsv',output/'comparison.tsv'))+')'
    command=['bash',str(PROJECT/'script/singularity_root_batch.sh'),str(mainconf),'root4star','-l','-b','-q',call]
    save(output/'command.json',command);print('NUCLEAR_RESPONSE_COMPARISON_RUN',output,flush=True)
    try:
        with (output/'root.log').open('w') as log:
            result=subprocess.run(command,cwd=PROJECT,stdout=log,stderr=subprocess.STDOUT,timeout=args.timeout_seconds)
    finally:
        preserve()
    text=(output/'root.log').read_text(errors='replace')
    if result.returncode!=0 or 'NUCLEAR_RESPONSE_COMPARISON_COMPLETE' not in text:
        raise RuntimeError('Comparison failed; existing data preserved. See '+str(output/'root.log'))
    audit=read_json(output/'native_comparison_audit.json')
    pdfs={name:validate_pdf(output/name,audit[key]) for name,key in
          [('response_parameters.pdf','summaryPages'),('saved_spectra_overlays.pdf','overlayPages')]}
    save(output/'pdf_geometry.json',pdfs)
    # Render one page from each report for readback/visual inspection, never modify PDFs.
    gs=shutil.which('gs')
    if gs:
        for name in pdfs:
            png=output/(Path(name).stem+'_page1.png')
            render=[gs,'-dSAFER','-dBATCH','-dNOPAUSE','-sDEVICE=png16m','-r110','-dFirstPage=1','-dLastPage=1',
                    '-sOutputFile='+str(png),str(output/name)]
            subprocess.run(render,check=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=60)
    report=['# Gaussian／log-normalの保存fit比較','',
        '既存のfit結果を同じnative入力・fit窓で比較しました。再fit、PID変更、本番採用は行っていません。','',
        '- 中心はGaussianの平均とlog-normalのmodeです。幅はGaussianの標準偏差とlog-normalの局所曲率幅で、同じ定義ではありません。',
        '- 赤点は核種同定未検証、または数値失敗・zero yield・軸外・境界などの注意がある点です。未支持の点は線で接続していません。',
        '- dと4Heのdevianceは同じTAG_D fitの値です。独立データとして二重計数しません。',
        '- Δdevianceは非nestedな関数モデルの比較指標であり、有意度ではありません。','',
        '| Background | TAG sample | Native slices | Gaussian D | Log-normal D | Log/Gauss D | Log-normal lower D |',
        '|---|---|---:|---:|---:|---:|---:|']
    for row in summaries_data:
        report.append('| {background} | {sample} | {nativeSlices} | {gaussianD:.6g} | {lognormalD:.6g} | {lognormalOverGaussianDeviance:.4g} | {lognormalLowerDevianceSlices} |'.format(**row))
    report+=['','## 出力','',
        '- [中心・幅・fit-window devianceの比較](response_parameters.pdf)',
        '- [保存スペクトル：全dE/dx軸とpeak ROI、linear/log](saved_spectra_overlays.pdf)',
        '- [全native区間の数値・警告](comparison.tsv)',
        '- [同じnative count/error/flow・軸・fit maskの照合](native_comparison_audit.json)',
        '- comparison_display_inputs.rootには未変更のnative観測分布と表示用に加算した分布を別名で保存しています。',
        '', '数値収束・モデル比較の改善だけでは核種同定、TOF→ALLの移植、intrinsic幅、purityの検証を意味しません。']
    (output/'README.md').write_text('\n'.join(report)+'\n')
    preserve()
    products={str(p.relative_to(output)):digest(p) for p in output.iterdir() if p.is_file() and p.name!='context.json'}
    save(output/'completed.json',dict(executionComplete=True,noFit=True,productionEnabled=False,
        physicsValidationComplete=False,inputFingerprintPreserved=True,protectedSourcesPreserved=True,
        rows=len(rows),nativeFitPairs=len(model_keys[0])*2,pdfReadbackGeometryValid=True,products=products,
        completedJst=datetime.now(timezone(timedelta(hours=9))).isoformat()))
    print('NUCLEAR_RESPONSE_COMPARISON_COMPLETE',output,flush=True)


if __name__=='__main__':
    main()
