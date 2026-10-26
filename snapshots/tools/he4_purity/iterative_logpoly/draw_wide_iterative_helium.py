#!/usr/bin/env python3
"""New DISPLAY-ONLY full-native-axis views of an immutable completed fit phase."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import yaml

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def fingerprint(path):
    s = path.stat()
    return dict(path=str(path), size=s.st_size, mtime_ns=s.st_mtime_ns, device=s.st_dev, inode=s.st_ino)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('study', 'phase-dir', 'mainconf', 'settings', 'output-dir', 'phase-label'):
        parser.add_argument('--'+name, required=True, type=str if name=='phase-label' else Path)
    args = parser.parse_args()
    study, phase = args.study.resolve(strict=True), args.phase_dir.resolve(strict=True)
    if phase != study and study not in phase.parents:
        raise ValueError('Phase must belong to the explicitly completed study.')
    context_path, completion_path = study/'context.json', study/'completed.json'
    context, completion = json.loads(context_path.read_text()), json.loads(completion_path.read_text())
    if context['productionEnabled'] or completion['productionEnabled'] or not completion['executionComplete']:
        raise ValueError('Only a completed study-only phase can be viewed.')
    mainconf = args.mainconf.resolve(strict=True)
    cfg = context['settings']
    if mainconf != (PROJECT/cfg['sample']['mainconf']).resolve(strict=True):
        raise ValueError('Explicit runtime mainconf must match the saved study.')
    runtime = yaml.safe_load(mainconf.read_text())
    analysis_path = (PROJECT/'config'/runtime['analysis']).resolve(strict=True)
    analysis = yaml.safe_load(analysis_path.read_text())
    tags = analysis['starTag']
    if (tags['collision'], tags['year'], tags['productionTag'], tags['libraryTag']) != ('auau3p85','year2021','P24iy','SL24y'):
        raise ValueError('ROOT5/SL24y 3p85/2021P24iy runtime is required.')
    source_fingerprint = context.get('source', context.get('input'))
    if source_fingerprint is None or fingerprint(Path(source_fingerprint['path'])) != source_fingerprint:
        raise ValueError('Original merged input fingerprint differs.')
    protected = {str(context_path):digest(context_path), str(completion_path):digest(completion_path),
                 str(mainconf):digest(mainconf), str(analysis_path):digest(analysis_path)}
    for encoded, expected in context['protected'].items():
        p = Path(encoded)
        p = p if p.is_absolute() else PROJECT/p
        protected[str(p)] = expected
    for encoded, expected in context['snapshots'].items():
        protected[encoded] = expected
    for relative, expected in completion['products'].items():
        p = (study/relative).resolve(strict=True)
        if study not in p.parents:
            raise ValueError('Completed product escapes study.')
        protected[str(p)] = expected
    def preserve():
        if fingerprint(Path(source_fingerprint['path'])) != source_fingerprint:
            raise ValueError('Original input changed.')
        for encoded, expected in protected.items():
            if digest(Path(encoded)) != expected:
                raise ValueError('Frozen source/product changed: '+encoded)
    preserve()
    settings_path = args.settings.resolve(strict=True)
    settings = yaml.safe_load(settings_path.read_text())
    if settings['schemaVersion'] != 1:
        raise ValueError('Unsupported wide-display settings schema.')
    manifest = phase/'views/slices.tsv'
    with manifest.open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    if [int(row['index']) for row in rows] != list(range(len(rows))):
        raise ValueError('Expected complete sorted prefix slices.')
    display = dict(cfg['display'], **settings['display'], expectedSlices=len(rows), allowSparseSlices=0)
    output = args.output_dir.resolve()
    allowed = (PROJECT/cfg['sample']['outputBase']).resolve()
    if output.exists() or allowed not in output.parents or study in output.parents:
        raise ValueError('Choose a new sibling display-output directory.')
    output.mkdir(parents=True)
    snapshots = output/'source'
    snapshots.mkdir()
    for path in (Path(__file__).resolve(), HERE/'WideIterativeHelium.C', settings_path):
        shutil.copy2(path, snapshots/path.name)
    frozen_renderer = study/'source/tools/he4_purity/iterative_logpoly/RenderIterativeHelium.C'
    shutil.copy2(frozen_renderer, snapshots/'FrozenRenderIterativeHelium.C')
    shutil.copy2(manifest, output/'slices.tsv')
    (output/'settings.tsv').write_text(''.join('%s\t%s\n'%(key,value) for key,value in display.items()))
    (output/'effective_display_settings.yaml').write_text(yaml.safe_dump(display,sort_keys=False))
    provenance = dict(schemaVersion=1, mode='DISPLAY ONLY', noRefit=True, productionChanged=False,
        inputFingerprint=source_fingerprint, sourceStudy=str(study), sourcePhase=str(phase),
        protectedPreRun=protected, frozenRendererSha256=digest(frozen_renderer),
        display=display, phaseLabel=args.phase_label,
        warning='Unadopted model-dependent isotope labels. Curves restricted to saved fit mask; outside is data only.')
    (output/'context.json').write_text(json.dumps(provenance,indent=2)+'\n')
    call = str(snapshots/'WideIterativeHelium.C')+'++('+','.join(json.dumps(str(p)) for p in
        (phase/'phase.root', output/'slices.tsv', output/'settings.tsv', output, args.phase_label))+')'
    command = ['bash',str(PROJECT/'script/singularity_root_batch.sh'),str(mainconf),'root4star','-l','-b','-q',call]
    (output/'command.json').write_text(json.dumps(command,indent=2)+'\n')
    with (output/'root.log').open('w') as log:
        result = subprocess.run(command,cwd=PROJECT,stdout=log,stderr=subprocess.STDOUT,timeout=settings['executionTimeoutSeconds'])
    if result.returncode or 'WIDE_ITERATIVE_HELIUM_COMPLETE' not in (output/'root.log').read_text():
        raise RuntimeError('Display-only drawing failed; see saved root.log.')
    preserve()
    audit = json.loads((output/'render_audit.json').read_text())
    if not all(row['visiblePhysicalRange']==[0,120] and row['visibleCounts']==row['allNormalCount'] for row in audit['slicesAudit']):
        raise ValueError('Not all normal physical bins are visible.')
    products = {str(p.relative_to(output)):digest(p) for p in output.rglob('*') if p.is_file() and 'source' not in p.relative_to(output).parts}
    source_hashes = {str(p.relative_to(snapshots)):digest(p) for p in snapshots.rglob('*') if p.is_file()}
    (output/'completed.json').write_text(json.dumps(dict(executionComplete=True,noRefit=True,productionChanged=False,
        protectedPostRunVerified=True,inputFingerprintPreserved=True,products=products,sourceSha256=source_hashes,
        sameNativeVisibleBins=True,curvesOnlyInSavedFitMask=True),indent=2)+'\n')
    print('WIDE_ITERATIVE_CONTROL_COMPLETE',output,flush=True)

if __name__ == '__main__':
    main()
