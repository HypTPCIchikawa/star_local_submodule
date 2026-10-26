#!/usr/bin/env python3
"""Validate saved report/figure/source hashes and unchanged ROOT fingerprint."""
import argparse
import hashlib
import json
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    out = args.output_dir.resolve(strict=True)
    state = json.loads((out / 'run_state.json').read_text())
    manifest = json.loads((out / 'macro_manifest.json').read_text())
    report = json.loads((out / 'report_metadata.json').read_text())
    if state['phase'] != 'complete':
        raise RuntimeError('Run is not complete: ' + state['phase'])
    original = state['provenance']['input']
    st = Path(original['path']).stat()
    current = dict(path=original['path'], size=st.st_size, mtime_ns=st.st_mtime_ns,
                   device=st.st_dev, inode=st.st_ino)
    if current != original:
        raise RuntimeError('Input ROOT fingerprint changed.')
    for path, expected in state['provenance']['configs'].items():
        if sha(path) != expected:
            raise RuntimeError('Configuration changed: ' + path)
    for macro in manifest['macros']:
        for field, key in [('source', 'source_sha256'), ('generated', 'generated_sha256')]:
            if sha(macro[field]) != macro[key]:
                raise RuntimeError('Source hash mismatch: ' + macro[field])
        completed = state['stages'][Path(macro['generated']).stem]
        if completed['phase'] != 'complete' or completed['signature'].split(':')[0] != macro['generated_sha256']:
            raise RuntimeError('Macro stage does not match generated source: ' + macro['name'])
    if sha(out / 'macros/KFSummarySupport.h') != manifest['support']['sha256']:
        raise RuntimeError('Generated helper hash mismatch.')
    if report['phase'] != 'complete' or report['inputRoot'] != original['path']:
        raise RuntimeError('Report input/state does not match the completed run.')
    checks = {out / 'inspect.json': report['inspect_sha256'],
              out / 'macro_manifest.json': report['manifest_sha256'],
              out / report['report']['name']: report['report']['sha256'],
              out / Path(report['report']['name']).with_suffix('.tex'): report['tex_sha256'],
              Path(__file__).with_name('build_report.py'): report['generator_sha256']}
    for name, record in report['figures'].items():
        checks[out / name] = record['sha256']
    for record in report['fit_audits']:
        checks[Path(record['path'])] = record['sha256']
    for path, expected in checks.items():
        if sha(path) != expected:
            raise RuntimeError('Report dependency hash mismatch: ' + str(path))
    print(json.dumps(dict(phase='validated', input_unchanged=True,
                         macros=len(manifest['macros']), figures=len(report['figures']),
                         report_pages=report['report']['pages'], output_dir=str(out)), indent=2))


if __name__ == '__main__':
    main()
