#!/usr/bin/env python3
"""Read-only PDF/input validation and generated SHA provenance (no fit)."""
import csv
import datetime
import hashlib
import json
import pathlib
import re
import subprocess

PROJECT = pathlib.Path('/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer')
BASE = PROJECT / 'rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015'
HERE = BASE / 'views_final_001'
EXPECTED = {'purity_cutoff_overview.pdf': 3, 'response_parameters.pdf': 5,
            'fit_spectra_d.pdf': 12, 'fit_spectra_t.pdf': 12,
            'fit_spectra_he3.pdf': 12, 'fit_spectra_he4.pdf': 16}


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    output = HERE / 'completed.json'
    if output.exists():
        raise RuntimeError('Refusing to overwrite a prior validation manifest')
    pages = {}
    for name, expected in EXPECTED.items():
        result = subprocess.run(['gs', '-dBATCH', '-dNOPAUSE', '-sDEVICE=nullpage',
                                 str(HERE / name)], stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, universal_newlines=True,
                                check=True)
        (HERE / (name + '.readback.log')).write_text(result.stdout)
        match = re.search(r'Processing pages 1 through ([0-9]+)\.', result.stdout)
        if not match or int(match.group(1)) != expected:
            raise RuntimeError('Wrong page count: ' + name)
        if 'Error:' in result.stdout or 'Unrecoverable' in result.stdout:
            raise RuntimeError('PDF parser failure: ' + name)
        pages[name] = expected

    model_inputs = {}
    for species in ('d', 't', 'he3', 'he4'):
        with (HERE / ('models_' + species + '.tsv')).open() as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        if len(rows) != 4:
            raise RuntimeError('Exactly four nominal models are required')
        for row in rows:
            path = PROJECT / row['path']
            model_inputs[row['path']] = {'size': path.stat().st_size, 'sha256': sha(path)}

    cache = PROJECT / 'rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_tof_20261014/inputs_001/cache.root'
    inputs = {str(cache.relative_to(PROJECT)): {'size': cache.stat().st_size, 'sha256': sha(cache)}}
    for directory in ('inputs', 'source'):
        for path in sorted((HERE / directory).glob('*')):
            if path.is_file() and path.suffix in ('.C', '.sh', '.tsv'):
                inputs[str(path.relative_to(HERE))] = {'size': path.stat().st_size, 'sha256': sha(path)}
    for relative in ('summary_free_sensitivity_001/completed.json',
                     'summary_free_sensitivity_001/summary.json',
                     'adopted_ranges_001/completed.json',
                     'adopted_ranges_001/adopted_ranges.json'):
        path = BASE / relative
        inputs[relative] = {'size': path.stat().st_size, 'sha256': sha(path)}
    artifacts = {}
    for path in sorted(HERE.iterdir()):
        if path.is_file() and path.name != 'completed.json':
            artifacts[path.name] = {'size': path.stat().st_size, 'sha256': sha(path)}
    manifest = {
        'status': 'completed', 'presentationOnly': True, 'refitted': False,
        'productionEnabled': False,
        'createdJST': datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat(),
        'runtime': 'STAR SL24y / sl73_x8664_gcc485 / ROOT 5.34.38',
        'pageCounts': pages, 'inputs': inputs, 'savedNominalModelRoots': model_inputs,
        'artifacts': artifacts,
        'visualQA': {
            'overview_page_1.png': 'PDF page 1: model envelopes, unassessed points and adopted/veto R ranges visible',
            'overview_page_2.png': 'PDF page 2: full dEdx ALL COLZ and unchanged center/window guides visible',
            'overview_page_3.png': 'PDF page 3: independent raw counts and adopted counts visible',
            'parameters_page_1.png': 'Response PDF page 4: He4 center/width and nominal models visible',
            'parameters_page_2.png': 'Response PDF page 5: ALL COLZ with supported fitted centers visible',
            'he4_full_page_1.png': 'He4 PDF page 1: native full range linear/log, no normalization',
            'he4_fit_page_1.png': 'He4 PDF page 3: native fit-window linear/log, no normalization'
        },
        'limits': [
            'Purity and response values are copied from the completed 20261015 study, not refitted.',
            'Model envelope is not a complete systematic uncertainty or a confidence interval.',
            'Unsupported numerical fits remain diagnostic; drawing a component does not validate it.',
            'Candidate retention concerns inclusive TPC tracks at R=0.4-2, not signal efficiency or CF pairs.',
            'Final PID coefficients and existing production mainconf are unchanged by rendering.'
        ]
    }
    required = set(manifest['visualQA'])
    if not required.issubset(artifacts):
        raise RuntimeError('Visual QA images are missing: ' + ', '.join(sorted(required - set(artifacts))))
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    print('FINAL_VIEWS_COMPLETE', json.dumps(pages, sort_keys=True))


if __name__ == '__main__':
    main()
