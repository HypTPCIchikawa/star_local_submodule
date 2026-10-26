#!/usr/bin/env python3
"""Freeze and export native diagnostic counts; no fit/PID or input mutation."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess

from prepare_tag_study import PROJECT, owned, quote, read, require, save, sha
from render_tag_study import verify_inputs


def run(fit_dir, output):
    plan = read(fit_dir/'plan.json')
    verified = verify_inputs(plan, fit_dir)
    require(not output.exists() and output.parent == fit_dir.parent,
            'Require a fresh direct study sibling output')
    cfg = plan['config']
    wrapper, main = owned(cfg['wrapper']), owned(cfg['mainconf'])
    output.mkdir()
    protected = dict(verified)
    sources = ('tools/nuclear_pid_restart_20261020/ExportRestartSpectra.C',
               'tools/nuclear_pid_restart_20261020/run_export_restart_spectra.C',
               'tools/nuclear_pid_restart_20261020/FitRestartSlices.C',
               'tools/nuclear_pid_tof/NuclearFitCommon.h')
    for name in sources:
        source, target = owned(name), output/'source'/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(source), str(target))
        protected[str(target)] = sha(source)
        require(sha(target) == protected[str(target)], 'Source freeze mismatch')
    for name in ('settings.tsv', 'slices.tsv'):
        protected[str(fit_dir/name)] = sha(fit_dir/name)
    protected[str(wrapper)] = sha(wrapper)
    protected[str(main)] = sha(main)
    runner = output/'source'/sources[1]
    macro = output/'source'/sources[0]
    table = output/'native_counts.tsv'
    invocation = str(runner)+'('+','.join(quote(p) for p in
        (macro, output/'build', plan['cache'], fit_dir/'settings.tsv',
         fit_dir/'slices.tsv', table))+')'
    command = ['bash', str(wrapper), str(main), 'root4star', '-l', '-b', '-q', invocation]
    save(output/'provenance.json', dict(sourceFit=str(fit_dir), command=command,
         protected=protected, productionEnabled=False, noFit=True,
         storageMappingHypothesis='ROOT default12 Float16; actual on-file precision not yet verified'))
    with (output/'export.log').open('x') as stream:
        subprocess.run(command, cwd=PROJECT, stdout=stream, stderr=subprocess.STDOUT,
                       timeout=300, check=True)
    require('RESTART_SPECTRA_EXPORT_COMPLETE' in (output/'export.log').read_text(),
            'Missing export/readback completion marker')
    for name, expected in protected.items():
        require(sha(Path(name)) == expected, 'Protected input/source changed: '+name)
    with table.open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    # Do not infer that all source values were integer fills from this representative sample.
    summary = dict(rowCount=len(rows), sliceCount=len(set(row['id'] for row in rows)),
        representativeMaximumCount=max(float(row['count']) for row in rows),
        nonintegerCountBins=sum(float(row['count']) != int(float(row['count'])) for row in rows),
        countErrorConsistent=all(abs(float(row['error'])**2-float(row['count'])) <=
                                1e-8*max(1., float(row['count'])) for row in rows))
    save(output/'completed.json', dict(status='completed', nativeReadbackExact=True,
         inputHashesUnchanged=True, noFit=True, productionEnabled=False,
         tableSha256=sha(table), summary=summary))
    print(json.dumps(summary, sort_keys=True), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fit-dir', required=True)
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    run(owned(args.fit_dir), owned(args.output_dir, False))
