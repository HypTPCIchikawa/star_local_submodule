#!/usr/bin/env python3
"""Run synthetic response tests against a frozen study source/settings snapshot."""
import argparse
from pathlib import Path
import subprocess

from prepare_tag_study import PROJECT, owned, quote, read, require, save, sha, unchanged


def run(fit_dir, output):
    plan = read(fit_dir/'plan.json')
    unchanged(plan)
    require(Path(plan['config']['fitSource']).name == 'FitSplitNormalSlices.C', 'Wrong response study')
    require(not output.exists() and output.parent == fit_dir, 'Require a fresh test child directory')
    sources = fit_dir/'source/tools/nuclear_pid_restart_20261020'
    runner, macro = sources/'run_test_split_normal.C', sources/'TestSplitNormalFit.C'
    require(runner.is_file() and macro.is_file(), 'Missing frozen test source')
    cfg = plan['config']
    output.mkdir()
    report = output/'test_report.json'
    invocation = str(runner)+'('+','.join(quote(p) for p in
        (macro, output/'build', fit_dir/'settings.tsv', report))+')'
    command = ['bash', str(owned(cfg['wrapper'])), str(owned(cfg['mainconf'])),
               'root4star', '-l', '-b', '-q', invocation]
    save(output/'command.json', command)
    with (output/'test.log').open('x') as stream:
        subprocess.run(command, cwd=PROJECT, stdout=stream, stderr=subprocess.STDOUT,
                       timeout=300, check=True)
    require('SPLIT_NORMAL_TEST_PASS' in (output/'test.log').read_text(), 'Missing test marker')
    result = read(report)
    require(result.get('status') == 'passed' and result.get('productionEnabled') is False and
            result.get('syntheticOnly') is True and result.get('coefficientsUpdated') is False,
            'Missing isolated successful synthetic test report')
    unchanged(plan)
    save(output/'completed.json', dict(status='passed', protectedInputsUnchanged=True,
         reportSha256=sha(report), productionEnabled=False))
    print('SPLIT_NORMAL_ORACLE_VERIFIED', output, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fit-dir', required=True)
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    run(owned(args.fit_dir), owned(args.output_dir, False))
