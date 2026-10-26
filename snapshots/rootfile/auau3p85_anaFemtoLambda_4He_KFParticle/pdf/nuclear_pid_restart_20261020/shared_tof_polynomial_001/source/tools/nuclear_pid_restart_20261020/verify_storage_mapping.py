#!/usr/bin/env python3
"""Run the Float16 mathematical oracle in a fresh ROOT5 build; no data fit."""
import argparse
from pathlib import Path
import shutil
import subprocess

from prepare_tag_study import PROJECT, owned, quote, require, save, sha


def run(output):
    require(not output.exists(), 'Require a fresh diagnostic output directory')
    wrapper = owned('script/singularity_root_batch.sh')
    main = owned('config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml')
    sources = ['tools/nuclear_pid_restart_20261020/'+name for name in
               ('Float16BinMapping.h', 'TestFloat16BinMapping.C', 'run_test_float16_mapping.C')]
    output.mkdir()
    protected = {str(path): sha(path) for path in (wrapper, main)}
    for name in sources:
        source, frozen = owned(name), output/'source'/name
        frozen.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(source), str(frozen))
        protected[str(frozen)] = sha(source)
        require(sha(frozen) == protected[str(frozen)], 'Source freeze mismatch')
    runner, macro = output/'source'/sources[2], output/'source'/sources[1]
    invocation = str(runner)+'('+','.join(quote(p) for p in (macro, output/'build'))+')'
    command = ['bash', str(wrapper), str(main), 'root4star', '-l', '-b', '-q', invocation]
    save(output/'provenance.json', dict(command=command, protected=protected,
         productionEnabled=False, actualOnFileEncodingAsserted=False))
    with (output/'test.log').open('x') as stream:
        subprocess.run(command, cwd=PROJECT, stdout=stream, stderr=subprocess.STDOUT,
                       timeout=300, check=True)
    require('FLOAT16_BIN_MAPPING_TEST_PASS' in (output/'test.log').read_text(), 'Missing test marker')
    for name, expected in protected.items():
        require(sha(Path(name)) == expected, 'Protected source changed: '+name)
    save(output/'completed.json', dict(status='passed', protectedInputsUnchanged=True,
         dataFit=False, productionEnabled=False))
    print('FLOAT16_ORACLE_VERIFIED', output, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', required=True)
    run(owned(parser.parse_args().output_dir, False))
