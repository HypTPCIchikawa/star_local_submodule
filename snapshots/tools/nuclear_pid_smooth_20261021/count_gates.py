#!/usr/bin/env python3
import argparse,json,shutil,subprocess
from study import P,OLD,OUT,sha,save,read
parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['candidate','original']);args=parser.parse_args()
src=OLD/'independent_gate_counts_001';done=read(src/'completed.json');dst=OUT/(args.mode+'_gate_counts');dst.mkdir()
for name in ['ExtractNuclearGateCounts.C','FitNuclearCutoff.C','NuclearFitCommon.h','NuclearLognormalFitCommon.h']:
 if sha(src/name)!=done['products'][name]:raise ValueError('Historical counter changed')
 shutil.copy2(src/name,dst/name)
run=OUT/args.mode/'low_he4_gaussian_linear_nominal'
for n in ['shapes.tsv','gates.tsv','settings.tsv']:shutil.copy2(run/n,dst/n)
cmd=['bash',str(P/'script/singularity_root_batch.sh'),str(P/'config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml'),'root4star','-l','-b','-q',str(dst/'ExtractNuclearGateCounts.C')+'+('+','.join(json.dumps(str(p)) for p in [done['cachePath'],dst/'shapes.tsv',dst/'gates.tsv',dst/'settings.tsv',dst/'gate_counts.tsv'])+')']
save(dst/'command.json',cmd)
with (dst/'run.log').open('w') as f:status=subprocess.run(cmd,cwd=P,stdout=f,stderr=subprocess.STDOUT).returncode
if status or 'NUCLEAR_GATE_COUNTS_COMPLETE rows=640' not in (dst/'run.log').read_text():raise RuntimeError('Independent count failed')
save(dst/'completed.json',dict(executionComplete=True,productionEnabled=False,fitIndependent=True,rows=640,products={p.name:sha(p) for p in dst.iterdir() if p.is_file() and p.suffix not in ['.so','.d']}))
print('GATE_COUNTS_COMPLETE',args.mode)
