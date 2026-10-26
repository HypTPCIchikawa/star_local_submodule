#!/bin/bash
# One-shot, user-authorized production. Run under tmux with explicit bash.
# A surviving *.started directory forbids submission again after interruption.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
run_dir="$PWD"
project_dir="$(cd ../.. && pwd)"
ana_name=auau3p85_anaFemtoLambda_4He_KFParticle
launch_name=production_3p85_2021P24iy_20261001
submit_log="$run_dir/${launch_name}_submit.log"
started_dir="$run_dir/${launch_name}.started"
trap 'rc=$?; if (( rc != 0 )); then echo "Workflow stopped (exit=$rc); inspect logs; NEVER resubmit automatically."; fi' EXIT

[[ ! -e "$submit_log" ]] || { echo "Refusing to overwrite submission log"; exit 1; }
mkdir "$started_dir"  # atomic, permanent duplicate-submission guard
TZ=Asia/Tokyo date -Is > "$started_dir/started_at.txt"
hostname > "$started_dir/host.txt"
echo "Starting one-shot 2021 P24iy production; submission log: $submit_log"
./submit.sh "../joblist/joblist_${ana_name}.xml" > "$submit_log" 2>&1
grep -q '^Scheduling successful' "$submit_log"
report_path=$(sed -n 's/^Wrote scheduling report to : *//p' "$submit_log" | tail -n 1 | tr -d '\r')
report_name=$(basename "$report_path" .report)
job_id="${report_name#${ana_name}}"
[[ "$job_id" =~ ^[0-9A-Fa-f]{32}$ ]] || { echo "Cannot identify this submission"; exit 1; }
runmeta="$run_dir/runmeta/runmeta_${ana_name}_${job_id}.json"
[[ -s "$runmeta" ]] || { echo "Missing runmeta: $runmeta"; exit 1; }

# Audit the actual assigned dataset before handing over to the strict watcher.
python3 -B - "$runmeta" "$ana_name" "$job_id" <<'PY'
import collections
import json
from pathlib import Path
import re
import sys
import xml.etree.ElementTree as ET

path, ana, jobid = sys.argv[1:]
meta = json.loads(Path(path).read_text())
if meta['anaName'] != ana or meta['jobid'] != jobid:
    raise SystemExit('Runmeta does not match this submission')
xml = ET.parse(meta['submit']['joblistSnapshot']['path']).getroot()
expected_url = ('catalog:star.bnl.gov?trgsetupname=production_3p85GeV_fixedTarget_2021,'
                'production=P24iy,filetype=daq_reco_picoDst,filename~st_physics,'
                'filename!~adc,storage=local')
inputs = xml.findall('input')
if len(inputs) != 1 or inputs[0].get('URL') != expected_url or inputs[0].get('nFiles') != 'all':
    raise SystemExit('Submitted catalog scope is not the authorized 2021 P24iy dataset')
run = Path(meta['projectRoot']) / 'job/run'
pattern = re.compile(re.escape(meta['jobPrefix']) + r'_(\d+)\.list$')
lists = {int(m.group(1)): p for p in run.iterdir() if (m := pattern.fullmatch(p.name))}
if not lists or sorted(lists) != list(range(len(lists))):
    raise SystemExit('Missing/non-contiguous subjob lists')
files, split_counts = [], collections.Counter()
catalog_events = 0
for index in sorted(lists):
    rows = [line.split() for line in lists[index].read_text().splitlines() if line.strip()]
    if not 1 <= len(rows) <= 5:
        raise SystemExit('Invalid file count in {}'.format(lists[index]))
    split_counts[len(rows)] += 1
    for row in rows:
        uri = row[0]
        if ('/production_3p85GeV_fixedTarget_2021/' not in uri or '/P24iy/2021/' not in uri
                or not re.fullmatch(r'st_physics_22\d{6}_[^/]+\.picoDst\.root', uri.rsplit('/', 1)[-1])
                or 'adc' in uri.rsplit('/', 1)[-1]):
            raise SystemExit('Unexpected input dataset: ' + uri)
        files.append(uri)
        if len(row) != 2 or not row[1].isdigit():
            raise SystemExit('Cannot audit catalog event count: ' + ' '.join(row))
        catalog_events += int(row[1])
if len(files) != len(set(files)) or len(files) != len({p.rsplit('/', 1)[-1] for p in files}):
    raise SystemExit('Duplicate production inputs')
print(json.dumps(dict(jobid=jobid, subjobs=len(lists), uniqueInputFiles=len(files),
                     filesPerSubjob=dict(split_counts), catalogEvents=catalog_events,
                     measuredPicoDstEntries=None, scope='2021 P24iy non-ADC local'), indent=2))
PY

state_dir="$run_dir/merge_notify/${ana_name}_${job_id}"
mkdir -p "$state_dir"
echo "$job_id" > "$started_dir/jobid.txt"
echo "Submission and input audit passed. Starting strict merge/archive/email watcher."
echo "Watcher log: $state_dir/monitor.log"
cd "$project_dir"
# No --allow-partial: missing/failed outputs must result in a warning, not a partial merge.
exec nice -n 10 python3 -B -u script/watch_job_merge_notify.py \
  --runmeta "$runmeta" \
  --event-log "/tmp/yichikawa/sched${job_id}.condor.log" \
  --root-runner script/singularity_root_batch.sh \
  --email yudai.ichikawa.d3@tohoku.ac.jp --watch > "$state_dir/monitor.log" 2>&1
