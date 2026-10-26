#!/bin/bash
# One user-authorized recovery only; run under tmux with explicit /bin/bash.
set -euo pipefail
recovery_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
if [[ -e "$recovery_dir/PAUSED_BY_USER.md" ]]; then
  echo "Stopped at user request. Inspect saved state and obtain a new resume instruction."
  exit 1
fi
project_dir="$(cd "$recovery_dir/../../.." && pwd -P)"
retry_ana=auau3p85_anaFemtoLambda_4He_KFParticle_recovery_20261003
started_dir="$recovery_dir/workflow.started"
mkdir "$started_dir"  # Permanent atomic guard; never remove to resubmit.
TZ=Asia/Tokyo date -Is > "$started_dir/started_at.txt"
hostname > "$started_dir/host.txt"
trap 'rc=$?; if (( rc != 0 )); then echo "Recovery stopped (exit=$rc). Inspect evidence; do not automatically resubmit."; fi' EXIT
cd "$project_dir"
python3 -B "$recovery_dir/check_original_runtime.py" > "$recovery_dir/runtime_match.json"
# This is an explicitly new attempt after the preserved interpreted-audit pilot.
# ROOT5's error-handler API needs the precompiled ACLiC macro.
[[ -s "$recovery_dir/checkRecoveryPico_C.so" ]]
python3 -B -u "$recovery_dir/prepare_inputs.py" --workers 3 --audit-suffix + --attempt-offset 1
python3 -B "$recovery_dir/verify_ready_groups.py"
[[ ! -e "$recovery_dir/submit.log" ]]
mkdir "$recovery_dir/submission.started"
cd "$project_dir/job/run"
./submit.sh "../joblist/joblist_${retry_ana}.xml" > "$recovery_dir/submit.log" 2>&1
grep -q '^Scheduling successful' "$recovery_dir/submit.log"
report_path=$(sed -n 's/^Wrote scheduling report to : *//p' "$recovery_dir/submit.log" | tail -n 1 | tr -d '\r')
report_name=$(basename "$report_path" .report)
job_id="${report_name#${retry_ana}}"
[[ "$job_id" =~ ^[0-9A-Fa-f]{32}$ ]]
runmeta="$project_dir/job/run/runmeta/runmeta_${retry_ana}_${job_id}.json"
[[ -s "$runmeta" ]]
echo "$job_id" > "$recovery_dir/submission.started/jobid.txt"
cd "$project_dir"
python3 -B "$recovery_dir/verify_ready_groups.py" --submitted-jobid "$job_id"
state_dir="$project_dir/job/run/merge_notify/${retry_ana}_${job_id}"
mkdir -p "$state_dir"
echo "Recovery submitted: SUMS $job_id. Starting strict retry-only merge/archive/mail."
# Retry output is isolated. Never merge the parent partial result twice.
exec nice -n 10 python3 -B -u script/watch_job_merge_notify.py \
  --runmeta "$runmeta" --event-log "/tmp/yichikawa/sched${job_id}.condor.log" \
  --root-runner script/singularity_root_batch.sh \
  --email yudai.ichikawa.d3@tohoku.ac.jp --watch > "$state_dir/monitor.log" 2>&1
