#!/usr/bin/env bash
# Local saved-histogram study only: no farm submission, PID writes, or mail.
set -euo pipefail
if [[ $# -ne 2 || "$1" != "--output-dir" ]]; then
  echo "Usage: bash run_pilot_suite.sh --output-dir PREPARED_NEW_PILOT" >&2
  exit 2
fi
task_project=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd -P)
cd "$task_project"
task_output="$2"
[[ -f "$task_output/plan.json" && ! -e "$task_output/suite.log" ]] || {
  echo "Require a prepared pilot with no previous suite.log" >&2
  exit 2
}
set -o noclobber
exec >"$task_output/suite.log" 2>&1
printf 'CORETAIL_SUITE_STARTED %s\n' "$(date -Is)"
python3 -B tools/nuclear_pid_cutoff/coretail_20261019/run_coretail_pilot.py run \
  --output-dir "$task_output" --job pending
python3 -B tools/nuclear_pid_cutoff/coretail_20261019/run_coretail_pilot.py check \
  --output-dir "$task_output"
printf 'CORETAIL_SUITE_CHECKED %s\n' "$(date -Is)"
