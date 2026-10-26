#!/bin/bash
# Presentation-only replay: no fitting, event processing, or production setting change.
set -euo pipefail
project=/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer
cd "$project"
study=rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015
views="$study/views_final_001"
mainconf=config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_allNucleiTof20261015.yaml
cache=rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_tof_20261014/inputs_001/cache.root
[[ $# == 1 ]] || { echo 'Expected overview|parameters|d|t|he3|he4'; exit 2; }
mkdir -p "$views/source" "$views/inputs"
snapshot() {
  if [[ -e "$2" ]]; then cmp "$1" "$2"; else cp -n "$1" "$2"; fi
}
for name in PlotNuclearCutoff.C PlotNuclearParameters.C PlotCutoffSpectra.C; do
  snapshot "tools/nuclear_pid_cutoff/$name" "$views/source/$name"
done
snapshot script/singularity_root_batch.sh "$views/source/singularity_root_batch.sh"
for name in summary.tsv native_parameters.tsv baseline_curves.tsv; do
  snapshot "$study/summary_free_sensitivity_001/$name" "$views/inputs/$name"
done
snapshot "$study/adopted_ranges_001/adopted_ranges.tsv" "$views/inputs/adopted_ranges.tsv"
case "$1" in
  overview)
    output="$views/purity_cutoff_overview.pdf"
    invocation="$views/source/PlotNuclearCutoff.C+(\"$cache\",\"$views/inputs/summary.tsv\",\"$views/inputs/baseline_curves.tsv\",\"$output\",0.4,2.0,0.5,\"$views/inputs/adopted_ranges.tsv\")"
    ;;
  parameters)
    output="$views/response_parameters.pdf"
    invocation="$views/source/PlotNuclearParameters.C+(\"$cache\",\"$views/inputs/native_parameters.tsv\",\"$views/inputs/baseline_curves.tsv\",\"$output\",0.4,2.0)"
    ;;
  d|t|he3|he4)
    output="$views/fit_spectra_$1.pdf"
    invocation="$views/source/PlotCutoffSpectra.C+(\"$views/models_$1.tsv\",\"$views/slices_$1.tsv\",\"$views/settings.tsv\",\"$output\")"
    ;;
  *) echo 'Unknown explicit figure name' >&2; exit 2;;
esac
[[ ! -e "$output" && ! -e "$views/$1.log" ]] || { echo 'Refusing to overwrite a prior rendering' >&2; exit 2; }
printf '%s\n' "bash script/singularity_root_batch.sh $mainconf root4star -l -b -q '$invocation'" > "$views/$1.command.txt"
bash script/singularity_root_batch.sh "$mainconf" root4star -l -b -q "$invocation" 2>&1 | tee "$views/$1.log"
[[ -s "$output" ]] || { echo 'Expected PDF not produced' >&2; exit 2; }
