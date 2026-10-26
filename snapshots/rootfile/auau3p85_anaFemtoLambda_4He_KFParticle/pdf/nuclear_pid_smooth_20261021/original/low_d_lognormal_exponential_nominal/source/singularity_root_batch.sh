#!/bin/bash
# Run ROOT utilities with the same ROOT5/STAR toolchain as FemtoLambda farm jobs.
# Usage: bash script/singularity_root_batch.sh MAINCONF hadd|root4star|root-config [ARGS...]
set -euo pipefail
[[ $# -ge 2 ]] || { echo "Usage: $0 MAINCONF hadd|root4star|root-config [ARGS...]" >&2; exit 2; }
SRB_MAINCONF="$1"
SRB_COMMAND="$2"
shift 2
case "$SRB_COMMAND" in hadd|root4star|root-config) ;; *) echo "ERROR: unsupported ROOT utility" >&2; exit 2;; esac
SRB_PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
case "$SRB_MAINCONF" in /*) ;; *) SRB_MAINCONF="$SRB_PROJECT/$SRB_MAINCONF";; esac
[[ -f "$SRB_MAINCONF" ]] || { echo "ERROR: missing explicit mainconf" >&2; exit 2; }
SRB_TAG="$(python3 "$SRB_PROJECT/script/analysis_info_helper.py" --library-tag --mainconf "$SRB_MAINCONF")"
[[ "$SRB_TAG" =~ ^[A-Za-z0-9_.-]+$ ]] || { echo "ERROR: invalid libraryTag" >&2; exit 2; }
SRB_STAR="/star/nfs4/AFS/star/packages/$SRB_TAG"
SRB_HOST_SYS=""
for SRB_CANDIDATE in sl73_x8664_gcc485 sl74_x8664_gcc485; do
  if [[ -d "$SRB_STAR/.$SRB_CANDIDATE" ]]; then SRB_HOST_SYS="$SRB_CANDIDATE"; break; fi
done
[[ -n "$SRB_HOST_SYS" ]] || { echo "ERROR: no batch-compatible STAR build" >&2; exit 2; }
shopt -s nullglob
SRB_ROOTS=(/cvmfs/star.sdcc.bnl.gov/star-spack/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/root-5.34.38-*/bin/root-config)
[[ ${#SRB_ROOTS[@]} -eq 1 ]] || { echo "ERROR: cannot uniquely resolve batch ROOT5" >&2; exit 2; }
SRB_ROOT_SYS="${SRB_ROOTS[0]%/bin/root-config}"
printf '[mainconf] %s\n[runtime] STAR=%s STAR_HOST_SYS=%s ROOTSYS=%s\n' "$SRB_MAINCONF" "$SRB_STAR" "$SRB_HOST_SYS" "$SRB_ROOT_SYS"
exec singularity exec --pwd "$SRB_PROJECT" \
  -B /gpfs:/gpfs -B /cvmfs:/cvmfs -B /star/nfs4/AFS:/star/nfs4/AFS -B /tmp:/tmp \
  --env "ROOTSYS=$SRB_ROOT_SYS" --env "STAR=$SRB_STAR" --env "STAR_HOST_SYS=$SRB_HOST_SYS" \
  --env "LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}" \
  /cvmfs/singularity.opensciencegrid.org/star-bnl/star-sw:latest /bin/bash -lc '
    set -euo pipefail
    project="$1"; shift
    export STAR_LIB="$STAR/.$STAR_HOST_SYS/lib" STAR_BIN="$STAR/.$STAR_HOST_SYS/bin"
    export PATH="$ROOTSYS/bin:$STAR_BIN:$STAR/mgr:$PATH"
    export LD_LIBRARY_PATH="$project/lib:$STAR_LIB:$ROOTSYS/lib:${LD_LIBRARY_PATH:-}"
    cd "$project"
    exec "$@"
  ' star-root-batch "$SRB_PROJECT" "$SRB_COMMAND" "$@"
