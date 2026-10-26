#!/usr/bin/env python3
"""Run immutable study job manifests through the existing ROOT5/StChain path.

Outputs/logs must not exist. No source build, setup.sh, farm submit or default
configuration modification. A requested event count silently clamped by the
analysis macro is explicitly rejected here.
"""
import argparse
import concurrent.futures
import json
import os
from pathlib import Path
import re
import subprocess

REPO = Path(__file__).resolve().parents[2]
WRAPPER = REPO / "rootfile/auau13p5_anaLambda_KFParticle/purity_study_20260927/provenance/helpers/run_root5_readonly.sh"


def quote(value):
    return json.dumps(str(value))


def execute(job, no_export=False):
    output, log = Path(job["output"]), Path(job["log"])
    if output.exists() or log.exists():
        raise RuntimeError("Refusing existing output/log for " + job["id"])
    log.parent.mkdir(parents=True, exist_ok=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    mounts = sorted(set(re.findall(r"/star/data\d+", Path(job["input"]).read_text())))
    if not mounts:
        raise RuntimeError("No explicit NFS input mounts found")
    extra = ",".join(path + ":" + path for path in mounts)
    env["SINGULARITY_BIND"] = ",".join(x for x in [env.get("SINGULARITY_BIND", ""), extra] if x)
    call = "analysis/run_anaLambda_KFParticle.C({},{},{},{},{})".format(
        quote(job["input"]), quote(output), quote(job["id"]), job["events"], quote(job["mainconf"]))
    print("START {} {} events={}".format(job["energy"], job["id"], job["events"]), flush=True)
    with log.open("x") as stream:
        result = subprocess.run(["bash", str(WRAPPER), call], cwd=REPO, env=env,
                                stdout=stream, stderr=subprocess.STDOUT, timeout=7200)
    text = log.read_text(errors="replace")
    matches = re.findall(r"KF job=\S+ requested=(\d+) completed=(\d+) reconstructed=(\d+) status=(\d+) RealTime=([\d.e+-]+)", text)
    if result.returncode or len(matches) != 1:
        raise RuntimeError("Reconstruction failed: {} exit={} log={}".format(job["id"], result.returncode, log))
    requested, completed, reco, status = map(int, matches[0][:4])
    if requested != job["events"] or completed != job["events"] or status or reco <= 0:
        raise RuntimeError("Invalid input/completion counts: " + str(matches[0]))
    summary = dict(job, completed=completed, reconstructed=reco, seconds=float(matches[0][4]), applicationStatus=status)
    if not no_export:
        export = Path(job["export"])
        if export.exists():
            raise RuntimeError("Refusing existing export " + str(export))
        export.parent.mkdir(parents=True, exist_ok=True)
        export_log = log.with_name(log.stem + "_export.log")
        call = "common/macro/exportLambdaKfUpstreamStudy.C({},{},{})".format(quote(output), quote(export), job["events"])
        with export_log.open("x") as stream:
            checked = subprocess.run(["bash", str(WRAPPER), call], cwd=REPO, env=env,
                                     stdout=stream, stderr=subprocess.STDOUT, timeout=1800)
        if checked.returncode or "PASS purity export:" not in export_log.read_text(errors="replace"):
            raise RuntimeError("Export/QA failed: " + str(export_log))
        summary["exportQaPassed"] = True
    with log.with_suffix(".json").open("x") as stream:
        json.dump(summary, stream, indent=2)
        stream.write("\n")
    print("DONE {} {} input={} reco={} seconds={}".format(job["energy"], job["id"], completed, reco, summary["seconds"]), flush=True)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=Path, required=True)
    parser.add_argument("--variants", default="all")
    parser.add_argument("--workers", type=int, default=1, choices=(1, 2))
    parser.add_argument("--no-export", action="store_true")
    args = parser.parse_args()
    jobs = json.loads(args.jobs.read_text())
    if args.variants != "all":
        wanted = set(args.variants.split(","))
        jobs = [job for job in jobs if job["id"] in wanted]
        if {job["id"] for job in jobs} != wanted:
            raise RuntimeError("Unknown requested variant")
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(execute, job, args.no_export) for job in jobs]
        for future in concurrent.futures.as_completed(futures):
            future.result()


if __name__ == "__main__":
    main()
