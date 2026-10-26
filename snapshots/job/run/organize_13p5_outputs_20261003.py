#!/usr/bin/env python3
"""One-run naming/layout change of the already validated 13p5 partial merge."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import sys

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "script"))
import archive_merged_subjobs as arch

ANA = "auau13p5_anaFemtoLambda_4He_KFParticle"
JOBID = "1213058CE62BCF113B2E3C75ACBD7053"
OUTPUT = REPO / "rootfile" / ANA
STATE = REPO / "job/run/merge_notify" / (ANA + "_" + JOBID)
RECORD = STATE / "layout_20261003"
OLD = OUTPUT / (ANA + "_" + JOBID + "_merge.root")
NEW = OUTPUT / (ANA + "_2020P24iy_partial_merge.root")
STAGE = OUTPUT / (".merge_staging_" + ANA + "_" + JOBID)
ARCHIVE = OUTPUT / "archive/hadd_work" / JOBID
FAILED_SUFFIX = "_failed_path_20261001_1926"


def digest(path):
    sha = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            sha.update(chunk)
    return sha.hexdigest()


def verify_subjobs():
    data = arch.read_json(OUTPUT / "subjobs" / JOBID / "archive_manifest.json")
    if data["jobid"] != JOBID or data["inputCount"] != 542:
        raise RuntimeError("Unexpected subjob archive")
    for item in data["files"]:
        target = Path(item["archivedPath"])
        arch.contained(target, OUTPUT)
        arch.verify_fingerprint(target, item["fingerprint"])
    return data


def inspect():
    result = arch.read_json(STATE / "result.json")
    if (result["phase"] != "merged" or not result["partial"] or
            len(result["selectedIndices"]) != 542 or result["expectedSubjobs"] != 545 or
            result["excludedIndices"] != [149, 250, 251] or result["merge"]["histogramCount"] != 282):
        raise RuntimeError("Not the previously validated 13p5 partial merge")
    if result["merge"]["path"] != str(OLD):
        raise RuntimeError("Wrong merge provenance")
    if arch.read_json(STATE / "archive_result.json")["phase"] != "archived":
        raise RuntimeError("Original subjobs were not archived")
    verify_subjobs()
    return result


def apply(result):
    manifest = RECORD / "relocation_manifest.json"
    if manifest.exists():
        plan = arch.read_json(manifest)
    else:
        for path in (NEW, ARCHIVE, Path(str(ARCHIVE) + FAILED_SUFFIX)):
            arch.contained(path, REPO)
            if os.path.lexists(str(path)):
                raise RuntimeError("Destination already exists: " + str(path))
        fingerprint = arch.regular(OLD)
        stage_final = STAGE / "round_1_part_0.root"
        if arch.regular(stage_final) != fingerprint or fingerprint["size"] != result["merge"]["size"]:
            raise RuntimeError("Published file differs from audited stage")
        if result["merge"].get("fingerprint"):
            arch.verify_fingerprint(OLD, result["merge"]["fingerprint"])
        staged = {p.name: arch.regular(p) for p in STAGE.iterdir()}
        failed = Path(str(STAGE) + FAILED_SUFFIX)
        if len(staged) != 7 or list(failed.iterdir()):
            raise RuntimeError("Unexpected staging contents")
        plan = dict(schemaVersion=1, createdAt=arch.now(), anaName=ANA, jobid=JOBID,
                    partial=True, inputs=542, expectedSubjobs=545, excludedIndices=[149, 250, 251],
                    originalMergedRoot=str(OLD), mergedRoot=str(NEW), fingerprint=fingerprint,
                    sha256=digest(OLD), originalStagingDirectory=str(STAGE), stagingArchive=str(ARCHIVE),
                    originalFailedStagingDirectory=str(failed),
                    failedStagingArchive=str(ARCHIVE) + FAILED_SUFFIX, stagedFiles=staged,
                    historicalEvidenceSha256={n:digest(STATE/n) for n in
                        ["result.json", "status.json", "input_manifest.json", "input_audit.tsv", "merged_audit.tsv", "archive_manifest.json", "archive_result.json"]})
        arch.verify_fingerprint(OLD, fingerprint)
        arch.create_text(manifest, json.dumps(plan, ensure_ascii=False, indent=2) + "\n")
        arch.sync_directory(RECORD)
    fingerprint = plan["fingerprint"]
    if OLD.exists():
        arch.verify_fingerprint(OLD, fingerprint)
        if not os.path.lexists(str(NEW)):
            os.link(str(OLD), str(NEW), follow_symlinks=False)
        arch.verify_fingerprint(NEW, fingerprint)
        arch.sync_directory(OUTPUT)
        os.unlink(str(OLD))  # A verified hardlink retains the same ROOT data.
    arch.verify_fingerprint(NEW, fingerprint)
    ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
    for old, new in [(STAGE, ARCHIVE),
                     (Path(plan["originalFailedStagingDirectory"]), Path(plan["failedStagingArchive"]))]:
        if old.exists():
            if os.path.lexists(str(new)):
                raise RuntimeError("Staging destination collision")
            os.rename(str(old), str(new))
        if not new.is_dir() or new.is_symlink():
            raise RuntimeError("Invalid staging archive")
    arch.sync_directory(ARCHIVE.parent)
    arch.sync_directory(OUTPUT)
    for name, expected in plan["stagedFiles"].items():
        arch.verify_fingerprint(ARCHIVE / name, expected)
    if list(Path(plan["failedStagingArchive"]).iterdir()):
        raise RuntimeError("The historical failed stage should be empty")
    verify_subjobs()
    for name, sha in plan["historicalEvidenceSha256"].items():
        if digest(STATE / name) != sha:
            raise RuntimeError("Historical evidence changed")
    if sorted(p.name for p in OUTPUT.iterdir() if p.suffix == ".root") != [NEW.name]:
        raise RuntimeError("Unexpected top-level ROOT files")
    complete = dict(phase="completed", completedAt=arch.now(), partial=True,
                    mergedRoot=str(NEW), mergedSha256=plan["sha256"], size=fingerprint["size"],
                    inputCount=542, excludedIndices=[149, 250, 251],
                    subjobDirectory=str(OUTPUT / "subjobs" / JOBID),
                    stagingArchive=str(ARCHIVE), archivedStagingFiles=7,
                    relocationManifest=str(manifest), historicalEvidencePreserved=True,
                    rootDataDeleted=False, analysisResubmitted=False)
    arch.save_json(RECORD / "result.json", complete)
    print(json.dumps(complete, ensure_ascii=False, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    result = inspect()
    if not args.apply:
        print("Validated 542 archived inputs and successful partial merge; new name:", str(NEW))
        return
    RECORD.mkdir(exist_ok=True)
    with (STATE / "monitor.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        apply(result)


if __name__ == "__main__":
    main()
