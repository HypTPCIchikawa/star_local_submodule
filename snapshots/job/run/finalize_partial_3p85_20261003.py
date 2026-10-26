#!/usr/bin/env python3
"""One-run, explicitly approved partial-merge acceptance and recoverable layout.

No re-hadd, analysis submission, global tolerance change, or ROOT-data deletion.
The failed automated result remains unchanged; acceptance is a separate record.
"""
import argparse
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "script"))
import archive_merged_subjobs as arch

ANA = "auau3p85_anaFemtoLambda_4He_KFParticle"
JOBID = "0BEFC64E62F19BBA534A129CABC9B106"
STATE = REPO / "job/run/merge_notify" / (ANA + "_" + JOBID)
OUTPUT = REPO / "rootfile" / ANA
ACCEPT = STATE / "manual_acceptance_20261003"
FINAL = OUTPUT / (ANA + "_2021P24iy_partial_merge.root")
SUBJOBS = OUTPUT / "subjobs" / JOBID
WORK = OUTPUT / "archive/hadd_work" / JOBID
ALLOWED = {"hRefMultWeight;1", "hVz;1"}


def digest(path):
    sha = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            sha.update(chunk)
    return sha.hexdigest()


def exact_json(path, obj):
    if path.exists():
        if arch.read_json(path) != obj:
            raise RuntimeError("Immutable record mismatch: " + str(path))
    else:
        arch.create_text(path, json.dumps(obj, ensure_ascii=False, indent=2) + "\n")
        arch.sync_directory(path.parent)


def report(path):
    files, hist = [], {}
    lines = path.read_text().splitlines()
    for line in lines:
        row = line.split("\t")
        if row[0] == "FILE":
            files.append(row[1])
        elif row[0] == "HIST":
            if len(row) != 6 or row[1] in hist:
                raise RuntimeError("Invalid audit row")
            values = [float(x) for x in row[3:]]
            if not all(math.isfinite(x) for x in values):
                raise RuntimeError("Nonfinite histogram summary")
            hist[row[1]] = dict(schema=row[2], values=values)
    if not lines or lines[-1] != "DONE\t" + str(len(files)):
        raise RuntimeError("Incomplete ROOT audit")
    return files, hist


def inspect():
    result = arch.read_json(STATE / "result.json")
    if (result["phase"] != "failed" or not result["partial"] or
            result["error"] != "Merged histogram is not the input sum: hRefMultWeight;1"):
        raise RuntimeError("Not the user-approved failed partial merge")
    selected, excluded = result["selectedIndices"], result["excludedIndices"]
    if (len(selected) != 26956 or len(excluded) != 38 or
            sorted(set(selected + excluded)) != list(range(26994)) or
            len(set(selected + excluded)) != len(selected + excluded)):
        raise RuntimeError("Unexpected partial input scope")
    manifest = arch.read_json(STATE / "input_manifest.json")["inputs"]
    if len(manifest) != len(selected):
        raise RuntimeError("Input manifest count mismatch")
    expected_files, expected = report(STATE / "input_audit.tsv")
    actual_files, actual = report(STATE / "merged_audit.tsv")
    if expected_files != [x["path"] for x in manifest] or len(actual_files) != 1:
        raise RuntimeError("Audits are for different input scope")
    if set(expected) != set(actual) or len(expected) != 282:
        raise RuntimeError("Histogram names/count differ")
    deviations = []
    for name, item in expected.items():
        other = actual[name]
        if item["schema"] != other["schema"] or item["values"][0] != other["values"][0]:
            raise RuntimeError("Schema or Entries changed: " + name)
        old, new = item["values"], other["values"]
        if all(math.isclose(a, b, rel_tol=1e-7, abs_tol=1e-7) for a, b in zip(old, new)):
            continue
        relative = [abs(b - a) / max(abs(a), 1.) for a, b in zip(old, new)]
        if name not in ALLOWED or not item["schema"].startswith("TH1F:") or max(relative) > 1.5e-7:
            raise RuntimeError("Deviation outside user authorization: " + name)
        deviations.append(dict(name=name, schema=item["schema"], input=old, merged=new,
                               relativeDifference=relative,
                               percentDifference=[x * 100 for x in relative]))
    if {x["name"] for x in deviations} != ALLOWED:
        raise RuntimeError("Measured exception set changed")
    stage_file = Path(actual_files[0])
    stage_dir = OUTPUT / (".merge_staging_" + ANA + "_" + JOBID)
    if stage_file != stage_dir / "round_2_part_0.root":
        raise RuntimeError("Unexpected final staging file")
    inputs = []
    for index, item in zip(selected, manifest):
        source = OUTPUT / (ANA + "_" + JOBID + "_" + str(index) + ".root")
        if item["path"] != str(source):
            raise RuntimeError("Input outside authorized exact job")
        inputs.append(dict(index=index, originalPath=str(source),
                           archivedPath=str(SUBJOBS / source.name),
                           fingerprint={k: item[k] for k in ("size", "mtimeNs", "inode", "device")}))
    records = ["result.json", "status.json", "input_manifest.json", "input_audit.tsv", "merged_audit.tsv"]
    evidence = {n: digest(STATE / n) for n in records}
    return dict(schemaVersion=1, anaName=ANA, jobid=JOBID, partial=True,
                expectedSubjobs=26994, inputCount=len(inputs), excludedIndices=excluded,
                finalRoot=str(FINAL), subjobDirectory=str(SUBJOBS),
                originalStagingDirectory=str(stage_dir), stagingArchive=str(WORK),
                originalStagedRoot=str(stage_file), histogramCount=len(actual),
                userAcceptedMaximumPercent=0.000015, acceptedDeviations=deviations,
                validation="ROOT recursive read passed; all 282 histogram Entries/schemas match; two measured sum deviations accepted explicitly by user",
                approval="0.000015%以下ならばいいことにします。これで完成という判断をして、directoryとファイル名の整理をお願いします。",
                caution="This accepts the measured differences, not a proof that rounding is their sole cause. No physics QA certification; 38 missing jobs excluded.",
                originalEvidenceSha256=evidence, files=inputs)


def check_inputs(plan):
    # Parent directories are fixed, canonical destinations. Each leaf is then
    # lstat-checked as a nonsymlink regular file by verify_fingerprint().
    for directory in (OUTPUT, SUBJOBS, WORK.parent):
        arch.contained(directory, REPO)
    count = 0
    for entry in plan["files"]:
        source, target = Path(entry["originalPath"]), Path(entry["archivedPath"])
        found = False
        for path in (source, target):
            if os.path.lexists(str(path)):
                arch.verify_fingerprint(path, entry["fingerprint"])
                found = True
        if not found:
            raise RuntimeError("Missing recorded input: " + str(source))
        count += 1
        if count % 5000 == 0:
            print(arch.now(), "checked input fingerprints", count, flush=True)


def finalize(plan):
    acceptance_path = ACCEPT / "acceptance.json"
    plan_path = ACCEPT / "relocation_manifest.json"
    if acceptance_path.exists():
        acceptance = arch.read_json(acceptance_path)
        if arch.read_json(plan_path) != plan:
            raise RuntimeError("Cannot resume changed accepted scope")
        fingerprint = acceptance["mergedFingerprint"]
    else:
        for path in (FINAL, SUBJOBS / "archive_manifest.json", WORK):
            if os.path.lexists(str(path)):
                raise RuntimeError("Preexisting destination: " + str(path))
        stage = Path(plan["originalStagedRoot"])
        fingerprint = arch.regular(stage)
        if fingerprint["mtimeNs"] > (STATE / "merged_audit.tsv").stat().st_mtime_ns:
            raise RuntimeError("Staged ROOT changed after audit")
        stage_entries = list(stage.parent.iterdir())
        if len(stage_entries) != 274 or not all(re.fullmatch(r"round_\d+_part_\d+\.root", p.name) for p in stage_entries):
            raise RuntimeError("Unexpected staging contents")
        work_files = {p.name: arch.regular(p) for p in stage_entries}
        print(arch.now(), "hashing accepted merged ROOT", fingerprint["size"], flush=True)
        sha = digest(stage)
        arch.verify_fingerprint(stage, fingerprint)
        acceptance = dict(phase="accepted_partial", acceptedAt=arch.now(),
                          scope={k:v for k,v in plan.items() if k != "files"},
                          mergedFingerprint=fingerprint, mergedSha256=sha,
                          stagedFiles=work_files,
                          originalFailurePreserved=str(STATE / "result.json"))
        exact_json(plan_path, plan)
        exact_json(acceptance_path, acceptance)
    # Full preflight before the first move, including resumable source/target pairs.
    check_inputs(plan)
    stage = Path(plan["originalStagedRoot"])
    if not stage.exists():
        stage = WORK / stage.name
    arch.verify_fingerprint(stage, fingerprint)
    if not FINAL.exists():
        os.link(str(stage), str(FINAL), follow_symlinks=False)
        arch.sync_directory(OUTPUT)
    arch.verify_fingerprint(FINAL, fingerprint)
    SUBJOBS.mkdir(parents=True, exist_ok=True)
    exact_json(SUBJOBS / "archive_manifest.json", plan)
    print(arch.now(), "published", str(FINAL), flush=True)
    for count, entry in enumerate(plan["files"], 1):
        source, target = Path(entry["originalPath"]), Path(entry["archivedPath"])
        if os.path.lexists(str(source)):
            arch.verify_fingerprint(source, entry["fingerprint"])
            if not os.path.lexists(str(target)):
                os.link(str(source), str(target), follow_symlinks=False)
            arch.verify_fingerprint(target, entry["fingerprint"])
            arch.sync_directory(SUBJOBS)
            os.unlink(str(source))  # Verified retained hardlink preserves data.
        arch.verify_fingerprint(target, entry["fingerprint"])
        if count % 1000 == 0:
            arch.save_json(ACCEPT / "progress.json", dict(phase="archiving", moved=count, expected=len(plan["files"]), checkedAt=arch.now()))
            print(arch.now(), "archived", count, "/", len(plan["files"]), flush=True)
    arch.sync_directory(OUTPUT)
    stage_dir = Path(plan["originalStagingDirectory"])
    WORK.parent.mkdir(parents=True, exist_ok=True)
    if stage_dir.exists():
        if os.path.lexists(str(WORK)):
            raise RuntimeError("Staging archive collision")
        os.rename(str(stage_dir), str(WORK))
        arch.sync_directory(WORK.parent)
        arch.sync_directory(OUTPUT)
    for name, expected in acceptance["stagedFiles"].items():
        arch.verify_fingerprint(WORK / name, expected)
    check_inputs(plan)
    remaining = sorted(p.name for p in OUTPUT.iterdir() if p.suffix == ".root")
    if remaining != [FINAL.name]:
        raise RuntimeError("Unexpected top-level ROOTs remain")
    for name, sha in plan["originalEvidenceSha256"].items():
        if digest(STATE / name) != sha:
            raise RuntimeError("Historical evidence was modified")
    arch.verify_fingerprint(FINAL, fingerprint)
    result = dict(phase="completed", completedAt=arch.now(), partial=True,
                  inputCount=plan["inputCount"], excludedIndices=plan["excludedIndices"],
                  mergedRoot=str(FINAL), mergedSha256=acceptance["mergedSha256"],
                  size=fingerprint["size"], histogramCount=plan["histogramCount"],
                  acceptance=str(acceptance_path), subjobDirectory=str(SUBJOBS),
                  archivedSubjobs=len(plan["files"]), stagingArchive=str(WORK),
                  archivedStagingFiles=len(acceptance["stagedFiles"]),
                  originalFailurePreserved=True, rootDataDeleted=False,
                  analysisResubmitted=False, globalToleranceChanged=False)
    arch.save_json(ACCEPT / "result.json", result)
    arch.save_json(ACCEPT / "progress.json", result)
    print(json.dumps(result, ensure_ascii=False, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true", help="Apply this explicitly approved one-run acceptance")
    args = parser.parse_args()
    plan = inspect()
    if not args.apply:
        print(json.dumps({k:v for k,v in plan.items() if k != "files"}, ensure_ascii=False, indent=2))
        return
    ACCEPT.mkdir(exist_ok=True)
    # Same lock as the original watcher; never finalize alongside an active merge.
    with (STATE / "monitor.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with (ACCEPT / "finalize.lock").open("a") as second:
            fcntl.flock(second, fcntl.LOCK_EX | fcntl.LOCK_NB)
            try:
                finalize(plan)
            except Exception as error:
                arch.save_json(ACCEPT / "progress.json", dict(phase="failed", checkedAt=arch.now(), error=str(error)))
                raise


if __name__ == "__main__":
    main()
