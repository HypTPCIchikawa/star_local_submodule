#!/usr/bin/env python3
"""Stage and fully read missing-job inputs, without modifying source PicoDst.

Only complete, validated original groups are eligible for retry. Empty-but-valid
PicoDst inputs retain their position and provenance in the group manifest.
"""
import argparse
import concurrent.futures
import datetime
import fcntl
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time
from urllib.parse import urlsplit

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
SOURCE_ID = "0BEFC64E62F19BBA534A129CABC9B106"
ANALYSIS = "auau3p85_anaFemtoLambda_4He_KFParticle"
CACHE = REPO / "rootfile" / ANALYSIS / "recoveries/20261003/input_cache"
MAINCONF = "config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml"


def now():
    return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()


def save(path, value):
    temporary = path.with_suffix(path.suffix + ".writing")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def run(command, log, timeout):
    with log.open("x") as stream:
        stream.write("COMMAND: " + repr(command) + "\n")
        stream.flush()
        try:
            result = subprocess.run(command, cwd=str(REPO), stdout=stream,
                                    stderr=subprocess.STDOUT, timeout=timeout)
            return result.returncode
        except subprocess.TimeoutExpired:
            stream.write("\nTIMEOUT\n")
            return 124


def process(item, replicas, audit_suffix, attempt_offset):
    filename = item["filename"]
    evidence = HERE / "input_audits" / filename
    evidence.mkdir(parents=True, exist_ok=True)
    summary = evidence / "result.json"
    if summary.exists():
        old = json.loads(summary.read_text())
        if old.get("status") in ("valid", "empty"):
            cached = Path(old["path"])
            if not cached.is_file() or digest(cached) != old["sha256"]:
                raise RuntimeError("Validated cache changed: " + str(cached))
            return old
        # Failed results are immutable for this attempt. A later recovery must
        # deliberately use a new scope instead of endlessly retrying sources.
        return old
    attempts = []
    # STAR data servers use dynamically assigned ports, NOT redirector :1095.
    # Deep locate returns live endpoints for exactly this original logical file.
    parsed = urlsplit(item["uri"])
    remote_path = "/" + parsed.path.lstrip("/")
    command = ["xrdfs", parsed.netloc, "locate", "-d", "-r", remote_path]
    sources = []
    location = dict(command=command)
    try:
        located = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 universal_newlines=True, timeout=30)
        location.update(returncode=located.returncode, stdout=located.stdout, stderr=located.stderr)
        if located.returncode == 0:
            for line in located.stdout.splitlines():
                fields = line.split()
                if len(fields) != 3 or fields[1:] != ["Server", "Read"]:
                    continue
                endpoint = fields[0]
                # xrdfs reports IPv4 as IPv4-compatible bracketed IPv6.
                match = re.fullmatch(r"\[::(\d+\.\d+\.\d+\.\d+)\]:(\d+)", endpoint)
                if match:
                    endpoint = match.group(1) + ":" + match.group(2)
                if not re.fullmatch(r"[A-Za-z0-9_.-]+:\d+", endpoint):
                    continue
                source = "root://" + endpoint + "/" + remote_path
                if source not in sources:
                    sources.append(source)
    except subprocess.TimeoutExpired as error:
        location["error"] = str(error)
    sources.append(item["uri"])
    location["sources"] = sources
    location["catalogReplicas"] = replicas
    save(evidence / ("locations_%02d.json" % attempt_offset), location)
    result = dict(filename=filename, uri=item["uri"], startedAtJST=now(),
                  checksumNote="Source server query checksum unsupported; local SHA256 plus full ROOT event read, and worker-copy SHA256 verification.",
                  attempts=attempts, status="invalid")
    for position, uri in enumerate(sources):
        number = position + attempt_offset
        attempt_dir = evidence / ("attempt_%02d" % number)
        if attempt_dir.exists():
            # Never overwrite a interrupted transfer or its audit evidence.
            raise RuntimeError("Unfinished attempt requires inspection: " + str(attempt_dir))
        attempt_dir.mkdir()
        staged = CACHE / (filename + ".attempt_%02d" % number)
        if staged.exists():
            raise RuntimeError("Refusing to overwrite cache candidate: " + str(staged))
        attempt = dict(uri=uri, startedAtJST=now(), staged=str(staged))
        attempts.append(attempt)
        copy_timeout = 900
        rc = run(["xrdcp", "--nopbar", uri, str(staged)],
                 attempt_dir / "copy.log", copy_timeout)
        attempt["copyExit"] = rc
        save(attempt_dir / "result.json", attempt)
        if rc != 0:
            attempt["finishedAtJST"] = now()
            save(attempt_dir / "result.json", attempt)
            continue
        audit = attempt_dir / "root_audit.json"
        macro = str(HERE / "checkRecoveryPico.C") + audit_suffix
        command = ["bash", "script/singularity_root_batch.sh", MAINCONF,
                   "root4star", "-l", "-b", "-q", '{}("{}","{}")'.format(macro, staged, audit)]
        rc = run(command, attempt_dir / "root_audit.log", 3600)
        attempt["auditExit"] = rc
        if audit.exists():
            try:
                attempt["audit"] = json.loads(audit.read_text())
            except ValueError as error:
                attempt["auditReportError"] = str(error)
                save(attempt_dir / "result.json", attempt)
                raise RuntimeError("Audit infrastructure did not produce valid JSON: " + str(audit))
        report = attempt.get("audit", {})
        if ((rc == 0 and report.get("status") == "valid") or
                (rc == 10 and report.get("status") == "empty")):
            target = CACHE / filename
            if target.exists():
                raise RuntimeError("Refusing to replace existing cache: " + str(target))
            staged.rename(target)
            result.update(status=report["status"], path=str(target),
                          entries=report["entries"], entriesRead=report["entriesRead"],
                          size=target.stat().st_size, sha256=digest(target),
                          recoveredFrom=uri)
        attempt["finishedAtJST"] = now()
        save(attempt_dir / "result.json", attempt)
        if result["status"] != "invalid":
            break
    result["finishedAtJST"] = now()
    save(summary, result)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--audit-suffix", default="+")
    parser.add_argument("--attempt-offset", type=int, default=0,
                        help="New attempt number after an explicitly inspected diagnostic failure")
    parser.add_argument("--only", nargs="*", type=int,
                        help="Original subjob numbers, for a bounded diagnostic pilot")
    args = parser.parse_args()
    if not 1 <= args.workers <= 4:
        raise RuntimeError("workers must be 1..4")
    if not 0 <= args.attempt_offset <= 99:
        raise RuntimeError("Invalid attempt offset")
    lock = (HERE / "prepare.lock").open("a")
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    data = json.loads((HERE / "catalog_replicas.json").read_text())
    if data["parentJobId"] != SOURCE_ID or data["summary"]["inputsWithErrors"]:
        raise RuntimeError("Catalog audit does not match the recovery")
    groups = data["groups"]
    if args.only is not None:
        groups = [g for g in groups if g["originalSubjob"] in args.only]
        if len(groups) != len(args.only):
            raise RuntimeError("Unknown or duplicate diagnostic group")
    lookup = {x["filename"]: x["replicas"] for x in data["catalog"]}
    items = [item for group in groups for item in group["inputs"]]
    for item in items:
        if ("/production_3p85GeV_fixedTarget_2021/" not in item["uri"] or
                "/P24iy/2021/" not in item["uri"]):
            raise RuntimeError("Out-of-scope source")
    CACHE.mkdir(parents=True, exist_ok=True)
    results = {}
    started = now()
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(process, item, lookup[item["filename"]], args.audit_suffix, args.attempt_offset): item for item in items}
        for future in concurrent.futures.as_completed(futures):
            item = futures[future]
            try:
                result = future.result()
            except Exception as error:
                result = dict(filename=item["filename"], uri=item["uri"], status="invalid", error=repr(error))
            results[item["filename"]] = result
            print("{} {}/{} {} {}".format(now(), len(results), len(items), item["filename"], result["status"]), flush=True)
            save(HERE / "prepare_progress.json", dict(startedAtJST=started, checkedAtJST=now(),
                 expected=len(items), finished=len(results),
                 valid=sum(r["status"] == "valid" for r in results.values()),
                 empty=sum(r["status"] == "empty" for r in results.values()),
                 invalid=sum(r["status"] == "invalid" for r in results.values())))
    ready, blocked = [], []
    manifests = HERE / "manifests"
    manifests.mkdir(exist_ok=True)
    for group in groups:
        checked = [results[i["filename"]] for i in group["inputs"]]
        if not all(r["status"] in ("valid", "empty") for r in checked) or not any(r["status"] == "valid" for r in checked):
            blocked.append(dict(originalSubjob=group["originalSubjob"], inputs=checked))
            continue
        record = dict(schemaVersion=1, sourceJobid=SOURCE_ID,
                      originalSubjob=group["originalSubjob"], originalList=group["originalList"],
                      inputs=[dict(path=r["path"], sha256=r["sha256"], entries=r["entries"],
                                   uri=i["uri"], order=i["order"], audit=str(HERE / "input_audits" / i["filename"] / "result.json"))
                              for i, r in zip(group["inputs"], checked)])
        path = manifests / ("group_%05d.json" % group["originalSubjob"])
        if path.exists() and json.loads(path.read_text()) != record:
            raise RuntimeError("Changed manifest: " + str(path))
        if not path.exists():
            with path.open("x") as stream:
                json.dump(record, stream, indent=2)
                stream.write("\n")
        ready.append(dict(originalSubjob=group["originalSubjob"], manifest=str(path), entries=sum(r["entries"] for r in checked)))
    report = dict(sourceJobid=SOURCE_ID, startedAtJST=started, finishedAtJST=now(),
                  requestedGroups=len(groups), ready=ready, blocked=blocked)
    report_path = HERE / ("prepare_report.json" if args.only is None else "prepare_pilot_report.json")
    save(report_path, report)
    if args.only is None:
        filelist = HERE / "manifest_tokens.list"
        with filelist.open("x") as stream:
            stream.writelines(x["manifest"] + "\n" for x in ready)
    print("READY_GROUPS={} BLOCKED_GROUPS={}".format(len(ready), len(blocked)), flush=True)
    return 0 if ready else 1


if __name__ == "__main__":
    raise SystemExit(main())
