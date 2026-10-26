#!/usr/bin/env python3
"""Preserve only validated hadd inputs in subjobs/<SUMS id>, without clobbering.

Run after watch_job_merge_notify.py has published a validated merge. This does
not merge, delete ROOT data, email, or rewrite historical provenance records.
An immutable relocation plan makes an interrupted hardlink/unlink resumable.
"""
import argparse
import datetime
import fcntl
import json
import os
from pathlib import Path
import re
import stat
import tempfile
import xml.etree.ElementTree as ET


def now():
    return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()


def regular(path):
    info = path.lstat()
    if not stat.S_ISREG(info.st_mode):
        raise ValueError("Not a regular nonsymlink file: {}".format(path))
    return dict(size=info.st_size, mtimeNs=info.st_mtime_ns,
                inode=info.st_ino, device=info.st_dev)


def read_json(path):
    regular(path)
    return json.loads(path.read_text())


def contained(path, root):
    """Reject traversal and symlinks below the canonical, trusted project root."""
    if not path.is_absolute() or ".." in path.parts:
        raise ValueError("Noncanonical path: {}".format(path))
    relative = path.relative_to(root)
    current = root
    for part in relative.parts:
        current = current / part
        if current.is_symlink():
            raise ValueError("Symlink in archive path: {}".format(current))
    return path


def create_text(path, content):
    with path.open("x") as stream:
        stream.write(content)
        stream.flush()
        os.fsync(stream.fileno())


def sync_directory(path):
    descriptor = os.open(str(path), os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def save_json(path, data):
    if os.path.lexists(str(path)):
        regular(path)
    descriptor, temporary = tempfile.mkstemp(prefix=".archive-", dir=str(path.parent))
    try:
        with os.fdopen(descriptor, "w") as stream:
            json.dump(data, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, str(path))
        sync_directory(path.parent)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def verify_fingerprint(path, expected):
    actual = regular(path)
    if actual != {key: expected[key] for key in actual}:
        raise ValueError("File fingerprint changed: {}".format(path))


def prepare(runmeta):
    meta = read_json(runmeta)
    ana, jobid = meta["anaName"], meta["jobid"]
    if not re.fullmatch(r"[A-Za-z0-9_-]+", ana) or not re.fullmatch(r"[A-Fa-f0-9]{32}", jobid):
        raise ValueError("Invalid analysis name or SUMS id")
    root = Path(meta["projectRoot"]).resolve(strict=True)
    state = contained(Path(__file__).resolve().parent, root)
    result = read_json(state / "result.json")
    if result.get("phase") != "merged":
        raise ValueError("Only phase=merged may be archived")
    if Path(result["stateDirectory"]).resolve() != state:
        raise ValueError("Result stateDirectory mismatch")
    xml_path = Path(meta["submit"]["joblistSnapshot"]["path"]).resolve(strict=True)
    contained(xml_path, root)
    xml = ET.parse(xml_path).getroot()
    output = xml.find("output").attrib
    if not output["toURL"].startswith("file:"):
        raise ValueError("Not a local output")
    directory = Path(output["toURL"][5:]).resolve(strict=True)
    contained(directory, root / "rootfile")
    template = Path(output["fromScratch"]).name
    if template != ana + "_$JOBID.root":
        raise ValueError("Output filename does not match analysis scope")
    merged = directory / (ana + "_" + jobid + "_merge.root")
    merge = result["merge"]
    if Path(merge["path"]) != merged or not merge.get("validation") or merge.get("histogramCount", 0) <= 0:
        raise ValueError("Missing or wrong validated merge record")
    merged_stat = regular(merged)
    if merged_stat["size"] != merge["size"]:
        raise ValueError("Published merge size changed")
    with merged.open("rb") as stream:
        if stream.read(4) != b"root":
            raise ValueError("Published merge has no ROOT header")
    if merge.get("fingerprint"):
        verify_fingerprint(merged, merge["fingerprint"])
    # Existing merges predate publication fingerprints. The audited staging
    # file must still be the same inode as the published hardlink.
    report_path = state / "merged_audit.tsv"
    regular(report_path)
    report = report_path.read_text().splitlines()
    files = [row.split("\t")[1] for row in report if row.startswith("FILE\t")]
    if (not report or report[-1] != "DONE\t1" or len(files) != 1
            or sum(row.startswith("HIST\t") for row in report) != merge["histogramCount"]):
        raise ValueError("Incomplete merged ROOT audit evidence")
    staged = contained(Path(files[0]), directory)
    if regular(staged) != merged_stat:
        raise ValueError("Published merge differs from audited staging file")
    published = datetime.datetime.fromisoformat(merge["publishedAt"]).timestamp()
    if merged_stat["mtimeNs"] / 1.e9 > published:
        raise ValueError("Published merge modified after validation")
    selected = result["selectedIndices"]
    if (not selected or any(type(i) is not int or i < 0 for i in selected)
            or selected != sorted(set(selected))
            or any(i >= result["expectedSubjobs"] for i in selected)
            or set(selected).intersection(result.get("excludedIndices", []))):
        raise ValueError("Invalid selected subjob indices")
    inputs = read_json(state / "input_manifest.json")["inputs"]
    if len(inputs) != len(selected) or merge["inputCount"] != len(selected):
        raise ValueError("Input count mismatch")
    archive = contained(directory / "subjobs" / jobid, root)
    entries = []
    for index, entry in zip(selected, inputs):
        source = directory / (ana + "_" + jobid + "_" + str(index) + ".root")
        if entry["path"] != str(source):
            raise ValueError("Input manifest outside exact selected scope")
        entries.append(dict(index=index, originalPath=str(source),
                            archivedPath=str(archive / source.name),
                            fingerprint={key: entry[key] for key in ("size", "mtimeNs", "inode", "device")}))
    plan = dict(schemaVersion=1, anaName=ana, jobid=jobid, mergedRoot=str(merged),
                directory=str(archive), inputCount=len(entries), files=entries)
    return root, state, archive, plan


def archive_merged_subjobs(runmeta):
    root, state, archive, plan = prepare(Path(runmeta).resolve(strict=True))
    lock_path = state / "archive.lock"
    descriptor = os.open(str(lock_path), os.O_CREAT | os.O_RDWR | os.O_NOFOLLOW, 0o600)
    with os.fdopen(descriptor, "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        # Reload after the lock; no stale selection may authorize relocation.
        root, state, archive, plan = prepare(Path(runmeta).resolve(strict=True))
        manifest_path = state / "archive_manifest.json"
        resume = os.path.lexists(str(manifest_path))
        if resume and read_json(manifest_path) != plan:
            raise ValueError("Archive plan differs from current merge provenance")
        # Check ALL paths and collisions before moving even the first input.
        for entry in plan["files"]:
            source = contained(Path(entry["originalPath"]), root)
            target = contained(Path(entry["archivedPath"]), root)
            source_exists, target_exists = os.path.lexists(str(source)), os.path.lexists(str(target))
            if target_exists and not resume:
                raise FileExistsError("Archive destination collision: {}".format(target))
            if not source_exists and not target_exists:
                raise FileNotFoundError("Missing merge input: {}".format(source))
            for path, exists in ((source, source_exists), (target, target_exists)):
                if exists:
                    verify_fingerprint(path, entry["fingerprint"])
        archive.mkdir(parents=True, exist_ok=True)
        contained(archive, root)
        if archive.stat().st_dev != Path(plan["mergedRoot"]).stat().st_dev:
            raise ValueError("Archive must be on the same filesystem")
        if not resume:
            create_text(manifest_path, json.dumps(plan, ensure_ascii=False, indent=2) + "\n")
        local_manifest = archive / "archive_manifest.json"
        if os.path.lexists(str(local_manifest)):
            if read_json(local_manifest) != plan:
                raise ValueError("Archive-local manifest collision")
        else:
            create_text(local_manifest, json.dumps(plan, ensure_ascii=False, indent=2) + "\n")
        readme = archive / "README.md"
        readme_text = ("# 統合済みsubjob ROOTの保存\n\n"
                       "- 解析: `{anaName}`\n- SUMS ID: `{jobid}`\n"
                       "- 統合ROOT: `{mergedRoot}`\n- 保存対象: 検査済み統合入力 {inputCount} ファイルのみ。\n"
                       "- ROOTデータは削除していません。未統合・失敗jobの出力は移動対象外です。\n"
                       "- `archive_manifest.json` に元のパス、現在のパス、移動前fingerprintを記録しています。\n"
                       "- 本ディレクトリを統合ROOTと同時にhaddしないでください（二重加算になります）。\n").format(**plan)
        if os.path.lexists(str(readme)):
            regular(readme)
            if readme.read_text() != readme_text:
                raise ValueError("Archive README collision")
        else:
            create_text(readme, readme_text)
        sync_directory(state)
        sync_directory(archive)
        record = dict(phase="archiving", directory=str(archive), inputCount=plan["inputCount"],
                      mergedRoot=plan["mergedRoot"], manifest=str(local_manifest), checkedAt=now())
        save_json(state / "archive_result.json", record)
        try:
            for entry in plan["files"]:
                source, target = Path(entry["originalPath"]), Path(entry["archivedPath"])
                if os.path.lexists(str(source)):
                    verify_fingerprint(source, entry["fingerprint"])
                    if not os.path.lexists(str(target)):
                        os.link(str(source), str(target), follow_symlinks=False)
                    # Persist the retained name before removing its original name.
                    sync_directory(target.parent)
                    verify_fingerprint(target, entry["fingerprint"])
                    verify_fingerprint(source, entry["fingerprint"])
                    os.unlink(str(source))  # The verified hardlink retains all ROOT data.
                verify_fingerprint(target, entry["fingerprint"])
            sync_directory(Path(plan["mergedRoot"]).parent)
            record.update(phase="archived", checkedAt=now())
            save_json(state / "archive_result.json", record)
        except Exception as error:
            record.update(phase="failed", error=str(error), checkedAt=now())
            save_json(state / "archive_result.json", record)
            raise
        return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runmeta", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(archive_merged_subjobs(args.runmeta), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
