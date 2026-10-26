#!/usr/bin/env python3
"""Stage one audited original three-file group without changing event order.

Only this recovery's JSON contract is accepted. All cached ROOTs must already
have passed the independent branch/full-entry audit; SHA256 binds that evidence
to the bytes used by the worker. Verified zero-entry files are logged and omitted
from the analysis list, never confused with unreadable inputs.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import stat


SOURCE_JOBID = "0BEFC64E62F19BBA534A129CABC9B106"
CHUNK_BYTES = 8 * 1024 * 1024


def regular(path):
    info = path.lstat()
    if not stat.S_ISREG(info.st_mode):
        raise ValueError("Not a regular nonsymlink file: {}".format(path))
    return info


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(CHUNK_BYTES), b""):
            result.update(block)
    return result.hexdigest()


def contained(path, directory):
    if not path.is_absolute() or ".." in path.parts:
        raise ValueError("Noncanonical path: {}".format(path))
    path.relative_to(directory)
    if path.resolve(strict=True) != path:
        raise ValueError("Symlink in input path: {}".format(path))
    return path


def stage_group(manifest_list, scratch, cache_root, output_list):
    scratch = scratch.resolve(strict=True)
    cache_root = cache_root.resolve(strict=True)
    if scratch == Path("/") or not scratch.is_dir() or not cache_root.is_dir():
        raise ValueError("Invalid scratch/cache directory")
    output_list = output_list.absolute()
    if output_list.parent != scratch or os.path.lexists(str(output_list)):
        raise ValueError("Output list must be a new direct child of scratch")
    rows = [line.split() for line in manifest_list.read_text().splitlines() if line.strip()]
    if len(rows) != 1 or len(rows[0]) != 1:
        raise ValueError("Exactly one path-only manifest token is required per job")
    manifest_path = Path(rows[0][0])
    regular(manifest_path)
    raw_manifest = manifest_path.read_bytes()
    manifest = json.loads(raw_manifest.decode("utf-8"))
    if manifest.get("schemaVersion") != 1 or manifest.get("sourceJobid") != SOURCE_JOBID:
        raise ValueError("Wrong manifest schema or source SUMS ID")
    index = manifest["originalSubjob"]
    if type(index) is not int or not 0 <= index < 26994:
        raise ValueError("Invalid original subjob index")
    inputs = manifest["inputs"]
    if not isinstance(inputs, list) or len(inputs) != 3:
        raise ValueError("Each recovery manifest must retain the original three inputs")
    names, paths = set(), set()
    for item in inputs:
        source = contained(Path(item["path"]), cache_root)
        regular(source)
        if source.name in names or source in paths:
            raise ValueError("Repeated input in original group")
        names.add(source.name)
        paths.add(source)
        if type(item["entries"]) is not int or item["entries"] < 0:
            raise ValueError("Invalid validated Entries")
        if not re.fullmatch(r"[0-9a-f]{64}", item["sha256"]):
            raise ValueError("Invalid validated SHA256")
        uri = item["uri"]
        if ("/production_3p85GeV_fixedTarget_2021/" not in uri
                or "/P24iy/2021/" not in uri
                or not source.name.startswith("st_physics_")
                or "adc" in source.name or not uri.endswith("/" + source.name)):
            raise ValueError("Input is outside the original nonADC 2021 P24iy scope")
    if not any(item["entries"] > 0 for item in inputs):
        raise ValueError("No nonempty audited PicoDst in original group")
    staged_dir = scratch / "recovery_input"
    staged_dir.mkdir()  # Never overwrite an earlier or uncertain attempt.
    selected, skipped = [], []
    for position, item in enumerate(inputs):
        source = Path(item["path"])
        before = regular(source)
        with source.open("rb") as stream:
            if stream.read(4) != b"root":
                raise ValueError("Cached input has no ROOT header: {}".format(source))
        if item["entries"] == 0:
            if digest(source) != item["sha256"]:
                raise ValueError("SHA256 mismatch for audited empty input: {}".format(source))
            skipped.append(dict(position=position, path=str(source), reason="validated_zero_entries"))
            print("RECOVERY VERIFIED EMPTY originalSubjob={} input={} path={}".format(
                index, position, source), flush=True)
        else:
            destination = staged_dir / source.name
            sha = hashlib.sha256()
            with source.open("rb") as src, destination.open("xb") as dst:
                for block in iter(lambda: src.read(CHUNK_BYTES), b""):
                    sha.update(block)
                    dst.write(block)
                dst.flush()
                os.fsync(dst.fileno())
            if sha.hexdigest() != item["sha256"] or digest(destination) != item["sha256"]:
                raise ValueError("SHA256 mismatch while staging: {}".format(source))
            selected.append(str(destination))
            print("RECOVERY STAGED originalSubjob={} input={} entries={} path={}".format(
                index, position, item["entries"], destination), flush=True)
        after = regular(source)
        if (before.st_size, before.st_mtime_ns, before.st_ino, before.st_dev) != (
                after.st_size, after.st_mtime_ns, after.st_ino, after.st_dev):
            raise ValueError("Cache input changed during staging: {}".format(source))
    with output_list.open("x") as out:
        out.write("".join(path + "\n" for path in selected))
        out.flush()
        os.fsync(out.fileno())
    report = dict(sourceJobid=SOURCE_JOBID, originalSubjob=index,
                  manifestSha256=hashlib.sha256(raw_manifest).hexdigest(),
                  manifest=manifest, stagedInputs=selected, omittedEmptyInputs=skipped,
                  expectedEntries=sum(item["entries"] for item in inputs))
    with (scratch / "recovery_staging.json").open("x") as out:
        json.dump(report, out, indent=2)
        out.write("\n")
        out.flush()
        os.fsync(out.fileno())
    print("RECOVERY GROUP " + json.dumps(report, sort_keys=True), flush=True)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest-list", type=Path, required=True)
    parser.add_argument("--scratch", type=Path, required=True)
    parser.add_argument("--cache-root", type=Path, required=True)
    parser.add_argument("--output-list", type=Path, required=True)
    args = parser.parse_args()
    stage_group(args.manifest_list, args.scratch, args.cache_root, args.output_list)


if __name__ == "__main__":
    main()
