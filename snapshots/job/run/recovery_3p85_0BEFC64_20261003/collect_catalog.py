#!/usr/bin/env python3
"""Read-only catalog replica audit for the 38 missing 2021 P24iy groups.

Does not probe replicas, read ROOT data, stage data, or submit jobs. Original
three-file group order is retained, so recovery can preserve mixing boundaries.
"""
import concurrent.futures
import datetime
import json
from pathlib import Path
import subprocess
import time


JOB_ID = "0BEFC64E62F19BBA534A129CABC9B106"
ANALYSIS = "auau3p85_anaFemtoLambda_4He_KFParticle"
MISSING = [362, 1967, 2001, 5273, 8944, 10263, 10340, 10470,
           10784, 10786, 11005, 11006, 11494, 12357, 21109,
           21385, 21400, 21888, 22032, 22258, 22762, 23201,
           23362, 23486, 23536, 23539, 23587, 23590, 23598,
           23628, 23672, 23673, 23678, 23924, 24043, 24063,
           25445, 25555]
OUTPUT_DIR = Path(__file__).resolve().parent
REPO = OUTPUT_DIR.parents[2]
OUTPUT = OUTPUT_DIR / "catalog_replicas.json"


def timestamp():
    return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()


def query(filename):
    command = ["get_file_list.pl", "-keys", "node,path,filename,events",
               "-cond", "filename=" + filename + ",production=P24iy",
               "-limit", "0"]
    started = time.monotonic()
    result = {"filename": filename, "command": command, "replicas": [], "errors": []}
    try:
        process = subprocess.run(command, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE, universal_newlines=True,
                                 timeout=30, check=False)
        result.update(returncode=process.returncode, stdout=process.stdout, stderr=process.stderr)
        if process.returncode:
            result["errors"].append("Catalog command returned nonzero status")
        for line in process.stdout.splitlines():
            if not line.strip():
                continue
            fields = line.split("::")
            if len(fields) != 4:
                result["errors"].append("Unrecognized catalog row: " + line)
                continue
            node, path, row_filename, events = fields
            if row_filename != filename or "/production_3p85GeV_fixedTarget_2021/" not in path or "/P24iy/2021/" not in path:
                result["errors"].append("Rejected out-of-scope catalog row: " + line)
                continue
            try:
                events = int(events)
            except ValueError:
                result["errors"].append("Invalid catalog events: " + line)
                continue
            result["replicas"].append(dict(node=node, path=path, filename=row_filename, catalogEvents=events))
        if not result["replicas"]:
            result["errors"].append("No in-scope 2021 P24iy replicas")
    except (subprocess.TimeoutExpired, OSError) as error:
        result["errors"].append(str(error))
    result["elapsedSeconds"] = round(time.monotonic() - started, 3)
    return result


def main():
    if OUTPUT.exists():
        raise RuntimeError("Refusing to overwrite existing audit: " + str(OUTPUT))
    groups = []
    filenames = []
    for number in MISSING:
        source = REPO / "job/run" / (ANALYSIS + JOB_ID + "_" + str(number) + ".list")
        inputs = []
        for order, line in enumerate(source.read_text().splitlines()):
            if not line.strip():
                continue
            fields = line.split()
            uri = fields[0]
            filename = uri.rsplit("/", 1)[-1]
            if "/production_3p85GeV_fixedTarget_2021/" not in uri or "/P24iy/2021/" not in uri or not filename.startswith("st_physics_22"):
                raise RuntimeError("Out-of-scope original input: " + uri)
            inputs.append(dict(order=order, uri=uri, filename=filename,
                               originalListCatalogEvents=int(fields[1]) if len(fields) > 1 else None))
            filenames.append(filename)
        if len(inputs) != 3:
            raise RuntimeError("Expected original three-file group: " + str(source))
        groups.append(dict(originalSubjob=number, originalList=str(source), inputs=inputs))
    if len(filenames) != len(set(filenames)):
        raise RuntimeError("Duplicate filenames in missing original groups")
    results = {}
    started = timestamp()
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
        futures = {executor.submit(query, name): name for name in filenames}
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            results[result["filename"]] = result
            print("catalog {}/{}: {} replicas={} errors={}".format(
                len(results), len(filenames), result["filename"],
                len(result["replicas"]), len(result["errors"])), flush=True)
    report = dict(schemaVersion=1, parentJobId=JOB_ID, analysis=ANALYSIS,
                  startedAtJST=started, finishedAtJST=timestamp(),
                  scope="2021 P24iy only; catalog entries are not measured ROOT Entries",
                  groups=groups, catalog=[results[name] for name in filenames],
                  summary=dict(originalGroups=len(groups), originalInputs=len(filenames),
                               replicaRows=sum(len(item["replicas"]) for item in results.values()),
                               inputsWithErrors=sum(bool(item["errors"]) for item in results.values())))
    with OUTPUT.open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(json.dumps(report["summary"], sort_keys=True), flush=True)
    print(str(OUTPUT), flush=True)


if __name__ == "__main__":
    main()
