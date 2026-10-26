#!/usr/bin/env python3
"""Validate the exact retry subset before submit and map SUMS indices afterward."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import yaml

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
ORIGINAL = "auau3p85_anaFemtoLambda_4He_KFParticle"
RETRY = ORIGINAL + "_recovery_20261003"
SOURCE_ID = "0BEFC64E62F19BBA534A129CABC9B106"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--submitted-jobid")
    args = parser.parse_args()
    report = json.loads((HERE / "prepare_report.json").read_text())
    original = json.loads((HERE / "catalog_replicas.json").read_text())
    groups = {g["originalSubjob"]: g for g in original["groups"]}
    assert report["sourceJobid"] == SOURCE_ID and report["requestedGroups"] == 38
    assert 1 <= len(report["ready"]) <= 38
    ready = {item["originalSubjob"]: item for item in report["ready"]}
    blocked = {item["originalSubjob"] for item in report["blocked"]}
    assert len(ready) == len(report["ready"]) and not (ready.keys() & blocked)
    assert set(ready) | blocked == set(groups)
    paths = (HERE / "manifest_tokens.list").read_text().splitlines()
    assert len(paths) == len(set(paths)) == len(ready)
    expected_paths = {item["manifest"] for item in ready.values()}
    assert set(paths) == expected_paths
    for path in paths:
        manifest = json.loads(Path(path).read_text())
        index = manifest["originalSubjob"]
        assert manifest["sourceJobid"] == SOURCE_ID and index in ready
        assert len(manifest["inputs"]) == 3
        assert [x["uri"] for x in manifest["inputs"]] == [x["uri"] for x in groups[index]["inputs"]]
        for item in manifest["inputs"]:
            evidence = json.loads(Path(item["audit"]).read_text())
            assert evidence["status"] in ("valid", "empty")
            assert all(item[k] == evidence[k] for k in ("path", "sha256", "entries", "uri"))
            assert evidence["entriesRead"] == item["entries"]
    old_main = yaml.safe_load((REPO / ("config/mainconf/main_" + ORIGINAL + ".yaml")).read_text())
    new_main = yaml.safe_load((REPO / ("config/mainconf/main_" + RETRY + ".yaml")).read_text())
    old_analysis = yaml.safe_load((REPO / "config" / old_main.pop("analysis")).read_text())
    new_analysis = yaml.safe_load((REPO / "config" / new_main.pop("analysis")).read_text())
    assert old_main == new_main
    # Only operational identity/destination and the token split count may differ.
    for config in (old_analysis, new_analysis):
        for key in ("anaName", "name", "logDir", "errDir", "jobName", "scratchSubdir", "outputFileStem", "filesPerJob"):
            config["analysis"].pop(key, None)
    assert old_analysis == new_analysis
    xml = ET.parse(REPO / ("job/joblist/joblist_" + RETRY + ".xml")).getroot()
    assert xml.get("minFilesPerProcess") == xml.get("maxFilesPerProcess") == "1"
    assert xml.find("input").get("URL") == "filelist:" + str(HERE / "manifest_tokens.list")
    assert "/recoveries/20261003/" in xml.find("output").get("toURL")
    assert hashlib.sha256((HERE / "stage_group.py").read_bytes()).hexdigest() in xml.find("command").text
    if args.submitted_jobid:
        assert re.fullmatch(r"[0-9A-Fa-f]{32}", args.submitted_jobid)
        prefix = RETRY + args.submitted_jobid
        pattern = re.compile(re.escape(prefix) + r"_(\d+)\.list$")
        actual = {int(m.group(1)): p for p in (REPO / "job/run").iterdir() if (m := pattern.fullmatch(p.name))}
        assert sorted(actual) == list(range(len(ready)))
        mapping = []
        seen = set()
        for index, path in sorted(actual.items()):
            rows = [line.split() for line in path.read_text().splitlines() if line.strip()]
            assert len(rows) == 1
            token = rows[0][0]
            if token.startswith("file:"):
                token = token[5:]
            assert token in expected_paths and token not in seen
            seen.add(token)
            manifest = json.loads(Path(token).read_text())
            mapping.append(dict(retrySubjob=index, originalSubjob=manifest["originalSubjob"],
                                manifest=token, expectedEntries=sum(x["entries"] for x in manifest["inputs"])))
        with (HERE / "retry_job_mapping.json").open("x") as stream:
            json.dump(dict(sourceJobid=SOURCE_ID, retryJobid=args.submitted_jobid,
                           mappings=mapping, blockedOriginalSubjobs=sorted(blocked)), stream, indent=2)
            stream.write("\n")
    print("Validated {} original groups for retry; blocked originals {}".format(len(ready), sorted(blocked)))


if __name__ == "__main__":
    main()
