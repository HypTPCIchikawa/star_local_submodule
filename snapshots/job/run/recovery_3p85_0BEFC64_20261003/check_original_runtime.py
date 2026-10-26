#!/usr/bin/env python3
"""Refuse retry if analysis sources or original physics YAMLs have changed."""
import hashlib
import json
from pathlib import Path
import tarfile
import yaml

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
ANALYSIS = "auau3p85_anaFemtoLambda_4He_KFParticle"
ID = "0BEFC64E62F19BBA534A129CABC9B106"


def main():
    meta = json.loads((REPO / "job/run/runmeta" / ("runmeta_" + ANALYSIS + "_" + ID + ".json")).read_text())
    bundle = Path(meta["runtime"]["runtimeBundle"]["path"])
    if hashlib.sha256(bundle.read_bytes()).hexdigest() != meta["runtime"]["runtimeBundle"]["sha256"]:
        raise RuntimeError("Original runtime snapshot checksum mismatch")
    checked, changed = [], []
    mainconf = meta["submit"]["mainconf"]
    with tarfile.open(str(bundle), "r:gz") as archive:
        frozen_main = archive.extractfile(mainconf).read()
        used_configs = {mainconf}
        used_configs.update("config/" + v for v in yaml.safe_load(frozen_main).values())
        for item in archive:
            if not item.isfile():
                continue
            name = item.name
            source = (name.startswith(("analysis/", "include/", "StMaker/", "StRoot/", "src/"))
                      and name.endswith((".C", ".cxx", ".cc", ".cpp", ".h", ".hpp")))
            if not source and name not in used_configs:
                continue
            data = archive.extractfile(item).read()
            path = REPO / name
            if not path.is_file() or path.read_bytes() != data:
                changed.append(name)
            checked.append(name)
    report = dict(sourceJobid=ID, checkedFiles=len(checked), changedFiles=changed,
                  runtimeBundleSha256=meta["runtime"]["runtimeBundle"]["sha256"],
                  note="Sources and original referenced YAMLs compared byte-for-byte; rebuilt libraries are not claimed byte-identical.")
    print(json.dumps(report, indent=2))
    if changed:
        raise RuntimeError("Original analysis/config changed; do not submit as an identical retry")


if __name__ == "__main__":
    main()
