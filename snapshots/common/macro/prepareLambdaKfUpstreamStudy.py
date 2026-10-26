#!/usr/bin/env python3
"""Generate NEW study configs. Never change defaults, build, or submit jobs."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import yaml

REPO = Path(__file__).resolve().parents[2]

def read(path):
    with Path(path).open() as stream:
        return yaml.safe_load(stream)

def write_new(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as stream:
        stream.write(content)

def dump_new(path, obj):
    write_new(path, yaml.safe_dump(obj, sort_keys=False))

def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def config_reference(path):
    # The legacy ConfigManager unconditionally prepends config/ to references.
    return os.path.relpath(path, REPO / "config")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    study = read(args.config)
    assert study["schemaVersion"] == 1
    original_main = REPO / study["baseMainconf"]
    base_main = read(original_main)
    sources = {key: (REPO / "config" / value).resolve() for key, value in base_main.items()}
    base_kf, base_info = read(sources["kf"]), read(sources["analysis"])
    guarded = set(sources.values()) | {original_main, REPO / ".current_mainconf",
        REPO / "StMaker/StLambdaMaker/StLambdaMaker.cxx",
        REPO / "StMaker/StLambdaKFParticleMaker/StLambdaKFParticleMaker.cxx",
        REPO / "StMaker/kfparticle/StPicoKFParticleInterface.cxx"}
    guarded.update((REPO / "lib").glob("*KFParticle*.so"))
    guarded.update((REPO / "lib").glob("*KfParticle*.so"))
    guarded.update(REPO / spec["centrality"] for spec in study["energies"].values())
    guards = "".join("{}  {}\n".format(sha(path), path.relative_to(REPO))
                     for path in sorted(guarded) if path.is_file())
    for energy, spec in study["energies"].items():
        out = REPO / "rootfile" / (energy + "_anaLambda_KFParticle") / study["studyTag"]
        config_dir = out / "configs"
        if config_dir.exists():
            raise RuntimeError("Refusing existing configs: " + str(config_dir))
        write_new(out / "provenance/before_SHA256SUMS", guards)
        dump_new(out / "provenance/study_configuration.yaml", study)
        dump_new(config_dir / "event.yaml", read(sources["event"]))
        dump_new(config_dir / "centrality.yaml", read(REPO / spec["centrality"]))
        write_new(config_dir / "hist.yaml", sources["hist"].read_text())
        jobs = []
        for variant, overrides in study["upstreamVariants"].items():
            cuts = dict(base_kf, **overrides)
            info = json.loads(json.dumps(base_info))
            name = energy + "_anaLambda_KFParticle_" + study["studyTag"] + "_" + variant
            info["starTag"].update(triggerSets=spec["triggerSets"], year=spec["year"], collision=energy, storage="NFS")
            info["dataset"].update(allPicoDstList=str(out / "provenance/input/scan.list"),
                                   runRange=dict(zip(("min", "max"), spec["runRange"])))
            info["analysis"].update(anaName=name, name=name, jobName=name, scratchSubdir=name,
                                   outputFileStem=name, workDir=str(REPO), nFiles="all", mode="fxtmult")
            info_path, kf_path = config_dir / (variant + "_analysis.yaml"), config_dir / (variant + "_kf.yaml")
            main_path = config_dir / (variant + "_main.yaml")
            dump_new(info_path, info)
            dump_new(kf_path, cuts)
            dump_new(main_path, {"event": config_reference(config_dir / "event.yaml"),
                "centrality": config_reference(config_dir / "centrality.yaml"), "kf": config_reference(kf_path),
                "hist": config_reference(config_dir / "hist.yaml"), "analysis": config_reference(info_path)})
            jobs.append(dict(id=variant, energy=energy, mainconf=str(main_path),
                input=str(out / "provenance/input/scan.list"), events=study["scanEvents"],
                output=str(out / "scan10k" / (variant + ".root")),
                log=str(out / "provenance/logs" / (variant + "_10k.log")),
                export=str(out / "exports" / variant), overrides=overrides))
        write_new(out / "configs/jobs.json", json.dumps(jobs, indent=2) + "\n")
        print("Prepared {}: {} configs; {}".format(energy, len(jobs), config_dir))

if __name__ == "__main__":
    main()
