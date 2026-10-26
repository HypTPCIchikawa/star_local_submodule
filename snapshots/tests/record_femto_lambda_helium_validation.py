#!/usr/bin/env python3
"""Record the completed 20261013 helium implementation test without changing data.

This is a task-specific validation recorder, not an analysis or submission tool.
Source/config snapshots and the manifest are new files; existing output is refused.
"""
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess

import yaml

PROJECT = Path(__file__).resolve().parents[1]
TODO = PROJECT / "Todo/20261013"
OUT = PROJECT / "rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/helium_calibrated_he3Control_20261013"
PDF_DIR = PROJECT / "rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/calibrated_selection"
CALIBRATION = PROJECT / "config/calibration/dedx_helium_auau3p85_2021P24iy_from040.yaml"
MAIN_BASE = "config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle"


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def metadata(path):
    stat = path.stat()
    return dict(path=str(path.relative_to(PROJECT)), size=stat.st_size, sha256=sha(path))


def audit(log):
    rows = [json.loads(line) for line in log.read_text().splitlines() if line.startswith('{"status"')]
    if len(rows) != 1 or rows[0]["status"] != "passed":
        raise ValueError("Missing successful output audit: " + str(log))
    return rows[0]


def main():
    manifest = OUT / "validation_manifest.json"
    snapshot = OUT / "snapshot"
    if manifest.exists() or snapshot.exists():
        raise ValueError("Refusing an existing validation manifest or snapshot")
    curve = yaml.safe_load(CALIBRATION.read_text())
    if curve["source"]["year"] != 2021 or curve["source"]["production"] != "P24iy":
        raise ValueError("Wrong dataset in the explicit calibration")
    numerical = (TODO / "validation_helium_pid_runtime_final_v2.log").read_text()
    if "PASS helium dE/dx calibration: 8840 checks" not in numerical:
        raise ValueError("Final numerical/PID validation is incomplete")
    independent = json.loads((PDF_DIR / "validation/runtime_final_001/validation.json").read_text())
    if independent["status"] != "PASS" or independent["runtime_curve_comparison"]["rows"] != 5370:
        raise ValueError("Independent runtime curve comparison is incomplete")
    current_audits = {name: audit(TODO / ("validation_helium_outputs_" + name + "_final_v2.log"))
                      for name in ("new_base", "new_tof")}
    processing = {}
    for name in ("new_base", "pid_legacy", "new_tof"):
        version="" if name=="pid_legacy" else "_v2"
        log = TODO / ("validation_helium_events_" + name + version + ".log")
        match = re.search(r"requested=(\d+) completed=(\d+) reconstructed=(\d+) status=(\d+)", log.read_text())
        if not match or tuple(map(int, match.groups())) != (1000, 1000, 965, 0):
            raise ValueError("Incomplete real-event run: " + name)
        processing[name] = dict(requested=1000, completed=1000, reconstructed=965, status=0,
                                root=metadata(OUT / (name + "_1000" + version + ".root")), log=metadata(log))
    sources = [PROJECT / name for name in (
        "Makefile", "StMaker/common/FemtoNuclearDedxCalibration.h",
        "StMaker/common/FemtoNuclearDedxCalibration.cxx", "StMaker/common/NuclearIdDeDxVsMom.h",
        "StMaker/common/FemtoLambdaLegacy.h", "StMaker/common/FemtoLambdaLegacy.cxx",
        "StMaker/StFemtoMaker/StFemtoMaker.h", "StMaker/StFemtoMaker/StFemtoMaker.cxx",
        "analysis/FemtoLambdaRunSupport.h", "analysis/anaFemtoLambda_4He.C",
        "analysis/run_anaFemtoLambda_4He.C", "tests/auditFemtoLambdaHeliumControls.C",
        "tests/test_femto_nuclear_dedx_calibration.C", "tests/run_femto_nuclear_dedx_calibration_test.C",
        "tests/test_femto_lambda_helium_control.py", "tests/test_femto_lambda_entries.py",
        "tools/he4_purity/export_runtime_calibration.py", "tools/he4_purity/PlotCalibrationWindows.C",
        "tools/he4_purity/check_runtime_calibration.py", "tests/record_femto_lambda_helium_validation.py")]
    sources.extend(PROJECT / name for name in (
        "tests/femto_lambda_helium_control.C", "tests/run_femto_lambda_helium_control.C",
        "tools/checkHistFemtoLambdaHeliumControl.C"))
    configs = {CALIBRATION}
    for suffix in ("", "_highpurity", "_tofExtended", "_pidLegacy",
                   "_highpurity_pidLegacy", "_tofExtended_pidLegacy"):
        entry = PROJECT / (MAIN_BASE + suffix + ".yaml")
        configs.add(entry)
        refs = yaml.safe_load(entry.read_text())
        configs.update(PROJECT / "config" / value for value in refs.values())
        info = yaml.safe_load((PROJECT / "config" / refs["analysis"]).read_text())
        if (info["starTag"]["collision"], info["starTag"]["year"], info["starTag"]["productionTag"]) != (
                "auau3p85", "year2021", "P24iy"):
            raise ValueError("Wrong dataset in " + str(entry))
    for path in sorted(set(sources) | configs):
        destination = snapshot / path.relative_to(PROJECT)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, destination)
        if sha(path) != sha(destination):
            raise ValueError("Source/config changed during snapshot: " + str(path))
    data = dict(schemaVersion=1, status="passed", dataset="auau3p85 / 2021 P24iy only",
        runtime="STAR SL24y / sl73_x8664_gcc485 / ROOT 5.34.38 / Singularity",
        sourceRevision=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=str(PROJECT), text=True).strip(),
        sourceState="Dirty worktree, with exact task source/config snapshots; no commit or submit",
        calibration=metadata(CALIBRATION), knotCounts={k: len(v["knots"]) for k, v in curve["species"].items()},
        build=dict(exitCode=0, log=metadata(TODO / "validation_helium_controls_build_v2.log")),
        pidValidation=dict(checks=8840, log=metadata(TODO / "validation_helium_pid_runtime_final_v2.log")),
        oldTofRegression={energy: metadata(TODO / filename) for energy, filename in (
            ("3p85_pidLegacy", "validation_he4_tof_after_calibration_v2.log"),
            ("13p5_unchanged", "validation_he4_tof_13p5_after_calibration_v2.log"))},
        processing=processing, outputAudits=current_audits,
        independent=metadata(PDF_DIR / "validation/runtime_final_001/validation.json"),
        pairClosure=metadata(TODO / "validation_femto_helium_pair_v2.log"),
        makerPolicyGuards=metadata(TODO / "validation_femto_helium_policy_guards_v2.log"),
        qaPdf=metadata(OUT / "helium_control_inputs.pdf"),
        qaLog=metadata(TODO / "validation_femto_helium_cf_input_qa_v2_english.log"),
        sources=[metadata(p) for p in sources], configs=[metadata(p) for p in sorted(configs)],
        libraries=[metadata(PROJECT / "lib" / name) for name in (
            "libStarAnaConfig.so", "libStCommon.so", "libStFemtoMaker.so", "libStKfParticleCommon.so", "libKFParticle.so")],
        input=metadata(PROJECT / "rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/recoveries/20261003/input_cache/st_physics_22123040_raw_5000023.picoDst.root"),
        limitations=["Small integration test, not physical CF closure or isotope purity validation",
                    "Calibration has 6 quality and 199 reference fit knots for each species",
                    "TPC inclusive model purity is not a joint TOF/pair-conditioned contamination coefficient",
                    "No original merged ROOT, previous fit, batch job or git remote was changed"])
    manifest.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    print("Validation manifest:", manifest)


if __name__ == "__main__":
    main()
