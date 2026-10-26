#!/usr/bin/env python3
"""Render saved fresh-study Gaussian diagnostics; never fit or adopt PID values."""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
from typing import Any

import yaml


PROJECT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
NODE_HEADER = ["id", "species", "weight", "mean", "sigma"]
METRICS_HEADER = ["index", "id", "target", "mode", "r", "fitLo", "fitHi",
                  "mean", "sigma", "chi2Ndf", "numerical", "reason"]
SPECIES = {"d", "t", "he3", "he4"}
MODE_FOR_SPECIES = {"d": "TAG_D", "t": "TAG_T", "he3": "TAG_HE3", "he4": "TAG_D"}


def require(condition: Any, message: str) -> None:
    if not condition:
        raise ValueError(message)


def digest(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def read_json(path: Path) -> Any:
    with path.open(encoding="utf-8") as stream:
        return json.load(stream)


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def within(path: Path, parent: Path) -> bool:
    """Path.is_relative_to equivalent for the project's Python 3.8 host."""
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def under_project(value: str) -> Path:
    require(isinstance(value, str) and value, "Expected an explicit project-relative path")
    path = (PROJECT / value).resolve()
    require(within(path, PROJECT), f"Path escapes the project: {value}")
    require(path.is_file(), f"Missing project file: {path}")
    return path


def absolute_file(value: str) -> Path:
    require(isinstance(value, str) and Path(value).is_absolute(), "Input manifest paths must be absolute")
    path = Path(value).resolve()
    require(path.is_file(), f"Missing immutable input: {path}")
    return path


def finite(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def verify_inputs(plan: dict[str, Any], fit_dir: Path) -> dict[str, str]:
    inputs = plan.get("inputs")
    require(isinstance(inputs, dict) and inputs, "Missing immutable input SHA manifest")
    verified: dict[str, str] = {}
    for name, expected in inputs.items():
        path = absolute_file(name)
        require(isinstance(expected, str) and len(expected) == 64, f"Invalid SHA256 for {path}")
        require(digest(path) == expected, f"Immutable fit input changed: {path}")
        verified[str(path)] = expected
    cache = absolute_file(plan["cache"])
    require(str(cache) in verified, "Native cache is not in the immutable SHA manifest")
    source_hashes = plan.get("sourceSha256")
    require(isinstance(source_hashes, dict) and source_hashes, "Missing fit source SHA manifest")
    for name, expected in source_hashes.items():
        # Live sources may have advanced. The historical executed copy is authoritative.
        original = Path(name)
        require(original.is_absolute(), "Source manifest paths must be absolute")
        relative = original.relative_to(PROJECT)
        snapshot = (fit_dir / "source" / relative).resolve()
        require(within(snapshot, fit_dir / "source"), "Source snapshot escapes its directory")
        require(snapshot.is_file() and digest(snapshot) == expected, f"Executed fit snapshot changed: {snapshot}")
        verified[str(snapshot)] = expected
    fingerprint = plan.get("inputFingerprint")
    require(isinstance(fingerprint, dict) and "path" in fingerprint, "Missing raw-source fingerprint")
    raw_source = absolute_file(fingerprint["path"])
    size_keys = [key for key in ("size", "sizeBytes", "size_bytes", "bytes") if key in fingerprint]
    time_keys = [key for key in ("mtime_ns", "mtimeNs", "mtimeNS") if key in fingerprint]
    require(len(size_keys) == 1 and len(time_keys) == 1, "Unrecognized raw-source stat fingerprint schema")
    stat = raw_source.stat()
    require(stat.st_size == fingerprint[size_keys[0]] and stat.st_mtime_ns == fingerprint[time_keys[0]],
            f"Original ROOT source stat fingerprint changed: {raw_source}")
    for key, actual in (("inode", stat.st_ino), ("device", stat.st_dev)):
        if key in fingerprint:
            require(actual == fingerprint[key], f"Original ROOT source {key} changed: {raw_source}")
    return verified


def verify_diagnostic(plan: dict[str, Any], fits: dict[str, Any], completed: dict[str, Any]) -> None:
    require(plan.get("schemaVersion") == 1, "Unknown fit plan schema")
    for value, label in ((plan, "plan"), (completed, "fit completion")):
        for flag in ("productionEnabled", "coefficientsUpdated", "purityAvailable"):
            require(value.get(flag) is False, f"{label}: {flag} must be explicitly false")
    require(completed.get("status") == "completed" and completed.get("inputsUnchanged") is True,
            "Fit process did not complete with immutable inputs")
    require(fits.get("productionEnabled") is False and fits.get("purityAvailable") is False,
            "The fit output is not an isolated peak diagnostic")
    require(fits.get("rawReadback") is True and fits.get("predictionReadback") is True,
            "Fit native-spectrum/prediction readback did not pass")
    require(plan.get("totalFits") == len(plan.get("slices", [])) == fits.get("totalSlices") == len(fits.get("results", [])),
            "Fit output does not match planned slice count")
    require(plan["totalFits"] > 0, "Empty diagnostic study")
    seen = set()
    for planned, result in zip(plan["slices"], fits["results"]):
        require(planned["id"] not in seen, "Duplicate slice id in fit plan")
        seen.add(planned["id"])
        for key in ("id", "target", "mode", "rBin"):
            require(result.get(key) == planned[key], f"Fit row order/identity mismatch: {key}")
        require(planned["target"] in SPECIES and planned["mode"] == MODE_FOR_SPECIES[planned["target"]],
                "Initial calibration must use the matching TOF tag; He4 must use TAG_D")
        require(result.get("fitCompleted") in (True, False) and isinstance(result.get("fitCompleted"), bool),
                "Missing explicit per-row fit completion flag")
        require(result.get("curveAdoptionAllowed") is False and result.get("purityAvailable") is False
                and result.get("purity") is None, "A row attempts scientific adoption or purity certification")
        if result["fitCompleted"]:
            for key in ("rLo", "rHi", "fitLo", "fitHi"):
                require(finite(result.get(key)) and abs(result[key] - planned[key]) <= 1e-9,
                        f"Fit row geometry differs from plan: {key}")
            components = result.get("components")
            require(isinstance(components, list) and len(components) == 4, "Missing four known PID components")
            require({item["species"] for item in components} == SPECIES, "Unexpected component species")
            for item in components:
                target = item["species"] == planned["target"]
                for flag in ("active", "freeMean", "freeSigma"):
                    require(item.get(flag) is target, "Only the own-peak target may be active/free in this calibration")


def seed_shapes(plan: dict[str, Any]) -> dict[str, dict[str, tuple[float, float]]]:
    candidates = []
    for name in plan["inputs"]:
        path = Path(name)
        if path.suffix != ".tsv":
            continue
        with path.open(encoding="utf-8", newline="") as stream:
            if next(csv.reader(stream, delimiter="\t"), None) == NODE_HEADER:
                candidates.append(path)
    require(len(candidates) == 1, "Expected exactly one immutable original-PID node table")
    nodes: dict[str, dict[str, list[tuple[float, float, float]]]] = {}
    with candidates[0].open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            require(set(row) == set(NODE_HEADER), "Malformed original-PID node table")
            weight, mean, sigma = (float(row[key]) for key in ("weight", "mean", "sigma"))
            require(all(math.isfinite(value) for value in (weight, mean, sigma)) and weight > 0 and sigma > 0,
                    "Invalid original-PID node")
            require(row["species"] in SPECIES, "Unknown original-PID species")
            nodes.setdefault(row["id"], {}).setdefault(row["species"], []).append((weight, mean, sigma))
    result: dict[str, dict[str, tuple[float, float]]] = {}
    require(set(nodes) == {row["id"] for row in plan["slices"]}, "Node rows differ from planned slices")
    for key, species in nodes.items():
        require(set(species) == SPECIES, f"Missing known PID species in node row: {key}")
        result[key] = {}
        for name, values in species.items():
            weight = sum(item[0] for item in values)
            result[key][name] = (sum(w * m for w, m, _ in values) / weight,
                                 sum(w * s for w, _, s in values) / weight)
    return result


def metric_rows(plan: dict[str, Any], fits: dict[str, Any]) -> tuple[list[list[Any]], dict[str, Any]]:
    seeds = seed_shapes(plan)
    rows = []
    numerical_success = 0
    failed = 0
    foreign = 0
    for index, (planned, result) in enumerate(zip(plan["slices"], fits["results"])):
        original_mean, original_sigma = seeds[planned["id"]][planned["target"]]
        mean, sigma = original_mean, original_sigma
        chi2 = "null"
        notes = []
        model = result.get('model', fits.get('model', 'gaussian'))
        notes.append('Response=' + str(model))
        if 'storageMapping' in fits:
            notes.append('StorageMapping=' + str(fits['storageMapping']) + '; native counts unchanged')
        if model != 'gaussian' and model != 'gaussian_linear':
            notes.append('sigma is the Gaussian core width, not total RMS')
        numerical = False
        if result["fitCompleted"]:
            component = next(item for item in result["components"] if item["species"] == planned["target"])
            require(abs(component["initialMean"] - original_mean) <= 1e-8 and
                    abs(component["initialSigma"] - original_sigma) <= 1e-8,
                    "Fit initial mean/width differ from immutable original PID nodes")
            if finite(component.get("mean")) and finite(component.get("sigma")) and component["sigma"] > 0:
                mean, sigma = component["mean"], component["sigma"]
            else:
                notes.append("Invalid fitted shape; display padding uses original Git PID")
            if finite(result.get("pearsonPerNdf")):
                chi2 = result["pearsonPerNdf"]
            numerical = True
            for flag in ("finite", "pointAccepted", "amplitudeConverged", "shapeCovarianceValid", "yieldCovarianceValid"):
                if result.get(flag) is not True:
                    numerical = False
                    notes.append(f"{flag}={result.get(flag)}")
            if result.get("status") != 0:
                numerical = False
                notes.append(f"status={result.get('status')}")
            if result.get("atBoundary") is not False:
                numerical = False
                notes.append(f"atBoundary={result.get('atBoundary')}")
            if result.get("foreignPeakInsideFit") is True:
                foreign += 1
                notes.append("Another original-PID peak lies inside target-only fit range")
            if mean == original_mean and sigma == original_sigma and any("Invalid fitted shape" in note for note in notes):
                numerical = False
        else:
            failed += 1
            notes.append("Fit failed: " + str(result.get("error", "unspecified failure")))
            notes.append("Display mean/width are original Git PID, not a fitted measurement")
        if numerical:
            numerical_success += 1
        notes.append("Scientific fit acceptance and purity are not evaluated")
        reason = "; ".join(notes).replace("\t", " ").replace("\n", " ").replace("\r", " ")
        rows.append([index, planned["id"], planned["target"], planned["mode"], planned["r"],
                     planned["fitLo"], planned["fitHi"], mean, sigma, chi2, int(numerical), reason])
    summary = {"totalRows": len(rows), "fitExecutionFailures": failed,
               "numericallySuccessfulRows": numerical_success,
               "numericallyUnsuccessfulRows": len(rows) - numerical_success,
               "foreignPeakInsideFitRows": foreign,
               "rowsBySpecies": dict(collections.Counter(row[2] for row in rows)),
               "scientificFitsAccepted": False, "purityAvailable": False,
               "productionEnabled": False, "coefficientsUpdated": False,
               "note": "Rendering/readback success is technical only, not approval of Gaussian fits, purity, or PID coefficients."}
    return rows, summary


def verify_runtime_yaml(plan: dict[str, Any], path: Path) -> None:
    """Allow later comment corrections, never a change in runtime YAML values."""
    expected = plan['sourceSha256'].get(str(path))
    require(isinstance(expected, str), 'Runtime YAML is absent from the fit source manifest')
    if digest(path) == expected:
        return
    suffix = '/source/' + str(path.relative_to(PROJECT))
    copies = [Path(name) for name, value in plan['inputs'].items()
              if name.endswith(suffix) and value == expected]
    require(len(copies) == 1 and digest(copies[0]) == expected,
            'Missing or changed historical runtime YAML snapshot: ' + str(path))
    current, original = (yaml.safe_load(item.read_text(encoding='utf-8')) for item in (path, copies[0]))
    require(isinstance(current, dict) and current == original,
            'Runtime YAML values changed since fit: ' + str(path))


def runtime_paths(plan: dict[str, Any]) -> tuple[Path, Path]:
    cfg = plan["config"]
    mainconf, wrapper = under_project(cfg["mainconf"]), under_project(cfg["wrapper"])
    verify_runtime_yaml(plan, mainconf)
    require(plan["sourceSha256"].get(str(wrapper)) == digest(wrapper), f"Runtime wrapper changed since fit: {wrapper}")
    main = yaml.safe_load(mainconf.read_text(encoding="utf-8"))
    require(isinstance(main, dict) and isinstance(main.get("nuclearid"), str), "Missing nuclearid reference in mainconf")
    nuclear = under_project(str(Path("config") / main["nuclearid"]))
    config = yaml.safe_load(nuclear.read_text(encoding="utf-8"))
    require(isinstance(config, dict), "Malformed nuclearid configuration")
    calibration = config.get("dedxCalibration")
    require(isinstance(calibration, dict) and calibration.get("enabled") is False,
            "Current runtime must explicitly disable historical interpolation calibration")
    verify_runtime_yaml(plan, nuclear)
    return mainconf, wrapper


def run_command(command: list[str], log: Path) -> None:
    with log.open("w", encoding="utf-8") as stream:
        result = subprocess.run(command, cwd=PROJECT, stdout=stream, stderr=subprocess.STDOUT, check=False)
    require(result.returncode == 0, f"Command failed ({result.returncode}); inspect {log}")


def render(fit_dir: Path, output_dir: Path, settings_path: Path) -> dict[str, Any]:
    plan_path, fits_path = fit_dir / "plan.json", fit_dir / "fits.json"
    fit_completed_path = fit_dir / "fit.log.completed.json"
    plan, fits, completed = (read_json(path) for path in (plan_path, fits_path, fit_completed_path))
    verify_diagnostic(plan, fits, completed)
    verified = verify_inputs(plan, fit_dir)
    mainconf, wrapper = runtime_paths(plan)
    root_path = fit_dir / "native_spectra.root"
    require(root_path.is_file(), "Missing saved native-spectrum ROOT output")
    for path in (plan_path, fits_path, fit_completed_path, root_path):
        verified[str(path)] = digest(path)
    settings_sha = digest(settings_path)
    settings = yaml.safe_load(settings_path.read_text(encoding="utf-8"))
    required_settings = {"viewPaddingSigma", "displayRebin", "linearHeadroom", "logHeadroom", "logMinimum"}
    require(isinstance(settings, dict) and set(settings) == required_settings, "Plot settings must explicitly contain exactly five supported keys")
    require(finite(settings["viewPaddingSigma"]) and settings["viewPaddingSigma"] >= 0,
            "viewPaddingSigma must be finite and nonnegative")
    require(type(settings["displayRebin"]) is int and settings["displayRebin"] >= 1, "displayRebin must be a positive integer")
    for key in ("linearHeadroom", "logHeadroom"):
        require(finite(settings[key]) and settings[key] > 1, f"{key} must exceed one")
    require(finite(settings["logMinimum"]) and settings["logMinimum"] > 0, "logMinimum must be positive")
    rows, summary = metric_rows(plan, fits)
    metrics = output_dir / "metrics.tsv"
    with metrics.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(METRICS_HEADER); writer.writerows(rows)
    frozen_settings = output_dir / "view_settings.yaml"
    shutil.copy2(settings_path, frozen_settings)
    require(digest(frozen_settings) == settings_sha, "Plot settings changed while being frozen")
    numeric_settings = output_dir / "view_settings.tsv"
    numeric_settings.write_text("".join(f"{key}\t{settings[key]}\n" for key in sorted(settings)), encoding="utf-8")
    source_root = output_dir / "source"
    frozen_macro = source_root / HERE.relative_to(PROJECT) / "DrawRestartFits.C"
    frozen_engine = source_root / "tools/nuclear_pid_tof/NuclearFitCommon.h"
    historical_engine = fit_dir / "source/tools/nuclear_pid_tof/NuclearFitCommon.h"
    macro = HERE / "DrawRestartFits.C"
    require(macro.is_file() and historical_engine.is_file(), "Missing renderer or historical Gaussian math source")
    frozen_macro.parent.mkdir(parents=True)
    frozen_engine.parent.mkdir(parents=True)
    shutil.copy2(macro, frozen_macro); shutil.copy2(historical_engine, frozen_engine)
    frozen_python = source_root / HERE.relative_to(PROJECT) / Path(__file__).name
    shutil.copy2(Path(__file__), frozen_python)
    pdf = output_dir / "tag_peak_fit_diagnostics.pdf"
    runner = output_dir / "run_render_tag_study.C"
    call = "DrawRestartFits(" + ",".join(json.dumps(str(path)) for path in
                                          (root_path, metrics, numeric_settings, pdf)) + ");"
    runner.write_text("void run_render_tag_study() {\n"
                      "  if (gSystem->Load(\"libMinuit\") < 0) { gSystem->Exit(1); return; }\n"
                      f"  if (!gSystem->CompileMacro({json.dumps(str(frozen_macro))}, \"k\")) {{ gSystem->Exit(1); return; }}\n"
                      f"  gROOT->ProcessLine({json.dumps(call)});\n"
                      "}\n", encoding="utf-8")
    for path in (metrics, frozen_settings, numeric_settings, frozen_macro, frozen_engine, frozen_python, runner):
        verified[str(path)] = digest(path)
    command = ["bash", str(wrapper), str(mainconf), "root4star", "-l", "-b", "-q", str(runner) + "()"]
    require(not any("\n" in part or "\r" in part for part in command), "Control character in runtime path")
    write_json(output_dir / "command.json", command)
    render_log = output_dir / "render.log"
    run_command(command, render_log)
    text = render_log.read_text(encoding="utf-8", errors="replace")
    marker = re.findall(r"RESTART_VIEW_COMPLETE rows=(\d+) pages=(\d+) panels=(\d+) expectedPages=(\d+) refit=0 normalization=0 purityAssessed=0", text)
    expected_pages = (len(rows) + 1) // 2
    require(marker == [(str(len(rows)), str(expected_pages), str(2 * len(rows)), str(expected_pages))],
            "Missing or inconsistent renderer completion marker")
    panels = re.findall(r"RESTART_VIEW_PANEL serial=(\d+) id=(\S+) target=(\S+) mode=(\S+) log=([01]) xLo=(\S+) xHi=(\S+) displayRebin=(\d+) prediction=([01])", text)
    require(len(panels) == 2 * len(rows), "Renderer panel count mismatch")
    for panel, row in zip(panels, [row for row in rows for _ in range(2)]):
        serial, identity, target, mode, logarithmic, lo, hi, rebin, prediction = panel
        expected_index = 2 * row[0] + int(logarithmic)
        require(int(serial) == expected_index and identity == row[1] and target == row[2] and mode == row[3],
                "Renderer panel order/sample mismatch")
        require(int(rebin) == settings["displayRebin"] and float(lo) < float(hi), "Invalid local panel/rebin log")
    require(pdf.is_file() and pdf.stat().st_size > 0, "Missing rendered PDF")
    run_command(["gs", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=nullpage", "-f", str(pdf)],
                output_dir / "pdf_readback.log")
    pdf_readback = (output_dir / "pdf_readback.log").read_text(encoding="utf-8", errors="replace")
    ranges = re.findall(r"^Processing pages (\d+) through (\d+)\.$", pdf_readback, flags=re.MULTILINE)
    require(ranges == [("1", str(expected_pages))], "Ghostscript PDF page-count readback mismatch")
    read_pages = re.findall(r"^Page\s+(\d+)\s*$", pdf_readback, flags=re.MULTILINE)
    require(read_pages == [str(page) for page in range(1, expected_pages + 1)],
            "Ghostscript did not report successful interpretation of every PDF page")
    verify_inputs(plan, fit_dir)
    for name, expected in verified.items():
        require(digest(Path(name)) == expected, f"Input changed during rendering: {name}")
    runtime_paths(plan)
    summary.update({"status": "completed", "pages": expected_pages, "panels": len(panels),
                    "pdfReadback": True, "pdfPageCountVerifiedBy": "ghostscript",
                    "inputSha256Unchanged": True,
                    "fitDirectory": str(fit_dir), "outputDirectory": str(output_dir)})
    write_json(output_dir / "summary.json", summary)
    outputs = [pdf, metrics, frozen_settings, numeric_settings, render_log, runner,
               output_dir / "command.json", output_dir / "pdf_readback.log",
               output_dir / "summary.json", frozen_macro, frozen_engine, frozen_python]
    write_json(output_dir / "provenance.json", {"inputSha256": verified,
               "outputSha256": {str(path): digest(path) for path in outputs},
               "historicalFitSourcesVerified": True, "liveUnrelatedSourceEqualityRequired": False,
               "rendererSettingsSource": str(settings_path), "rendererSettingsSha256": settings_sha})
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fit-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--settings", required=True, type=Path)
    args = parser.parse_args()
    created = False
    output_dir = args.output_dir.resolve()
    try:
        fit_dir = args.fit_dir.resolve()
        require(fit_dir.is_dir(), "Fit directory does not exist")
        require(output_dir.parent == fit_dir and output_dir != fit_dir,
                "Output must be a fresh direct child of the fit directory")
        require(not output_dir.exists(), "Output directory already exists; historical views will not be overwritten")
        settings_path = args.settings.resolve()
        require(settings_path.is_file(), "Explicit plot settings file is missing")
        output_dir.mkdir(); created = True
        summary = render(fit_dir, output_dir, settings_path)
        write_json(output_dir / "completed.json", {"status": "completed", "technicalRenderingOnly": True,
                   "inputsUnchanged": True, "productionEnabled": False, "coefficientsUpdated": False,
                   "purityAvailable": False, "scientificFitsAccepted": False,
                   "pdf": str(output_dir / "tag_peak_fit_diagnostics.pdf")})
        print(json.dumps(summary, ensure_ascii=False, allow_nan=False))
        return 0
    except Exception as error:
        if created:
            write_json(output_dir / "completed.json", {"status": "failed", "error": str(error),
                       "productionEnabled": False, "coefficientsUpdated": False, "purityAvailable": False})
        print(f"RESTART_RENDER_FAILED: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
