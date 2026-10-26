#!/usr/bin/env python3
"""Compare immutable response fits; choose display solutions only, never refit/adopt."""
from __future__ import annotations

import argparse
import copy
import csv
import json
import math
from pathlib import Path
import re
import shutil
import sys

import yaml

import render_tag_study as view


PROJECT, HERE = view.PROJECT, view.HERE
FAMILIES = ("gaussian", "gaussian_left_exp", "gaussian_right_exp")
SOURCE_HEADER = ("caseIndex", "model", "root", "metrics")
ALLOWED_NUMERIC_DIFFERENCES = {"responseModel", "responseAlphaStart"}
OPTIMIZER_FIELDS = ("responseOptimizerDiagnosticsAvailable", "responseMinimizeReturned",
                    "responseMinimizeConverged", "responseMigradStatus", "responseMigradEdm",
                    "responseHesseAttempted", "responseHesseSucceeded", "responseHesseStatus")


def comparable_config(config):
    result = copy.deepcopy(config)
    result.pop("selectedResponseVariant", None)
    result.pop("responseVariants", None)
    for key in ALLOWED_NUMERIC_DIFFERENCES:
        result["numeric"].pop(key, None)
    return result


def read_table(path):
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        rows = list(reader)
    view.require(rows and all(None not in row and all(value is not None for value in row.values())
                              for row in rows), "Empty/malformed TSV: " + str(path))
    return rows


def read_numeric_settings(path):
    rows = read_table(path)
    view.require(set(rows[0]) == {"key", "value"}, "Invalid numeric settings TSV")
    result = {}
    for row in rows:
        view.require(row["key"] not in result, "Duplicate numeric settings key")
        if row["key"] not in ALLOWED_NUMERIC_DIFFERENCES:
            result[row["key"]] = row["value"]
    return result


def load_run(directory):
    plan, fits, completed = (view.read_json(directory / name) for name in
                             ("plan.json", "fits.json", "fit.log.completed.json"))
    view.verify_diagnostic(plan, fits, completed)
    protected = view.verify_inputs(plan, directory)
    mainconf, wrapper = view.runtime_paths(plan)
    for name in ("plan.json", "fits.json", "fit.log.completed.json", "native_spectra.root"):
        path = directory / name
        view.require(path.is_file(), "Missing saved fit product: " + str(path))
        protected[str(path)] = view.digest(path)
    family = plan["config"]["numeric"]["responseModel"]
    alpha_start = plan["config"]["numeric"]["responseAlphaStart"]
    view.require(family in FAMILIES and view.finite(alpha_start) and alpha_start > 0,
                 "Missing explicit response family/alpha start")
    view.require(fits.get("model") == family, "Top-level saved response family differs from plan")
    for row in fits["results"]:
        view.require(row.get("model") == family, "Saved row response family differs from plan")
    metrics, rendering_summary = view.metric_rows(plan, fits)
    original_shapes = view.seed_shapes(plan)
    gates = read_table(directory / "slices.tsv")
    expected = {"id", "mode", "target", "rBin", "fitLo", "fitHi", "gateLo", "gateHi"}
    view.require(len(gates) == len(plan["slices"]) and all(set(row) == expected for row in gates),
                 "Native slice/gate table differs from plan length/schema")
    for planned, row, gate in zip(plan["slices"], fits["results"], gates):
        view.require(gate["id"] == planned["id"] and gate["target"] == planned["target"] and
                     gate["mode"] == planned["mode"] and int(gate["rBin"]) == planned["rBin"],
                     "Immutable gate table identity differs from plan")
        for key in ("fitLo", "fitHi"):
            view.require(float(gate[key]) == planned[key], "Immutable fit-window edge differs from plan")
        if row["fitCompleted"]:
            for key in ("gateLo", "gateHi"):
                view.require(row.get(key) == float(gate[key]), "Saved gate differs from immutable slice table")
    return dict(directory=directory, plan=plan, fits=fits, metrics=metrics, gates=gates,
                family=family, alphaStart=alpha_start,
                variant=plan["config"].get("selectedResponseVariant"),
                protected=protected, mainconf=mainconf, wrapper=wrapper,
                originalShapes=original_shapes,
                numeric=read_numeric_settings(directory / "settings.tsv"),
                renderingSummary=rendering_summary)


def verify_same_data(runs):
    reference = runs[0]
    for run in runs[1:]:
        for key in ("cache", "inputFingerprint", "slices", "smoke"):
            view.require(run["plan"].get(key) == reference["plan"].get(key),
                         "Response variants changed immutable plan data: " + key)
        cache = str(Path(reference["plan"]["cache"]).resolve())
        view.require(run["protected"][cache] == reference["protected"][cache], "Response cache SHA differs")
        view.require(comparable_config(run["plan"]["config"]) == comparable_config(reference["plan"]["config"]),
                     "Non-response configuration differs between variants")
        view.require(run["numeric"] == reference["numeric"], "Non-response numeric settings differ")
        view.require(run["mainconf"] == reference["mainconf"] and run["wrapper"] == reference["wrapper"],
                     "Response variants use different actual runtime entries")
        view.require(run["gates"] == reference["gates"], "Response variants changed native slices/gates")
        for name in ("nodes.tsv", "slices.tsv", "original_curves.tsv", "baseline_settings.tsv"):
            view.require(view.digest(run["directory"] / name) == view.digest(reference["directory"] / name),
                         "Response variants changed immutable original seeds/geometry: " + name)
        old_sources, new_sources = reference["plan"]["sourceSha256"], run["plan"]["sourceSha256"]
        view.require(set(old_sources) == set(new_sources), "Response fit source sets differ")
        for original, expected in old_sources.items():
            if new_sources[original] == expected:
                continue
            # Only the study YAML itself may differ in the explicitly allowed
            # response-selection fields.  C++/Python fit sources cannot differ.
            relative = Path(original).relative_to(PROJECT)
            view.require(relative.suffix in (".yaml", ".yml"), "Response engine/helper source changed")
            before = yaml.safe_load((reference["directory"] / "source" / relative).read_text())
            after = yaml.safe_load((run["directory"] / "source" / relative).read_text())
            view.require(isinstance(before, dict) and isinstance(after, dict) and "numeric" in before and
                         "numeric" in after and comparable_config(before) == comparable_config(after) ==
                         comparable_config(reference["plan"]["config"]),
                         "A non-response source/configuration changed")
    full_counts = []
    for index in range(len(reference["plan"]["slices"])):
        counts = [run["fits"]["results"][index].get("fullNativeCounts") for run in runs]
        available = [value for value in counts if view.finite(value)]
        view.require(all(value >= 0 for value in available) and
                     (not available or all(value == available[0] for value in available)),
                     "Response variants have different saved native counts")
        full_counts.append(dict(caseIndex=index, savedFullNativeCounts=available[0] if available else None,
                                countRowsAvailable=len(available), countRowsMissing=len(counts)-len(available),
                                exactRawIdentityAlsoEstablishedBy="same cache SHA/tag/R bin and verified raw native readback"))
    return full_counts


def component(row, target):
    if not row["fitCompleted"]:
        return None
    return next(item for item in row["components"] if item["species"] == target)


def finite_nonnegative(value):
    return view.finite(value) and value >= 0


def numerical_ok(run, index):
    row = run["fits"]["results"][index]
    if not bool(run["metrics"][index][10]) or not finite_nonnegative(row.get("dataDeviance")):
        return False
    if not finite_nonnegative(row.get("pearson")) or not finite_nonnegative(row.get("pearsonPerNdf")):
        return False
    if not view.finite(row.get("engineConditionalNdf")) or row["engineConditionalNdf"] <= 0:
        return False
    shape = component(row, run["plan"]["slices"][index]["target"])
    if shape is None or not all(view.finite(shape.get(key)) and shape[key] > 0
                                for key in ("mean", "sigma", "meanError", "sigmaError")):
        return False
    if run["family"] != "gaussian":
        true_flags = ("responseOptimizerDiagnosticsAvailable", "responseMinimizeReturned",
                      "responseMinimizeConverged", "responseHesseAttempted", "responseHesseSucceeded")
        if not all(row.get(key) is True for key in true_flags):
            return False
        if row.get("responseMigradStatus") != 0 or row.get("responseHesseStatus") != 0:
            return False
        if not view.finite(row.get("responseAlpha")) or row["responseAlpha"] <= 0 or not \
                view.finite(row.get("responseAlphaError")) or row["responseAlphaError"] <= 0:
            return False
    return True


def span(values):
    valid = [value for value in values if view.finite(value)]
    return dict(count=len(valid), min=min(valid) if valid else None,
                max=max(valid) if valid else None,
                span=max(valid)-min(valid) if valid else None)


def candidate_record(run, index):
    row = run["fits"]["results"][index]
    target = run["plan"]["slices"][index]["target"]
    shape = component(row, target)
    def saved(name):
        value = row.get(name)
        return value if view.finite(value) else None
    def shape_value(name):
        value = shape.get(name) if shape is not None else None
        return value if view.finite(value) else None
    deviance, nbin, ndf = saved("dataDeviance"), saved("nFitBins"), saved("engineConditionalNdf")
    count = nbin-ndf if nbin is not None and ndf is not None else None
    if count is not None and (count < 0 or count != int(count)):
        count = None
    aic = deviance+2*count if finite_nonnegative(deviance) and count is not None else None
    original = run["originalShapes"][run["plan"]["slices"][index]["id"]][target]
    return dict(fitDirectory=str(run["directory"]), variant=run["variant"], family=run["family"],
                alphaStart=run["alphaStart"], fitCompleted=row["fitCompleted"],
                numericalOK=numerical_ok(run, index), dataDeviance=deviance, AICReferenceOnly=aic,
                AICParameterCount=count, nFitBins=nbin, ndf=ndf,
                pearsonPerNdf=saved("pearsonPerNdf"), mean=shape_value("mean"),
                meanError=shape_value("meanError"), coreSigma=shape_value("sigma"),
                coreSigmaError=shape_value("sigmaError"), alpha=saved("responseAlpha"),
                alphaError=saved("responseAlphaError"), originalMean=original[0],
                originalSigma=original[1], status=row.get("status"),
                atBoundary=row.get("atBoundary"), foreignPeakInsideFit=row.get("foreignPeakInsideFit"),
                optimizerFlags={key: row.get(key) for key in OPTIMIZER_FIELDS})


def select_display(runs, index):
    candidates = [(run, candidate_record(run, index)) for run in runs]
    numerical = [(run, value) for run, value in candidates if value["numericalOK"]]
    finite = [(run, value) for run, value in candidates if value["fitCompleted"] and
              finite_nonnegative(value["dataDeviance"])]
    pool, selection = (numerical, "minimum_deviance_numerical_solution_for_display") if numerical else \
        ((finite, "minimum_finite_deviance_failed_solution_for_display") if finite else
         (candidates[:1], "raw_only_fallback"))
    selected_run, selected = min(pool, key=lambda item: item[1]["dataDeviance"]
                                 if item[1]["dataDeviance"] is not None else math.inf)
    selected = copy.deepcopy(selected)
    selected["displaySelection"] = selection
    selected["scientificAdoptionAllowed"] = False
    selected["candidates"] = [value for _, value in candidates]
    selected["seedSpans"] = {}
    for label, values in (("allSavedCandidates", [item for _, item in candidates]),
                           ("numericalCandidates", [item for _, item in numerical])):
        selected["seedSpans"][label] = {
            key: span([value[key] for value in values])
            for key in ("mean", "coreSigma", "alpha", "dataDeviance")}
    metrics = list(selected_run["metrics"][index])
    view.require(metrics[0] == index, "Selected metrics lost their real native-spectrum index")
    metrics[10] = int(selected["numericalOK"])
    # Full phase/seed diagnostics remain in comparison.tsv and summary.json.
    # Keep the plotted note short enough not to obscure the actual peak.
    notes = ["Model=" + selected_run["family"],
             "alpha=" + (format(selected["alpha"], '.4g') if selected['alpha'] is not None else 'NA') +
             "+/-" + (format(selected["alphaError"], '.3g') if selected['alphaError'] is not None else 'NA'),
             "seed=" + str(selected_run["variant"]), "core sigma; DISPLAY ONLY, not scientifically adopted"]
    if not selected['numericalOK']:
        notes.append('FAILED NUMERICAL CHECK; see comparison.tsv')
    if selection == "raw_only_fallback":
        notes.append("RAW_ONLY_FALLBACK=1")
    metrics[11] = "; ".join(notes).replace("\t", " ").replace("\r", " ").replace("\n", " ")
    return selected_run, selected, metrics


def validate_view_settings(path):
    settings = yaml.safe_load(path.read_text())
    view.require(isinstance(settings, dict) and set(settings) ==
                 {"viewPaddingSigma", "displayRebin", "linearHeadroom", "logHeadroom", "logMinimum"},
                 "Use the exact existing five-key view settings schema")
    view.require(view.finite(settings["viewPaddingSigma"]) and settings["viewPaddingSigma"] >= 0,
                 "Invalid view padding")
    view.require(type(settings["displayRebin"]) is int and settings["displayRebin"] >= 1,
                 "Invalid display rebin")
    view.require(all(view.finite(settings[key]) and settings[key] > 1
                     for key in ("linearHeadroom", "logHeadroom")) and
                 view.finite(settings["logMinimum"]) and settings["logMinimum"] > 0,
                 "Invalid explicit headroom/log-minimum settings")
    return settings


def table(path, header, rows):
    with path.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows([["null" if value is None else value for value in row] for row in rows])


def compare(fit_dirs, output_dir, settings_path, execution_state=None):
    view.require(3 <= len(fit_dirs) <= 5 and len(set(fit_dirs)) == len(fit_dirs),
                 "Supply three to five distinct completed response variants")
    parent = fit_dirs[0].parent
    view.require(all(path.parent == parent for path in fit_dirs) and output_dir.parent == parent and
                 view.within(output_dir, PROJECT) and not output_dir.exists(),
                 "Comparison must be a fresh direct sibling under the common study parent")
    runs = [load_run(path) for path in fit_dirs]
    view.require({run["family"] for run in runs} == set(FAMILIES), "All three response families are required")
    full_counts = verify_same_data(runs)
    settings = validate_view_settings(settings_path)
    protected = {str(settings_path): view.digest(settings_path)}
    for run in runs:
        for path, expected in run["protected"].items():
            view.require(path not in protected or protected[path] == expected, "Conflicting protected source SHA")
            protected[path] = expected
    output_dir.mkdir()
    if execution_state is not None:
        execution_state["outputCreated"] = True
    metrics_by_family = {family: [] for family in FAMILIES}
    sources, cases, comparison_rows = [], [], []
    metrics_paths = {family: output_dir / (family + "_metrics.tsv") for family in FAMILIES}
    for index, planned in enumerate(runs[0]["plan"]["slices"]):
        case = {key: planned[key] for key in ("id", "target", "mode", "r", "rLo", "rHi", "fitLo", "fitHi")}
        case["caseIndex"], case["families"] = index, {}
        for family in FAMILIES:
            selected_run, selected, metric = select_display([run for run in runs if run["family"] == family], index)
            metrics_by_family[family].append(metric)
            sources.append((index, family, str(selected_run["directory"] / "native_spectra.root"),
                            str(metrics_paths[family])))
            case["families"][family] = selected
            comparison_rows.append([index, planned["id"], planned["target"], planned["r"], family,
                selected["variant"], selected["alphaStart"], selected["displaySelection"], selected["numericalOK"],
                selected["originalMean"], selected["originalSigma"], selected["mean"], selected["meanError"],
                selected["coreSigma"], selected["coreSigmaError"], selected["alpha"], selected["alphaError"],
                selected["dataDeviance"], selected["pearsonPerNdf"], selected["ndf"], selected["AICParameterCount"],
                selected["AICReferenceOnly"], selected["status"], selected["atBoundary"],
                json.dumps(selected["optimizerFlags"], separators=(",", ":")),
                json.dumps(selected["seedSpans"], separators=(",", ":"))])
        cases.append(case)
    for family in FAMILIES:
        table(metrics_paths[family], view.METRICS_HEADER, metrics_by_family[family])
    sources_path = output_dir / "sources.tsv"
    table(sources_path, SOURCE_HEADER, sources)
    comparison_path = output_dir / "comparison.tsv"
    table(comparison_path, ("caseIndex", "id", "target", "r", "family", "variant", "alphaStart", "displaySelection",
          "numericalOK", "originalMean", "originalSigma", "mean", "meanError", "coreSigma", "coreSigmaError",
          "alpha", "alphaError", "dataDeviance", "chi2Ndf", "ndf", "AICParameterCount", "AICReferenceOnly",
          "status", "atBoundary", "optimizerFlags", "seedSpans"), comparison_rows)
    frozen_settings = output_dir / "view_settings.yaml"
    shutil.copy2(settings_path, frozen_settings)
    numeric_settings = output_dir / "view_settings.tsv"
    numeric_settings.write_text("".join(f"{key}\t{settings[key]}\n" for key in sorted(settings)))
    source_root = output_dir / "source"
    source_files = [HERE / "DrawResponseComparison.C", HERE / "DrawRestartFits.C",
                    HERE / "compare_response_studies.py", HERE / "render_tag_study.py"]
    frozen_sources = []
    for path in source_files:
        view.require(path.is_file(), "Missing comparison source: " + str(path))
        destination = source_root / path.relative_to(PROJECT)
        destination.parent.mkdir(parents=True, exist_ok=True)
        expected = view.digest(path);shutil.copy2(path, destination)
        view.require(view.digest(destination) == expected, "Comparison source changed while freezing")
        protected[str(path)] = expected;frozen_sources.append(destination)
    engine_relative = Path("tools/nuclear_pid_tof/NuclearFitCommon.h")
    engine = runs[0]["directory"] / "source" / engine_relative
    frozen_engine = source_root / engine_relative
    frozen_engine.parent.mkdir(parents=True, exist_ok=True);shutil.copy2(engine, frozen_engine)
    view.require(view.digest(frozen_engine) == view.digest(engine), "Historical engine snapshot copy changed")
    frozen_sources.append(frozen_engine)
    macro = source_root / (HERE / "DrawResponseComparison.C").relative_to(PROJECT)
    pdf = output_dir / "response_comparison.pdf"
    runner = output_dir / "run_compare_response_studies.C"
    build = output_dir / "build";build.mkdir()
    call = "DrawResponseComparison(" + ",".join(json.dumps(str(path)) for path in
                                                   (sources_path, numeric_settings, pdf)) + ");"
    runner.write_text("void run_compare_response_studies() {\n"
        "  if (gSystem->Load(\"libMinuit\") < 0) { gSystem->Exit(1); return; }\n"
        f"  gSystem->SetBuildDir({json.dumps(str(build))}, kFALSE);\n"
        f"  if (!gSystem->CompileMacro({json.dumps(str(macro))}, \"k\")) {{ gSystem->Exit(1); return; }}\n"
        "  Int_t error = 0;\n"
        f"  gROOT->ProcessLine({json.dumps(call)}, &error);\n"
        "  if (error) gSystem->Exit(1);\n}\n")
    for path in [sources_path, comparison_path, frozen_settings, numeric_settings, runner,
                 *metrics_paths.values(), *frozen_sources]:
        protected[str(path)] = view.digest(path)
    command = ["bash", str(runs[0]["wrapper"]), str(runs[0]["mainconf"]),
               "root4star", "-l", "-b", "-q", str(runner) + "()"]
    view.require(not any("\n" in part or "\r" in part for part in command), "Control character in runtime command")
    view.write_json(output_dir / "command.json", command)
    view.run_command(command, output_dir / "render.log")
    text = (output_dir / "render.log").read_text(errors="replace")
    marker = re.findall(r"RESPONSE_COMPARISON_COMPLETE cases=(\d+) pages=(\d+) panels=(\d+) refit=0 normalization=0", text)
    count = len(cases)
    view.require(marker == [(str(count), str(count), str(6*count))], "Comparison completion marker/count mismatch")
    view.require(pdf.is_file() and pdf.stat().st_size > 0, "Missing closed comparison PDF")
    view.run_command(["gs", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-dPDFSTOPONERROR", "-sDEVICE=nullpage", "-f", str(pdf)],
                     output_dir / "pdf_readback.log")
    text = (output_dir / "pdf_readback.log").read_text(errors="replace")
    view.require(re.findall(r"^Processing pages (\d+) through (\d+)\.$", text, re.MULTILINE) == [("1", str(count))],
                 "Ghostscript PDF page count differs")
    view.require(re.findall(r"^Page\s+(\d+)\s*$", text, re.MULTILINE) == [str(i) for i in range(1, count+1)],
                 "Ghostscript did not interpret every comparison page in order")
    for run in runs:
        view.verify_inputs(run["plan"], run["directory"]);view.runtime_paths(run["plan"])
    for path, expected in protected.items():
        view.require(view.digest(Path(path)) == expected, "Comparison input/source changed: " + path)
    report = dict(status="completed", cases=cases, caseCount=count, pages=count, panels=6*count,
        runDirectories=[str(path) for path in fit_dirs], outputDirectory=str(output_dir),
        productionEnabled=False, coefficientsUpdated=False, purityAvailable=False, scientificAdoptionAllowed=False,
        selectionMeaning="Lowest-deviance numerically valid seed is selected for display only; finite failed fits/raw are explicit fallbacks.",
        sigmaMeaning="Gaussian core width, NOT total RMS; native-R mixture broadening is not a new PID width.",
        AICMeaning="Reference dataDeviance + 2*(nFitBins-engineConditionalNdf), including active yields; not a scientific acceptance threshold.",
        identicalNativeDataVerified=True, nativeCounts=full_counts, inputsUnchanged=True,
        pdfReadback=True, refit=False, normalization=False)
    view.write_json(output_dir / "summary.json", report)
    outputs = [pdf, output_dir / "summary.json", output_dir / "render.log", output_dir / "pdf_readback.log",
               output_dir / "command.json", sources_path, comparison_path, frozen_settings, numeric_settings,
               runner, *metrics_paths.values(), *frozen_sources]
    view.write_json(output_dir / "provenance.json", dict(inputAndSourceSha256=protected,
        outputSha256={str(path): view.digest(path) for path in outputs}, historicalFitSourcesVerified=True,
        rendererSettingsSource=str(settings_path), rendererSettingsSha256=view.digest(settings_path)))
    view.write_json(output_dir / "completed.json", dict(status="completed", inputsUnchanged=True,
        productionEnabled=False, coefficientsUpdated=False, purityAvailable=False, scientificAdoptionAllowed=False,
        cases=count, pages=count, panels=6*count, refit=False, normalization=False, pdfReadback=True))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fit-dir", action="append", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--settings", required=True, type=Path)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    execution_state = {}
    try:
        result = compare([path.resolve(strict=True) for path in args.fit_dir], output,
                         args.settings.resolve(strict=True), execution_state=execution_state)
        print("RESPONSE_STUDIES_COMPARISON_SAVED", output, "cases=", result["caseCount"], flush=True)
        return 0
    except Exception as error:
        if execution_state.get("outputCreated") and output.is_dir() and \
                not (output / "completed.json").exists() and not (output / "failure.json").exists():
            view.write_json(output / "failure.json", dict(status="failed", error=str(error),
                productionEnabled=False, coefficientsUpdated=False, purityAvailable=False,
                scientificAdoptionAllowed=False))
        print("RESPONSE_STUDIES_COMPARISON_FAILED", error, file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
