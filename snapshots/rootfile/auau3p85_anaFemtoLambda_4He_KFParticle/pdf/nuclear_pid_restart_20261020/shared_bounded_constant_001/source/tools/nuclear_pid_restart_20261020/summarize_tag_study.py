#!/usr/bin/env python3
"""Summarize saved TOF-tag fits without fitting, ROOT reads or curve adoption.

The requested chi-square sensitivity levels and representative R values are
reporting coordinates only.  None is a scientific acceptance criterion.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import statistics


SPECIES = ("d", "t", "he3", "he4")
PROJECT = Path(__file__).resolve().parents[2]
SENSITIVITY_LEVELS = (2, 3, 5, 10)
REPRESENTATIVE_R = (0.875, 1.775)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def number(value):
    return type(value) in (int, float) and math.isfinite(value)


def nonnegative(value):
    return number(value) and value >= 0


def read_json(path):
    value = json.loads(path.read_text())
    require(isinstance(value, dict), "Expected JSON object: " + str(path))
    return value


def verify_protected(values, label):
    require(isinstance(values, dict) and values, "Missing protected " + label)
    verified = {}
    for value, expected in values.items():
        require(isinstance(value, str) and Path(value).is_absolute(),
                "Protected paths must be absolute: " + label)
        require(isinstance(expected, str) and re.fullmatch(r"[0-9a-f]{64}", expected),
                "Invalid protected SHA256: " + value)
        path = Path(value).resolve(strict=True)
        require(path.is_file(), "Protected input is not a file: " + value)
        actual = digest(path)
        require(actual == expected, "Protected " + label + " changed: " + value)
        verified[str(path)] = actual
    return verified


def verify_source_snapshot(values, fit_dir, project):
    """The executed frozen source, not a subsequently edited live file, is authoritative."""
    require(isinstance(values, dict) and values, "Missing protected source snapshot")
    project = Path(project).resolve(strict=True)
    snapshot_root = (Path(fit_dir) / "source").resolve(strict=True)
    verified = {}
    for value, expected in values.items():
        require(isinstance(value, str) and Path(value).is_absolute(),
                "Protected source paths must be absolute")
        require(isinstance(expected, str) and re.fullmatch(r"[0-9a-f]{64}", expected),
                "Invalid protected source SHA256: " + value)
        original = Path(value).resolve()  # The historical live file need not still exist.
        try:
            relative = original.relative_to(project)
        except ValueError:
            raise ValueError("Protected source is outside the explicit repository: " + value)
        snapshot = (snapshot_root / relative).resolve(strict=True)
        require(snapshot_root in snapshot.parents and snapshot.is_file(),
                "Protected source snapshot escapes its frozen tree: " + value)
        actual = digest(snapshot)
        require(actual == expected, "Protected source snapshot changed: " + str(snapshot))
        verified[value] = dict(snapshotPath=str(snapshot), sha256=actual,
                               authoritative="executed_snapshot", liveSourceCompared=False)
    return verified


def metric(values):
    """An empty population is unassessed, never a fabricated zero."""
    return dict(count=len(values), min=min(values) if values else None,
                median=statistics.median(values) if values else None,
                max=max(values) if values else None)


def sensitivity(values):
    return dict(population=len(values), thresholdsAreAcceptanceCriteria=False,
                countsAtOrBelow={str(level): sum(value <= level for value in values)
                                 for level in SENSITIVITY_LEVELS})


def validate_and_normalize(planned, row, split):
    required_plan = ("id", "target", "mode", "rBin", "rLo", "rHi", "r",
                     "fitLo", "fitHi", "fitWindowAxisClipped")
    require(all(key in planned for key in required_plan), "Incomplete planned slice")
    require(planned["target"] in SPECIES, "Unknown planned target species")
    require(all(number(planned[key]) for key in ("rLo", "rHi", "r", "fitLo", "fitHi"))
            and 0 < planned["rLo"] < planned["rHi"]
            and planned["rLo"] <= planned["r"] <= planned["rHi"]
            and planned["fitLo"] < planned["fitHi"], "Invalid planned slice coordinates")
    require(type(planned["fitWindowAxisClipped"]) is bool, "Missing planned clipping flag")
    for key in ("id", "target", "mode", "rBin"):
        require(row.get(key) == planned[key], "Saved fit row order/identity mismatch: " + key)
    require(row.get("curveAdoptionAllowed") is False and
            row.get("purityAvailable") is False and row.get("purity", "missing") is None,
            "A saved row permits curve/purity adoption")
    require(type(row.get("fitCompleted")) is bool, "Missing fitCompleted flag")
    normalized = dict(id=planned["id"], species=planned["target"], mode=planned["mode"],
                      rBin=planned["rBin"], r=planned["r"], rLo=planned["rLo"],
                      rHi=planned["rHi"], region="lowR" if planned["r"] < split else "highR",
                      fitWindowAxisClipped=planned["fitWindowAxisClipped"],
                      fitCompleted=row["fitCompleted"], numericalFlagsOK=False,
                      numericalOK=False, assessment="unassessed", reasons=[],
                      foreignPeakInsideFit=None, fullNativeCounts=None,
                      fitNativeCounts=None, chi2PerNdf=None, status=None,
                      originalMean=None, originalSigma=None, fittedMean=None,
                      fittedSigma=None, meanError=None, sigmaError=None,
                      engineConditionalNdf=None, fittedShapeParameterCount=None)
    if not row["fitCompleted"]:
        require(isinstance(row.get("error"), str) and row["error"],
                "An incomplete fit must record its error")
        normalized["reasons"] = ["fit_not_completed"]
        normalized["error"] = row["error"]
        return normalized
    for key in ("rLo", "rHi", "fitLo", "fitHi"):
        require(row.get(key) == planned[key], "Saved slice coordinates changed: " + key)
    required = ("fullNativeCounts", "fitNativeCounts", "nFitBins", "status", "finite",
                "pointAccepted", "amplitudeConverged", "atBoundary", "yieldCovarianceValid",
                "shapeCovarianceValid", "foreignPeakInsideFit", "engineConditionalNdf",
                "fittedShapeParameterCount", "pearson", "pearsonPerNdf", "components")
    require(all(key in row for key in required), "Completed fit has incomplete saved diagnostics")
    flags = ("finite", "pointAccepted", "amplitudeConverged", "atBoundary",
             "yieldCovarianceValid", "shapeCovarianceValid", "foreignPeakInsideFit")
    require(all(type(row[key]) is bool for key in flags), "Fit flags must be explicit booleans")
    require(isinstance(row["components"], list), "Saved components must be an array")
    components = row["components"]
    require(len(components) == len(SPECIES) and
            {component.get("species") for component in components} == set(SPECIES),
            "Saved component array must contain each species exactly once")
    for component in components:
        active = component["species"] == planned["target"]
        require(all(component.get(key) is active for key in ("active", "freeMean", "freeSigma")),
                "Saved component activity/free-shape policy differs from target-only pilot")
    component = next(item for item in components if item["species"] == planned["target"])
    required_shape = ("initialMean", "initialSigma", "mean", "sigma", "meanError", "sigmaError")
    require(all(key in component for key in required_shape), "Incomplete saved target shape")
    normalized.update(status=row["status"], foreignPeakInsideFit=row["foreignPeakInsideFit"],
                      engineConditionalNdf=row["engineConditionalNdf"],
                      fittedShapeParameterCount=row["fittedShapeParameterCount"])
    fields = (("initialMean", "originalMean"), ("initialSigma", "originalSigma"),
              ("mean", "fittedMean"), ("sigma", "fittedSigma"),
              ("meanError", "meanError"), ("sigmaError", "sigmaError"))
    valid_shape = True
    for source, destination in fields:
        value = component[source]
        valid = number(value) and (value >= 0 if source.endswith("Error") else value > 0)
        normalized[destination] = value if valid else None
        valid_shape = valid_shape and valid
    if not valid_shape:
        normalized["reasons"].append("invalid_or_partial_target_shape")
    for key in ("fullNativeCounts", "fitNativeCounts"):
        if nonnegative(row[key]):
            normalized[key] = row[key]
        else:
            normalized["reasons"].append("invalid_" + key)
    valid_chi2 = (nonnegative(row["pearson"]) and nonnegative(row["pearsonPerNdf"])
                  and number(row["engineConditionalNdf"]) and row["engineConditionalNdf"] > 0)
    if valid_chi2:
        normalized["chi2PerNdf"] = row["pearsonPerNdf"]
    else:
        normalized["reasons"].append("chi2_unassessed")
    status_ok = type(row["status"]) is int and row["status"] == 0
    normalized["numericalFlagsOK"] = status_ok and all(
        row[key] for key in ("finite", "pointAccepted", "amplitudeConverged",
                             "shapeCovarianceValid", "yieldCovarianceValid")) and not row["atBoundary"]
    if not normalized["numericalFlagsOK"]:
        normalized["reasons"].append("numerical_flags_not_all_passed")
    normalized["numericalOK"] = normalized["numericalFlagsOK"] and not normalized["reasons"]
    if normalized["numericalOK"]:
        normalized["assessment"] = "numerically_available_not_scientifically_adopted"
    return normalized


def group_summary(rows):
    completed = [row for row in rows if row["fitCompleted"]]
    numerical = [row for row in rows if row["numericalOK"]]
    chi2_all = [row["chi2PerNdf"] for row in completed if row["chi2PerNdf"] is not None]
    chi2_numerical = [row["chi2PerNdf"] for row in numerical if row["chi2PerNdf"] is not None]
    counts = [row["fullNativeCounts"] for row in completed if row["fullNativeCounts"] is not None]
    return dict(total=len(rows), fitCompleted=len(completed), numericalOK=len(numerical),
                unassessed=len(rows) - len(numerical),
                foreignPeakInsideFit=sum(row["foreignPeakInsideFit"] is True for row in completed),
                foreignPeakStatusUnknown=sum(row["foreignPeakInsideFit"] is None for row in rows),
                nativeFullCounts=metric(counts),
                chi2PerNdfNumericalOK=dict(statistics=metric(chi2_numerical),
                                           sensitivity=sensitivity(chi2_numerical)),
                chi2PerNdfAllCompleted=dict(statistics=metric(chi2_all),
                                            sensitivity=sensitivity(chi2_all)),
                scientificAcceptanceSelected=False, curveAdoptionAllowed=False)


def summarize(fit_dir, project=None):
    fit_dir = Path(fit_dir).resolve(strict=True)
    paths = {"plan": fit_dir / "plan.json", "fits": fit_dir / "fits.json",
             "completion": fit_dir / "fit.log.completed.json",
             "nativeRoot": fit_dir / "native_spectra.root"}
    hashes_before = {key: digest(path) for key, path in paths.items()}
    plan, fits, completion = (read_json(paths[key]) for key in ("plan", "fits", "completion"))
    require(plan.get("schemaVersion") == 1, "Unsupported plan schema")
    expected_completion = dict(status="completed", inputsUnchanged=True, productionEnabled=False,
                               coefficientsUpdated=False, purityAvailable=False)
    require(completion == expected_completion, "Fit completion is not the required diagnostic-only state")
    require(fits.get("rawReadback") is True and fits.get("predictionReadback") is True,
            "Saved raw/prediction ROOT readback was not verified")
    require(fits.get("productionEnabled") is False and fits.get("purityAvailable") is False,
            "Saved fits enable production or purity adoption")
    total = plan.get("totalFits")
    require(type(total) is int and total > 0 and isinstance(plan.get("slices"), list)
            and isinstance(fits.get("results"), list) and len(plan["slices"]) == total
            and len(fits["results"]) == total and fits.get("totalSlices") == total,
            "Plan/results slice count mismatch")
    ranges = plan["config"]["baselinePlots"]
    split = ranges["rHighMin"]
    require(number(split) and split > 0 and ranges["rLowMax"] == split,
            "The explicit low/high R boundary is invalid or discontinuous")
    verified_inputs = verify_protected(plan.get("inputs"), "input")
    project = PROJECT if project is None else Path(project)
    verified_sources = verify_source_snapshot(plan.get("sourceSha256"), fit_dir, project)
    require("inputFingerprint" in plan and isinstance(plan.get("cache"), str),
            "Plan does not record its input fingerprint/cache")
    rows = [validate_and_normalize(planned, row, split)
            for planned, row in zip(plan["slices"], fits["results"])]
    identities = [(row["id"], row["mode"], row["species"], row["rBin"]) for row in rows]
    require(len(set(identities)) == total, "Duplicate planned/saved slice identity")
    by_species = {}
    for species in SPECIES:
        selected = [row for row in rows if row["species"] == species]
        by_species[species] = dict(allR=group_summary(selected),
                                   lowR=group_summary([row for row in selected if row["region"] == "lowR"]),
                                   highR=group_summary([row for row in selected if row["region"] == "highR"]))
    cases = [row for row in rows if row["r"] in REPRESENTATIVE_R]
    hashes_after = {key: digest(path) for key, path in paths.items()}
    require(hashes_before == hashes_after, "Fit products changed during summary")
    require(verify_protected(plan["inputs"], "input") == verified_inputs and
            verify_source_snapshot(plan["sourceSha256"], fit_dir, project) == verified_sources,
            "Protected inputs/sources changed during summary")
    return dict(schemaVersion=1, status="summary_completed", technicalExecutionCompleted=True,
                allSlicesFitCompleted=all(row["fitCompleted"] for row in rows),
                noFitPerformed=True, noRootReadbackPerformedBySummary=True,
                savedRawAndPredictionReadbackVerified=True, productionEnabled=False,
                coefficientsUpdated=False, curveAdoptionAllowed=False, purityAvailable=False,
                scientificAcceptanceSelected=False,
                lowHighBoundary=dict(value=split, low="R < boundary", high="R >= boundary"),
                chi2SensitivityLevels=list(SENSITIVITY_LEVELS),
                chi2SensitivityMeaning="Descriptive counts only; no threshold was selected for adoption.",
                initialShapeMeaning="Saved GL3-node-weighted initialMean/initialSigma; not recalculated.",
                inputFingerprint=plan["inputFingerprint"], cache=plan["cache"],
                fitDirectory=str(fit_dir), sourceFilesVerified=verified_sources,
                sourceAuthority="Executed frozen fit-dir/source tree; later live edits are not compared.",
                inputFilesVerified=verified_inputs,
                hashesAtSummaryCreation={key: dict(path=str(paths[key]), sha256=value)
                                         for key, value in hashes_after.items()},
                productHashMeaning="Current hashes at summary creation; completion did not freeze output hashes.",
                totals=group_summary(rows), species=by_species,
                representativeRRequested=list(REPRESENTATIVE_R), representativeCases=cases, rows=rows)


def save_report(report, output):
    with Path(output).open("x") as stream:
        json.dump(report, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fit-dir", required=True, type=Path)
    parser.add_argument("--output-file", required=True, type=Path)
    args = parser.parse_args()
    require(not args.output_file.exists(), "Summary output must be new")
    report = summarize(args.fit_dir)
    save_report(report, args.output_file)
    print("TAG_STUDY_SUMMARY_SAVED", args.output_file, "total=", report["totals"]["total"],
          "numericalOK=", report["totals"]["numericalOK"], "scientificAdoption=disabled", flush=True)


if __name__ == "__main__":
    main()
