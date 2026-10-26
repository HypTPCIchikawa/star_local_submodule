#!/usr/bin/env python3
"""Read-only replay of the independent resume audit; write only a new report."""
import argparse
import collections
import csv
import datetime
import hashlib
import itertools
import json
import math
from pathlib import Path

import yaml


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()
    project = args.project.resolve()
    base = project / "rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015"
    report_path = args.report.resolve()
    own_folder = Path(__file__).resolve().parent
    if report_path.parent != own_folder or report_path.exists():
        parser.error("Report must be a new file in the replay script's own audit directory")
    report = {
        "schemaVersion": 1,
        "createdJst": datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat(),
        "project": str(project),
        "studyBase": str(base),
        "readOnlyAnalysisInputs": True,
        "sourceRootRehashed": False,
        "newFitOrEventAnalysisPerformed": False,
        "source": str(Path(__file__).resolve()),
        "sourceSha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "checks": {},
    }
    counts = {"hashChecks": 0, "bytesHashed": 0}
    failures = []

    def check(condition, label):
        report["checks"][label] = bool(condition)
        if not condition:
            failures.append(label)

    def read_json(path):
        return json.loads(path.read_text())

    def verify_hash(path, expected):
        if path.stat().st_size > 100000000:
            raise RuntimeError("Refusing unexpectedly large file: " + str(path))
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1048576), b""):
                digest.update(chunk)
                counts["bytesHashed"] += len(chunk)
        counts["hashChecks"] += 1
        if digest.hexdigest() != expected:
            failures.append("SHA mismatch: " + str(path))

    try:
        for folder_name in ("adopted_ranges_001", "summary_free_sensitivity_001"):
            folder = base / folder_name
            completed = read_json(folder / "completed.json")
            check(completed["executionComplete"] is True, folder_name + "_executionComplete")
            check(completed["productionEnabled"] is False, folder_name + "_notProductionAdoption")
            for name, expected in completed["products"].items():
                verify_hash(folder / name, expected)

        summary = read_json(base / "summary_free_sensitivity_001/summary.json")
        context = read_json(base / "summary_free_sensitivity_001/context.json")
        adopted = read_json(base / "adopted_ranges_001/adopted_ranges.json")
        runs = context["inputRunManifestSignatures"]
        for run in runs:
            folder = Path(run["folder"])
            verify_hash(folder / "completed.json", run["completedSha256"])
            for name, expected in run["criticalSha256"].items():
                verify_hash(folder / name, expected)

        expected_coverage = {
            "_".join(values)
            for values in itertools.product(
                ("d", "t", "he3", "he4"),
                ("gaussian", "lognormal"),
                ("linear", "exponential"),
                ("nominal", "mean_minus", "mean_plus", "width_minus", "width_plus"),
            )
        }
        check(set(summary["inputStudies"]) == expected_coverage, "full_4species_4models_5cases_coverage")
        check(len(runs) == 80 and len({r["folder"] for r in runs}) == 80, "80_unique_runs")
        check(
            {r["folder"] for r in runs} == set(summary["inputStudies"].values()),
            "manifest_and_summary_input_sets_equal",
        )
        completions = [read_json(Path(r["folder"]) / "completed.json") for r in runs]
        for flag in ("executionComplete", "sourceInputPreserved", "protectedSnapshotsVerified"):
            check(all(c[flag] is True for c in completions), "all_80_" + flag)
        check(all(c["nativeRows"] == 160 for c in completions), "all_80_native160")
        check(len(summary["aggregateRows"]) == 128, "128_aggregate_rows")
        check(len(summary["detailedRows"]) == 2560, "2560_model_rows")
        check(summary["declaredShapeSensitivityComplete"] is True, "declared_sensitivity_complete")
        check(summary["fullUncertaintyValidated"] is False, "not_full_uncertainty_validated")
        report["coverage"] = {
            "runs": len(runs),
            "perSpecies": dict(collections.Counter(k.split("_")[0] for k in summary["inputStudies"])),
            "nativeFits": sum(c["nativeRows"] for c in completions),
        }

        with (base / "adopted_ranges_001/adopted_ranges.tsv").open() as stream:
            tsv_rows = list(csv.DictReader(stream, delimiter="\t"))
        report["species"] = {}
        for species, domain in adopted["species"].items():
            rows = [r for r in summary["aggregateRows"] if r["species"] == species]
            keep = [
                r for r in rows
                if any(
                    r["rLo"] >= lo - 1e-9 and r["rHi"] <= hi + 1e-9
                    for lo, hi in domain["adoptedLowTpcRanges"]
                )
            ]
            denominator = sum(r["rawSelected"] for r in rows)
            numerator = sum(r["rawSelected"] for r in keep)
            low = min(
                model["purityLowerIncludingOmittedBound"]
                for r in keep for model in r["modelDiagnostics"].values()
            )
            for name, computed, saved in (
                ("baselineRawCount", denominator, domain["baselineRawCount"]),
                ("adoptedRawCount", numerator, domain["adoptedRawCount"]),
                ("rawTrackRetention", numerator / denominator, domain["rawTrackRetention"]),
                ("minModelPurityIncludingOmittedBound", low, domain["minModelPurityIncludingOmittedBound"]),
            ):
                check(math.isclose(computed, saved, rel_tol=1e-12), species + "_" + name)
            check(all(r["quality"] and r["modelsSupported"] for r in keep), species + "_kept_bins_supported")
            check(low >= 0.5, species + "_model_low_ge_50percent")
            tab_ranges = [(float(r["min"]), float(r["max"])) for r in tsv_rows if r["species"] == species]
            check(
                len(tab_ranges) == len(domain["adoptedLowTpcRanges"])
                and all(
                    math.isclose(a, c, abs_tol=1e-9) and math.isclose(b, d, abs_tol=1e-9)
                    for (a, b), (c, d) in zip(tab_ranges, domain["adoptedLowTpcRanges"])
                ),
                species + "_json_tsv_ranges_equal",
            )
            report["species"][species] = {
                "min": domain["adoptedMin"],
                "veto": domain["adoptedVeto"],
                "kept005Bins": len(keep),
                "rawTrackRetention": numerator / denominator,
                "minModelPurityIncludingOmittedBound": low,
            }

        config = project / "config"
        cut_name = "cuts/nuclearid/nuclearid_femtoLambda_auau3p85_purity50_20261015.yaml"
        main_name = "mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_purity50_20261015.yaml"
        old_name = "cuts/nuclearid/nuclearid_femtoLambda_auau3p85_allNucleiTof20261015.yaml"
        nuclear = yaml.safe_load((config / cut_name).read_text())
        previous = yaml.safe_load((config / old_name).read_text())
        mainconf = yaml.safe_load((config / main_name).read_text())
        analysis_info = yaml.safe_load((config / mainconf["analysis"]).read_text())
        check(
            {k: v for k, v in nuclear.items() if k not in {"speciesMinPOverQ", "speciesPOverQVeto"}} == previous,
            "old_fields_unchanged_except_new_min_and_veto_blocks",
        )
        for species, domain in adopted["species"].items():
            check(nuclear["speciesMinPOverQ"][species] == domain["adoptedMin"], species + "_yaml_min")
            check(nuclear["speciesPOverQVeto"][species] == domain["adoptedVeto"], species + "_yaml_veto")
        check(mainconf["nuclearid"] == cut_name, "mainconf_points_to_final_nuclear_config")
        check(
            mainconf["analysis"] == "analysis/analysis_info_auau3p85_anaFemtoLambda_4He_KFParticle.yaml",
            "canonical_analysis_info_reused",
        )
        check(
            analysis_info["starTag"]["triggerSets"] == "production_3p85GeV_fixedTarget_2021"
            and analysis_info["starTag"]["year"] == "year2021"
            and analysis_info["starTag"]["productionTag"] == "P24iy",
            "2021_P24iy_only",
        )
        check(analysis_info["analysis"]["anaName"] == "auau3p85_anaFemtoLambda_4He_KFParticle", "canonical_output_identity")
        check(all((config / value).is_file() for value in mainconf.values()), "all_mainconf_references_exist")
        check(nuclear["maxPOverQ"] == 2.5, "old_R25_cap_preserved")
        check(nuclear["m2_selection"] is False and "he4HighRigidityTof" not in nuclear, "no_old_TOF_conflicts")
        expected_windows = {
            "d": {"mass2Min": 3.268773, "mass2Max": 3.693147},
            "t": {"mass2Min": 7.2564275, "mass2Max": 8.2816925},
            "he3": {"mass2Min": 1.78285, "mass2Max": 2.06485},
            "he4": {"mass2Min": 3.0, "mass2Max": 4.0},
        }
        high_tof = nuclear["highRigidityTof"]
        check(high_tof["enabled"] is True and high_tof["thresholdPOverQ"] == 2.0, "all4_TOF_threshold_R2")
        check(high_tof["species"] == expected_windows, "all4_TOF_windows_match")
        check(high_tof["extendHe4BeyondMaxPOverQ"] is True, "only_explicit_He4_cap_extension_configured")

        frozen = read_json(base / "baseline_001/baseline.json")
        for name in (
            "StMaker/common/NuclearIdDeDxVsMom.h",
            "config/calibration/dedx_helium_auau3p85_2021P24iy_from040.yaml",
            "config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml",
            "config/cuts/nuclearid/nuclearid_femtoLambda_auau3p85_calibrated.yaml",
        ):
            verify_hash(project / name, frozen["sourceSha256"][name])
        calibration = yaml.safe_load((config / nuclear["dedxCalibration"]["calibrationFile"]).read_text())
        check(calibration == frozen["current"]["calibration"], "current_calibration_equals_frozen_study_coefficients")
        report["finalConfigSha256"] = {
            name: hashlib.sha256((config / name).read_bytes()).hexdigest()
            for name in (cut_name, main_name)
        }
    except Exception as error:
        failures.append(type(error).__name__ + ": " + str(error))

    report.update(counts)
    report["failures"] = failures
    report["passed"] = not failures
    report["limitations"] = [
        "Audit checks saved products, provenance, numeric aggregation and final configuration; no new fit or physics validation.",
        "Retention is inclusive QA TPC track-count retention for R=0.4..2, not true signal efficiency or CF pair retention.",
        "Purity is the low model envelope in 0.05 GeV/c bins, not all-systematics or statistical confidence coverage.",
        "Runtime TOF boundaries and no-reassignment are covered by the separately saved domain tests, not re-executed here.",
        "Original 10GB ROOT is not rehashed; saved per-run sourceInputPreserved/protectedSnapshotsVerified evidence is checked.",
    ]
    with report_path.open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(json.dumps({k: report[k] for k in ("passed", "hashChecks", "bytesHashed", "failures")}, indent=2))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
