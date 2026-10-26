#!/usr/bin/env python3
"""Audit a scoped Git-PID restore without modifying analysis configuration.

The settings and the pre-edit repository tree are mandatory inputs.  Numerical
samples are a reproducibility check, not a statement about PID validity.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess

import yaml


SPECIES = ("d", "t", "he3", "he4")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha_bytes(value):
    return hashlib.sha256(value).hexdigest()


def checked_relative(value):
    require(isinstance(value, str) and value, "A nonempty relative path is required")
    path = Path(value)
    require(not path.is_absolute() and ".." not in path.parts,
            "Only repository-relative paths are accepted: " + value)
    return path


def source_path(root, value):
    relative = checked_relative(value)
    result = (root / relative).resolve(strict=True)
    require(root == result or root in result.parents,
            "Source escapes its repository tree: " + str(value))
    require(result.is_file(), "Source is not a file: " + str(value))
    return result


def read_yaml(path):
    value = yaml.safe_load(path.read_text())
    require(isinstance(value, dict), "Expected a YAML mapping: " + str(path))
    return value


def only_restore_switch(before, after):
    """Require the entire nuclear selection to differ in one boolean only."""
    old, new = before.get("dedxCalibration"), after.get("dedxCalibration")
    require(isinstance(old, dict) and isinstance(new, dict),
            "Both cuts must have an explicit dedxCalibration mapping")
    require(old.get("enabled") is True and new.get("enabled") is False,
            "The only accepted calibration transition is explicit true -> false")
    require(set(new) <= {"enabled", "calibrationFile"},
            "A disabled calibration must not inject knots, coefficients or shapes")
    expected = dict(before)
    expected["dedxCalibration"] = dict(old, enabled=False)
    require(after == expected,
            "Nuclear cuts changed beyond dedxCalibration.enabled")


def verify_main_transition(before, after, old_nuclear, new_nuclear):
    require(before.get("nuclearid") == old_nuclear,
            "Before mainconf has the wrong nuclearid reference")
    require(after.get("nuclearid") == new_nuclear,
            "Current mainconf has the wrong nuclearid reference")
    expected = dict(before, nuclearid=new_nuclear)
    require(after == expected, "Mainconf changed beyond its nuclearid reference")
    for name, reference in after.items():
        require(isinstance(reference, str) and reference,
                "Mainconf contains a non-path reference: " + str(name))


def tables_from_header(contents):
    result = {}
    for label, symbol, count in (("mean", "kMean", 4),
                                  ("sigma", "kSigma", 4),
                                  ("domain", "kPRange", 2)):
        match = re.search(r"\b" + symbol + r"(?:\s*\[[^\]]+\])+\s*=\s*\{(.*?)\};",
                          contents, re.S)
        require(match is not None, "Missing Git PID table: " + symbol)
        rows = [[float(item.strip()) for item in row.split(",")]
                for row in re.findall(r"\{([^{}]+)\}", match.group(1))]
        require(len(rows) == 4 and all(len(row) == count for row in rows),
                "Malformed Git PID table: " + symbol)
        require(all(math.isfinite(number) for row in rows for number in row),
                "Nonfinite Git PID table: " + symbol)
        result[label] = rows
    require(all(row[3] == 0 for row in result["sigma"]),
            "Git width is no longer the required quadratic")
    require(all(0 < lo < hi for lo, hi in result["domain"]),
            "Invalid Git PID domain")
    return result


def polynomial(coefficients, x):
    result = 0.0
    for coefficient in reversed(coefficients):
        result = result * x + coefficient
    return result


def sample_curves(tables, rigidities, multipliers):
    require(isinstance(rigidities, list) and rigidities,
            "rigiditySamples must be explicit and nonempty")
    require(isinstance(multipliers, list) and multipliers,
            "dedxProbeMultipliers must be explicit and nonempty")
    require(all(type(r) in (int, float) and math.isfinite(r) and r > 0
                for r in rigidities), "Invalid sample rigidity")
    require(all(type(d) in (int, float) and math.isfinite(d)
                for d in multipliers), "Invalid dE/dx probe multiplier")
    rows = []
    for species, name in enumerate(SPECIES):
        lo, hi = tables["domain"][species]
        for rigidity in rigidities:
            # No constant-width fallback or extrapolation is accepted here.
            if not lo <= rigidity <= hi:
                rows.append(dict(species=name, rigidity=rigidity,
                                 supportedByGitDomain=False, mean=None,
                                 sigma=None, pulls=None))
                continue
            x = math.log10(rigidity)
            mean = polynomial(tables["mean"][species], x)
            sigma = polynomial(tables["sigma"][species], x)
            require(math.isfinite(mean) and mean > 0 and
                    math.isfinite(sigma) and sigma > 0,
                    "Git curve has invalid mean or width at " + name + " R=" + str(rigidity))
            probes = [mean + multiplier * sigma for multiplier in multipliers]
            require(all(math.isfinite(probe) and probe > 0 for probe in probes),
                    "The explicit mean + k*sigma probe is nonpositive")
            rows.append(dict(species=name, rigidity=rigidity,
                             supportedByGitDomain=True, mean=mean, sigma=sigma,
                             pulls=[dict(dedx=probe, pull=(probe - mean) / sigma)
                                    for probe in probes]))
    return rows


def audit(project, settings_path, before_root):
    project = Path(project).resolve(strict=True)
    before_root = Path(before_root).resolve(strict=True)
    settings_path = Path(settings_path).resolve(strict=True)
    settings = read_yaml(settings_path)
    required = {"schemaVersion", "gitRevision", "legacyHeader", "mainconfs",
                "rigiditySamples", "dedxProbeMultipliers", "comparisonTolerance"}
    require(set(settings) == required and settings["schemaVersion"] == 1,
            "Unexpected or missing restore settings fields")
    tolerances = settings["comparisonTolerance"]
    require(isinstance(tolerances, dict) and set(tolerances) == {"relative", "absolute"}
            and all(type(value) in (int, float) and math.isfinite(value) and value >= 0
                    for value in tolerances.values()),
            "Explicit finite nonnegative comparison tolerances are required")
    revision = settings["gitRevision"]
    require(isinstance(revision, str) and revision and not revision.startswith("-"),
            "Explicit Git revision is required")
    resolved_revision = subprocess.check_output(
        ["git", "rev-parse", "--verify", revision + "^{commit}"], cwd=project,
        text=True).strip()
    header_relative = str(checked_relative(settings["legacyHeader"]))
    git_header = subprocess.check_output(
        ["git", "show", resolved_revision + ":" + header_relative], cwd=project)
    header = source_path(project, header_relative)
    require(header.read_bytes() == git_header,
            "Current PID header is not byte-identical to the explicit Git revision")
    tables = tables_from_header(git_header.decode("utf-8"))
    entries = settings["mainconfs"]
    require(isinstance(entries, list) and entries, "Explicit mainconf list is required")
    seen = set()
    records = []
    protected = {str(header.relative_to(project)): sha_bytes(header.read_bytes())}
    for entry in entries:
        require(isinstance(entry, dict) and set(entry) ==
                {"mainconf", "oldNuclearid", "newNuclearid"},
                "Malformed mainconf restore entry")
        relative = str(checked_relative(entry["mainconf"]))
        require(relative not in seen, "Duplicate mainconf: " + relative)
        seen.add(relative)
        before_path = source_path(before_root, relative)
        current_path = source_path(project, relative)
        before, current = read_yaml(before_path), read_yaml(current_path)
        verify_main_transition(before, current, entry["oldNuclearid"], entry["newNuclearid"])
        before_cuts = source_path(before_root, "config/" + entry["oldNuclearid"])
        historical_cuts = source_path(project, "config/" + entry["oldNuclearid"])
        require(historical_cuts.read_bytes() == before_cuts.read_bytes(),
                "The dated original nuclear YAML was modified")
        current_cuts = source_path(project, "config/" + entry["newNuclearid"])
        only_restore_switch(read_yaml(before_cuts), read_yaml(current_cuts))
        references = {}
        for concern, ref in current.items():
            live = source_path(project, "config/" + ref)
            if concern != "nuclearid":
                old = source_path(before_root, "config/" + before[concern])
                require(live.read_bytes() == old.read_bytes(),
                        "A non-nuclear mainconf dependency changed: " + concern)
            key = str(live.relative_to(project))
            references[concern] = dict(path=key, sha256=sha_bytes(live.read_bytes()))
            protected[key] = references[concern]["sha256"]
        protected[str(historical_cuts.relative_to(project))] = sha_bytes(historical_cuts.read_bytes())
        protected[relative] = sha_bytes(current_path.read_bytes())
        records.append(dict(mainconf=relative, nuclearid=entry["newNuclearid"],
                            calibrationEnabled=False, references=references,
                            mainconfTransitionVerified=True,
                            nonPidSelectionUnchanged=True))
    return dict(schemaVersion=1, status="passed", gitRevision=resolved_revision,
                legacyHeader=header_relative, legacyHeaderSha256=sha_bytes(git_header),
                legacyHeaderByteIdentical=True,
                settingsSha256=sha_bytes(settings_path.read_bytes()),
                beforeSnapshot=str(before_root), mainconfs=records,
                samples=sample_curves(tables, settings["rigiditySamples"], settings["dedxProbeMultipliers"]),
                comparisonTolerance=tolerances,
                protectedSources=protected, runtimeProviderVerified=False,
                productionFitsExecuted=False,
                note="Static configuration/byte audit only; run the provider macro separately.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--settings", required=True, type=Path)
    parser.add_argument("--before-dir", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--oracle-tsv", type=Path)
    args = parser.parse_args()
    report = audit(args.project, args.settings, args.before_dir)
    encoded = json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    if args.output:
        with args.output.open("x") as stream:
            stream.write(encoded)
    if args.oracle_tsv:
        with args.oracle_tsv.open("x") as stream:
            stream.write("mainconf\tspecies\trigidity\tdedx\texpectedPull\n")
            for mainconf in report["mainconfs"]:
                for sample in report["samples"]:
                    if not sample["supportedByGitDomain"]:
                        continue
                    for probe in sample["pulls"]:
                        values = (mainconf["mainconf"], SPECIES.index(sample["species"]),
                                  sample["rigidity"], probe["dedx"], probe["pull"])
                        stream.write("\t".join(str(value) for value in values) + "\n")
    print(encoded, end="")


if __name__ == "__main__":
    main()
