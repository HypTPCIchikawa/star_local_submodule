#!/usr/bin/env python3
"""Export explicit production domains, preserving frozen low-R ranges."""
import argparse
import csv
import json
from pathlib import Path
from combine_tpc_tof import SPECIES, quality, read_table, sha, write_table


def merge(intervals):
    result = []
    for lo, hi in sorted(intervals):
        lo, hi = round(lo, 8), round(hi, 8)
        if result and abs(lo - result[-1][1]) < 1e-7:
            result[-1][1] = hi
        elif result and lo < result[-1][1]:
            raise ValueError("Overlapping accepted domains")
        else:
            result.append([lo, hi])
    return result


def complement(intervals, maximum):
    result, cursor = [], intervals[0][0]
    for lo, hi in intervals:
        if lo > cursor:
            result.append(dict(min=cursor, max=lo))
        cursor = hi
    if cursor < maximum:
        result.append(dict(min=cursor, max=maximum))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("low-table", "high-table", "high-completed", "low-ranges", "output-dir"):
        parser.add_argument("--" + name, required=True)
    for name in ("threshold", "split", "maximum"):
        parser.add_argument("--" + name, type=float, required=True)
    args = parser.parse_args()
    paths = {k: Path(getattr(args, k)).resolve() for k in ("low_table", "high_table", "high_completed", "low_ranges")}
    protected = {str(p): sha(p) for p in paths.values()}
    done = json.loads(paths["high_completed"].read_text())
    if not done["executionComplete"] or not done["protectedInputsVerified"]:
        raise ValueError("Unverified high-R study")
    product = done["products"][paths["high_table"].name]
    expected = product["sha256"] if isinstance(product, dict) else product
    if expected != sha(paths["high_table"]):
        raise ValueError("Final summary SHA mismatch")
    if not 0 < args.threshold < 1 or not args.split < args.maximum:
        raise ValueError("Invalid explicit threshold/domain")
    low, high = read_table(paths["low_table"]), read_table(paths["high_table"])
    with paths["low_ranges"].open() as stream:
        frozen = list(csv.DictReader(stream, delimiter="\t"))
    domains, decisions, result = {}, [], []
    for species in SPECIES:
        accepted_low = [(float(r["min"]), float(r["max"])) for r in frozen if r["species"] == species]
        if not accepted_low:
            raise ValueError("Missing frozen low-R domain")
        accepted_high = []
        for regime, rows in (("TPC", low), ("TOF", high)):
            subset = sorted((r for r in rows if r["species"] == species), key=lambda r: float(r["rLo"]))
            if not subset:
                raise ValueError("Missing species/regime")
            for row in subset:
                lo, hi = round(float(row["rLo"]), 8), round(float(row["rHi"]), 8)
                valid = quality(row["quality"])
                passed = valid and float(row["purityLow"]) >= args.threshold
                if regime == "TPC":
                    selected = any(lo >= a - 1e-7 and hi <= b + 1e-7 for a, b in accepted_low)
                    if selected and not passed:
                        raise ValueError("Frozen low domain fails model-low purity")
                else:
                    if lo < args.split or hi > args.maximum:
                        raise ValueError("TOF domain outside explicit bounds")
                    selected = passed
                    if selected:
                        accepted_high.append((lo, hi))
                reason = "accepted" if selected else ("unassessed" if not valid else "below_threshold_or_frozen_low_veto")
                decisions.append(dict(species=species, regime=regime, rLo=lo, rHi=hi,
                    selected=int(selected), reason=reason, purityLow=row["purityLow"], rawSelected=row["rawSelected"]))
        intervals = merge(accepted_low + accepted_high)
        domains[species] = dict(minimum=intervals[0][0], accepted=intervals,
            veto=complement(intervals, args.maximum), maximum=args.maximum)
        result.extend(dict(species=species, min=lo, max=hi) for lo, hi in intervals)
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=False)
    write_table(output / "accepted_ranges.tsv", result, ("species", "min", "max"))
    write_table(output / "bin_decisions.tsv", decisions, tuple(decisions[0]))
    (output / "production_domains.json").write_text(json.dumps(dict(
        purityThreshold=args.threshold, split=args.split, maximum=args.maximum, domains=domains,
        protectedInputs=protected, unassessedIsZeroPurity=False,
        note="Frozen low-R ranges retained; TOF requires all models supported and model-low purity >= threshold"), indent=2) + "\n")
    if any(sha(p) != value for p, value in protected.items()):
        raise ValueError("Protected source changed")
    (output / "completed.json").write_text(json.dumps(dict(executionComplete=True,
        protectedInputsVerified=True, protectedInputs=protected,
        products={p.name: sha(p) for p in output.iterdir() if p.is_file()}), indent=2) + "\n")
    print(json.dumps(domains, indent=2))


if __name__ == "__main__":
    main()
