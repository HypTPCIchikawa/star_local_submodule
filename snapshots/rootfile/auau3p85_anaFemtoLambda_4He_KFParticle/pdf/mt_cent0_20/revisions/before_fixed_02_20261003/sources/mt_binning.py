#!/usr/bin/env python3
"""Choose reproducible mT slices from saved, nonnegative TH2 counts.

This module performs no ROOT I/O, CF calculation, fitting or purity correction.
All intervals are half-open and use existing histogram edges.  ``base`` is the
loaded summary-settings mapping; ``settings`` is the loaded mT-settings mapping.
The only optimization is the squared deviation of cumulative signal-SE counts
from equally spaced quantiles. Fixed-edge mode performs no optimization at all.
CF values and method differences are never used.
"""

import math
import sys


MODES = ("SE", "ME")
REGIONS = ("signal", "leftSB", "rightSB", "signalLow", "signalHigh")


class MtFlowError(ValueError):
    """An automatic split would silently omit relevant mT flow contents."""

    def __init__(self, message, audit):
        super().__init__(message)
        self.audit = audit


def _number(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("{} must be a finite number".format(name))
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("{} must be a finite number".format(name))
    return value


def _edges(values, name):
    if not isinstance(values, (list, tuple)) or len(values) < 2:
        raise ValueError("{} must contain at least two bin edges".format(name))
    result = [_number(v, name) for v in values]
    if any(b <= a for a, b in zip(result, result[1:])):
        raise ValueError("{} must be strictly increasing".format(name))
    return result


def _same_number(a, b):
    # Numerical serialization tolerance, not an analysis/count threshold.
    return math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-14)


def _same_edges(a, b):
    return len(a) == len(b) and all(_same_number(x, y) for x, y in zip(a, b))


def _matrix(values, nx, ny, name):
    if not isinstance(values, (list, tuple)) or len(values) != nx:
        raise ValueError("{} must have {} x-major rows".format(name, nx))
    result = []
    for ix, row in enumerate(values):
        if not isinstance(row, (list, tuple)) or len(row) != ny:
            raise ValueError("{} row {} must have {} columns".format(name, ix, ny))
        parsed = [_number(v, name) for v in row]
        if any(v < 0.0 for v in parsed):
            raise ValueError("{} contains negative raw counts/variance".format(name))
        result.append(parsed)
    return result


def _edge_index(edges, value, name):
    value = _number(value, name)
    matches = [i for i, edge in enumerate(edges) if _same_number(edge, value)]
    if len(matches) != 1:
        raise ValueError("{}={} is not a unique saved bin edge".format(name, value))
    return matches[0]


def _range_indices(edges, bounds, name):
    if not isinstance(bounds, (list, tuple)) or len(bounds) != 2:
        raise ValueError("{} must contain [low, high]".format(name))
    first = _edge_index(edges, bounds[0], name + " low")
    stop = _edge_index(edges, bounds[1], name + " high")
    if stop <= first:
        raise ValueError("{} must have positive width".format(name))
    return list(range(first, stop))


def _sum_cells(channel, x_indices, y_indices, flow=False):
    counts = channel["flowCounts" if flow else "counts"]
    error2 = channel["flowError2" if flow else "error2"]
    xs, ys = list(x_indices), list(y_indices)
    return {"count": math.fsum(counts[x][y] for x in xs for y in ys),
            "error2": math.fsum(error2[x][y] for x in xs for y in ys)}


def _nonzero(summary):
    return summary["count"] != 0.0 or summary["error2"] != 0.0


def _prefix(values):
    result = [0.0]
    for value in values:
        result.append(result[-1] + value)
    return result


def _valid_interval(prefixes, first, stop):
    return first < stop and all(p[stop] - p[first] > 0.0 for p in prefixes)


def _better(candidate, previous):
    if previous is None:
        return True
    cost, path = candidate
    old_cost, old_path = previous
    # Equal costs can differ at floating-point roundoff.  Prefer lower edges,
    # deterministically, without using any physics-dependent tie-breaker.
    tolerance = 32.0 * sys.float_info.epsilon * max(1.0, abs(cost), abs(old_cost))
    return cost < old_cost - tolerance or (abs(cost - old_cost) <= tolerance and path < old_path)


def _boundaries(first, stop, n_slices, split, norm_se, norm_me):
    if stop - first < n_slices:
        raise ValueError("fewer occupied-extent mT bins than requested slices")
    prefixes = [_prefix(values) for values in (split, norm_se, norm_me)]
    total = math.fsum(split[first:stop])
    if total <= 0.0:
        raise ValueError("signal SE has no counts in split_kstar")
    if any(p[stop] - p[first] <= 0.0 for p in prefixes[1:]):
        raise ValueError("signal SE/ME normalization integral must be positive")
    cumulative = prefixes[0]
    states = {first: (0.0, (first,))}
    for cut in range(1, n_slices):
        next_states = {}
        for end in range(first + cut, stop - (n_slices - cut) + 1):
            fraction = (cumulative[end] - cumulative[first]) / total
            local_cost = (fraction - cut / float(n_slices)) ** 2
            best = None
            for start, (cost, path) in states.items():
                if _valid_interval(prefixes, start, end):
                    candidate = (cost + local_cost, path + (end,))
                    if _better(candidate, best):
                        best = candidate
            if best is not None:
                next_states[end] = best
        states = next_states
    best = None
    for start, (cost, path) in states.items():
        if _valid_interval(prefixes, start, stop):
            candidate = (cost, path + (stop,))
            if _better(candidate, best):
                best = candidate
    if best is None:
        raise ValueError("no {} nonoverlapping mT slices have positive split SE and normalization SE/ME".format(n_slices))
    cost, path = best
    return list(path), cost, total


def _audit_flow(channels, nx, ny, relevant, split, norm, display):
    audit = {"policy": "Never reassign flow counts to regular mT bins",
             "relevant_kstar_indices": relevant, "channels": {},
             "relevant_mt_flow_nonzero": False, "any_flow_nonzero": False}
    for mode in MODES:
        audit["channels"][mode] = {}
        for region in REGIONS:
            channel = channels[mode][region]
            entry = {}
            for label, indices in (("relevant", relevant), ("split", split),
                                   ("normalization", norm), ("display", display),
                                   ("regular_kstar", range(nx))):
                root_x = [x + 1 for x in indices]
                entry[label] = {
                    "mt_underflow": _sum_cells(channel, root_x, [0], True),
                    "mt_overflow": _sum_cells(channel, root_x, [ny + 1], True)}
            entry["kstar_underflow_all_mt"] = _sum_cells(channel, [0], range(ny + 2), True)
            entry["kstar_overflow_all_mt"] = _sum_cells(channel, [nx + 1], range(ny + 2), True)
            entry["mt_flow_all_kstar"] = _sum_cells(channel, range(nx + 2), [0, ny + 1], True)
            if any(_nonzero(v) for v in entry["relevant"].values()):
                audit["relevant_mt_flow_nonzero"] = True
            if any(_nonzero(entry[key]) for key in
                   ("kstar_underflow_all_mt", "kstar_overflow_all_mt", "mt_flow_all_kstar")):
                audit["any_flow_nonzero"] = True
            audit["channels"][mode][region] = entry
    return audit


def _statistics(channels, first, stop, nx, group, rebinned_edges,
                split_indices, norm_indices, display_indices):
    result = {}
    coarse_split = sorted(set(i // group for i in split_indices))
    coarse_norm = sorted(set(i // group for i in norm_indices))
    coarse_display = sorted(set(i // group for i in display_indices))
    for mode in MODES:
        result[mode] = {}
        for region in REGIONS:
            channel = channels[mode][region]
            rebinned = [_sum_cells(channel, range(x, x + group), range(first, stop))
                        for x in range(0, nx, group)]
            counts = [item["count"] for item in rebinned]
            error2 = [item["error2"] for item in rebinned]
            result[mode][region] = {
                "split": _sum_cells(channel, split_indices, range(first, stop)),
                "normalization": _sum_cells(channel, norm_indices, range(first, stop)),
                "all_regular_kstar": _sum_cells(channel, range(nx), range(first, stop)),
                "rebinned": {"edges": list(rebinned_edges), "counts": counts, "error2": error2},
                "zero_stats": {
                    "all": [i for i, count in enumerate(counts) if count == 0.0],
                    "split": [i for i in coarse_split if counts[i] == 0.0],
                    "normalization": [i for i in coarse_norm if counts[i] == 0.0],
                    "display": [i for i in coarse_display if counts[i] == 0.0]}}
    return result


def choose_slices(data, settings, base):
    """Return count-selected or fixed mT intervals, without computing CFs.

    Required JSON is ``channels[SE/ME][region]`` with per-channel ``kstarEdges``,
    ``mtEdges``, x-major ``counts``/``error2`` and ROOT-indexed ``flowCounts``/
    ``flowError2``.  ``first``/``stop`` in the result are zero-based [start,stop).
    Numerical summaries contain ``count`` and ``error2`` (stored Sumw2 sum).

    Required settings: n_slices, split_kstar, split_strategy, mt_extent and
    mass_region. Required base settings: rebin_group, normalization, plot_x_range.
    A relevant mT under/overflow raises :class:`MtFlowError`; its ``audit``
    attribute preserves the flow report.  No minimum-statistics threshold is
    imposed beyond strictly positive split-SE and normalization-SE/ME integrals
    for quantile selection. Fixed mode instead requires ``fixed_mt_edges`` and
    ``mt_extent: configured_range``. It keeps empty slices with count-status
    warnings; unavailable CF normalizations must be masked by the caller.
    """
    strategy = settings["split_strategy"]
    fixed = strategy == "fixed_saved_mt_edges"
    if strategy not in ("nearest_se_signal_quantiles_at_saved_mt_edges", "fixed_saved_mt_edges"):
        raise ValueError("unsupported split_strategy")
    required_extent = "configured_range" if fixed else "occupied_regular_bins_all_channels"
    if settings["mt_extent"] != required_extent:
        raise ValueError("{} requires mt_extent: {}".format(strategy, required_extent))
    if settings["mass_region"] != "signal":
        raise ValueError("split selection currently requires mass_region: signal")
    n_slices, group = settings["n_slices"], base["rebin_group"]
    for value, name in ((n_slices, "n_slices"), (group, "rebin_group")):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError("{} must be a positive integer".format(name))
    source = data["channels"]
    reference = source["SE"]["signal"]
    k_edges = _edges(reference["kstarEdges"], "kstarEdges")
    mt_edges = _edges(reference["mtEdges"], "mtEdges")
    if fixed:
        requested_edges = _edges(settings["fixed_mt_edges"], "fixed_mt_edges")
        if len(requested_edges) != n_slices + 1:
            raise ValueError("fixed_mt_edges must contain n_slices + 1 boundaries")
        fixed_boundaries = [_edge_index(mt_edges, edge, "fixed_mt_edges") for edge in requested_edges]
        if any(b <= a for a, b in zip(fixed_boundaries, fixed_boundaries[1:])):
            raise ValueError("fixed_mt_edges must select distinct increasing saved edges")
    nx, ny = len(k_edges) - 1, len(mt_edges) - 1
    if nx % group:
        raise ValueError("rebin_group must divide the saved regular kstar bin count")
    coarse_edges = k_edges[::group]
    split = _range_indices(k_edges, settings["split_kstar"], "split_kstar")
    norm = _range_indices(k_edges, base["normalization"], "normalization")
    display = _range_indices(k_edges, base["plot_x_range"], "plot_x_range")
    # Selecting partially rebinned cells would change the requested norm/range.
    for bounds, name in ((settings["split_kstar"], "split_kstar"),
                         (base["normalization"], "normalization"),
                         (base["plot_x_range"], "plot_x_range")):
        _range_indices(coarse_edges, bounds, name + " (rebinned)")
    channels = {}
    occupied = set()
    for mode in MODES:
        channels[mode] = {}
        for region in REGIONS:
            raw = source[mode][region]
            name = mode + "/" + region
            if not _same_edges(k_edges, _edges(raw["kstarEdges"], name + " kstarEdges")) or not _same_edges(
                    mt_edges, _edges(raw["mtEdges"], name + " mtEdges")):
                raise ValueError("histogram axes differ for " + name)
            channel = {key: _matrix(raw[key], nx + extra, ny + extra, name + " " + key)
                       for key, extra in (("counts", 0), ("error2", 0),
                                          ("flowCounts", 2), ("flowError2", 2))}
            for x in range(nx):
                for y in range(ny):
                    for regular, flow in (("counts", "flowCounts"), ("error2", "flowError2")):
                        if not _same_number(channel[regular][x][y], channel[flow][x + 1][y + 1]):
                            raise ValueError("regular/flow array interior mismatch for " + name)
                    if channel["counts"][x][y] != 0.0 or channel["error2"][x][y] != 0.0:
                        occupied.add(y)
            channels[mode][region] = channel
    relevant = sorted(set(split + norm + display))
    flow = _audit_flow(channels, nx, ny, relevant, split, norm, display)
    if flow["relevant_mt_flow_nonzero"]:
        raise MtFlowError("nonzero mT flow in split/normalization/display kstar ranges; automatic split rejected", flow)
    if not occupied and not fixed:
        raise ValueError("all regular mT bins are empty")
    split_se = [_sum_cells(channels["SE"]["signal"], split, [y])["count"] for y in range(ny)]
    norm_se = [_sum_cells(channels["SE"]["signal"], norm, [y])["count"] for y in range(ny)]
    norm_me = [_sum_cells(channels["ME"]["signal"], norm, [y])["count"] for y in range(ny)]
    if fixed:
        boundaries, cost = fixed_boundaries, None
        first, stop = boundaries[0], boundaries[-1]
        total = math.fsum(split_se[first:stop])
    else:
        first, stop = min(occupied), max(occupied) + 1
        boundaries, cost, total = _boundaries(first, stop, n_slices, split_se, norm_se, norm_me)
    slices = []
    cumulative = 0.0
    for i, (start, end) in enumerate(zip(boundaries, boundaries[1:])):
        stats = _statistics(channels, start, end, nx, group, coarse_edges, split, norm, display)
        count = stats["SE"]["signal"]["split"]["count"]
        cumulative += count
        slices.append({"index": i, "first": start, "stop": end,
                       "low": mt_edges[start], "high": mt_edges[end],
                       "split_fraction": count / total if total > 0.0 else None,
                       "cumulative_split_fraction": cumulative / total if total > 0.0 else None,
                       "split_signal": {mode: stats[mode]["signal"]["split"] for mode in MODES},
                       "normalization_signal": {mode: stats[mode]["signal"]["normalization"] for mode in MODES},
                       "channels": stats,
                       "signal_cf_count_audit": {
                           "display_bins": len(set(x // group for x in display)),
                           "zero_se_display": len(stats["SE"]["signal"]["zero_stats"]["display"]),
                           "zero_me_display": len(stats["ME"]["signal"]["zero_stats"]["display"]),
                           "positive_me_display": sum(stats["ME"]["signal"]["rebinned"]["counts"][x] > 0.0
                                                      for x in sorted(set(j // group for j in display)))}})
    totals = _statistics(channels, first, stop, nx, group, coarse_edges, split, norm, display)
    closure = {}
    for mode in MODES:
        closure[mode] = {}
        for region in REGIONS:
            full = totals[mode][region]["rebinned"]
            differences = {}
            for key in ("counts", "error2"):
                summed = [math.fsum(s["channels"][mode][region]["rebinned"][key][x] for s in slices)
                          for x in range(len(coarse_edges) - 1)]
                differences[key] = [value - expected for value, expected in zip(summed, full[key])]
            closure[mode][region] = {
                "max_abs_count_difference": max(abs(v) for v in differences["counts"]),
                "max_abs_error2_difference": max(abs(v) for v in differences["error2"]),
                "rebinned_difference": differences,
                "outside_extent": {
                    "below": _sum_cells(channels[mode][region], range(nx), range(first)),
                    "above": _sum_cells(channels[mode][region], range(nx), range(stop, ny))}}
    warnings = []
    if flow["any_flow_nonzero"]:
        warnings.append("Flow exists outside the relevant regular kstar ranges; retained in the audit, not assigned to mT slices.")
    result = {"schema_version": 1, "slices": slices, "boundary_bins": boundaries,
            "mt_edges": [mt_edges[b] for b in boundaries],
            "extent": {"first": first, "stop": stop, "low": mt_edges[first], "high": mt_edges[stop]},
            "quantile_cost": cost, "split_signal_total": total,
            "target_cumulative_fractions": [i / float(n_slices) for i in range(1, n_slices)],
            "selection": {"split_kstar": list(settings["split_kstar"]),
                          "normalization": list(base["normalization"]),
                          "plot_x_range": list(base["plot_x_range"]),
                          "rebin_group": group,
                          "objective": "sum squared cumulative signal-SE quantile deviations",
                          "tie_breaker": "lexicographically lower saved mT edges",
                          "positivity_constraints": ["split signal SE", "normalization signal SE", "normalization signal ME"],
                          "minimum_statistics_threshold": None},
            "totals": totals, "closure": closure, "flow": flow, "warnings": warnings}
    if fixed:
        result["target_cumulative_fractions"] = []
        result["selection"].update({
            "split_strategy": strategy,
            "fixed_mt_edges": requested_edges,
            "objective": "fixed configured saved mT edges; no optimization",
            "tie_breaker": None,
            "positivity_constraints": [],
            "empty_interval_policy": "retain fixed edges, report counts, mask unavailable CF normalizations downstream"})
        result["closure_scope"] = "configured_range"
        result["occupied_extent"] = ({"first": min(occupied), "stop": max(occupied) + 1,
                                      "low": mt_edges[min(occupied)], "high": mt_edges[max(occupied) + 1]}
                                     if occupied else None)
        result["excluded_tails"] = {
            "below": _statistics(channels, 0, first, nx, group, coarse_edges, split, norm, display),
            "above": _statistics(channels, stop, ny, nx, group, coarse_edges, split, norm, display)}
        for entry in slices:
            stats = entry["channels"]
            norm_positive = {mode: stats[mode]["signal"]["normalization"]["count"] > 0.0
                             for mode in MODES}
            sideband_norm_positive = {
                mode: math.fsum(stats[mode][region]["normalization"]["count"]
                                for region in ("leftSB", "rightSB")) > 0.0 for mode in MODES}
            entry["count_status"] = {
                "split_signal_se_positive": entry["split_signal"]["SE"]["count"] > 0.0,
                "normalization_signal_positive": norm_positive,
                "normalization_sideband_positive": sideband_norm_positive,
                "signal_cf_normalization_valid": all(norm_positive.values()),
                "sideband_cf_normalization_valid": all(sideband_norm_positive.values()),
                "empty_all_regular_channels": not any(
                    stats[mode][region]["all_regular_kstar"]["count"] > 0.0
                    for mode in MODES for region in REGIONS)}
            status = entry["count_status"]
            if not status["split_signal_se_positive"] or not status["signal_cf_normalization_valid"] or not status["sideband_cf_normalization_valid"]:
                warnings.append("Fixed slice {} retains its configured edges despite empty split-SE or nonpositive signal/sideband normalization; see count_status and mask unavailable CFs.".format(entry["index"]))
        if any(_nonzero(totals_by_mode[mode][region]["all_regular_kstar"])
               for totals_by_mode in result["excluded_tails"].values()
               for mode in MODES for region in REGIONS):
            warnings.append("Populated regular mT tails lie outside the configured fixed range; excluded counts are audited, not merged into the end slices.")
    return result
