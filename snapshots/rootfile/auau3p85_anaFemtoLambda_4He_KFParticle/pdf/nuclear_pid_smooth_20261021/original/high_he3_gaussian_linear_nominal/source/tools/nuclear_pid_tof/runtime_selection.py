"""Read the actual legacy tables and current opt-in TPC selector, without edits."""
import bisect
import hashlib
import math
from pathlib import Path
import re

import yaml

SPECIES = ('d', 't', 'he3', 'he4')
TOL = 1.e-12  # Matches the existing C++ provider's numerical equality tolerance.


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def parse_legacy_tables(header):
    contents = Path(header).read_text()
    tables = {}
    for key, name, count in (('mean', 'kMean', 4), ('sigma', 'kSigma', 4),
                             ('domain', 'kPRange', 2)):
        match = re.search(r'\b'+name+r'\s*\[[^\]]+\]\s*=\s*\{(.*?)\};', contents, re.S)
        if not match:
            raise ValueError('Cannot read actual C++ table: '+name)
        rows = [[float(value.strip()) for value in row.split(',')]
                for row in re.findall(r'\{([^{}]+)\}', match.group(1))]
        if len(rows) != 4 or any(len(row) != count or not all(map(math.isfinite, row)) for row in rows):
            raise ValueError('Malformed actual C++ table: '+name)
        tables[key] = rows
    if any(row[3] != 0 for row in tables['sigma']):
        raise ValueError('The current sigma table is no longer quadratic.')
    if any(lo <= 0 or hi <= lo for lo, hi in tables['domain']):
        raise ValueError('Bad actual C++ rigidity domain.')
    tables['source'] = str(Path(header).resolve())
    tables['sourceSha256'] = digest(header)
    return tables


def polynomial(coefficients, x):
    value = 0.
    for coefficient in reversed(coefficients):
        value = value*x+coefficient
    return value


def legacy_parameters(tables, species, rigidity, ignore_domain=False):
    if species not in range(4) or not math.isfinite(rigidity) or rigidity <= 0:
        return 0., 0.
    lo, hi = tables['domain'][species]
    if not ignore_domain and not lo <= rigidity <= hi:
        return 0., 0.
    x = math.log10(rigidity)
    return polynomial(tables['mean'][species], x), polynomial(tables['sigma'][species], x)


def load_current_selection(project, mainconf):
    project = Path(project).resolve()
    mainconf = Path(mainconf).resolve(strict=True)
    main = yaml.safe_load(mainconf.read_text())
    nuclear_path = (project/'config'/main['nuclearid']).resolve(strict=True)
    cuts = yaml.safe_load(nuclear_path.read_text())
    if cuts['m2_selection']:
        raise ValueError('The ALL dE/dx histogram alone cannot reproduce active joint m2 selection.')
    if cuts.get('he4HighRigidityTof', {}).get('enabled', False):
        raise ValueError('The ALL histogram cannot reproduce the active He4 high-R TOF selector.')
    if not cuts['requireBestSpecies']:
        raise ValueError('Only the actual best-species current profile is supported.')
    for key in ('maxNSigmaNuclear', 'maxPOverQ'):
        if not math.isfinite(cuts[key]) or cuts[key] <= 0:
            raise ValueError('Missing or invalid current TPC cut: '+key)
    tables = parse_legacy_tables(project/'StMaker/common/NuclearIdDeDxVsMom.h')
    calibration = None
    files = {str(mainconf): digest(mainconf), str(nuclear_path): digest(nuclear_path),
             tables['source']: tables['sourceSha256']}
    if 'dedxCalibration' in cuts:
        enabled = cuts['dedxCalibration']['enabled']
        if not isinstance(enabled, bool):
            raise ValueError('Explicit calibration switch must be boolean.')
        if enabled:
            calibration_path = (project/'config'/cuts['dedxCalibration']['calibrationFile']).resolve(strict=True)
            calibration = yaml.safe_load(calibration_path.read_text())
            if (calibration['schemaVersion'], calibration['interpolation'], calibration['edgePolicy'],
                    calibration['outsideRangePolicy'], calibration['rigidityUnit'], calibration['dedxUnit']) != (
                    1, 'linear', 'clampWithinDomain', 'legacy', 'GeV/c', 'keV/cm'):
                raise ValueError('Unsupported current runtime calibration schema.')
            if (calibration['source']['collision'], calibration['source']['year'],
                    calibration['source']['production']) != ('auau3p85', 2021, 'P24iy'):
                raise ValueError('Current calibration has the wrong dataset.')
            if set(calibration['species']) != {'he3', 'he4'}:
                raise ValueError('Current runtime provider supports the explicit helium pair only.')
            for name in ('he3', 'he4'):
                item = calibration['species'][name]
                positions = [k['rigidity'] for k in item['knots']]
                if len(positions) < 2 or any(b <= a for a, b in zip(positions, positions[1:])):
                    raise ValueError('Unordered current calibration.')
                if any(not all(math.isfinite(k[p]) and k[p] > 0 for p in ('rigidity', 'mean', 'sigma'))
                       for k in item['knots']):
                    raise ValueError('Invalid current calibration knot.')
            files[str(calibration_path)] = digest(calibration_path)
    return dict(tables=tables, calibration=calibration, cuts=cuts, protectedFiles=files,
                runtimeMainconf=str(mainconf), nuclearConfig=str(nuclear_path))


def current_parameters(runtime, species, rigidity):
    legacy = legacy_parameters(runtime['tables'], species, rigidity)
    calibration = runtime['calibration']
    if species < 2 or calibration is None or not math.isfinite(rigidity) or rigidity <= 0:
        return (*legacy, 'legacy')
    item = calibration['species'][SPECIES[species]]
    lo, hi = item['domain']
    if rigidity < lo-TOL or rigidity > hi+TOL:
        return (*legacy, 'legacy_outside_domain')
    knots = item['knots']
    positions = [k['rigidity'] for k in knots]
    i = bisect.bisect_left(positions, rigidity)
    if i < len(knots) and abs(positions[i]-rigidity) <= TOL:
        k = knots[i]
        return k['mean'], k['sigma'], 'calibrated_knot'
    if i > 0 and abs(positions[i-1]-rigidity) <= TOL:
        k = knots[i-1]
        return k['mean'], k['sigma'], 'calibrated_knot'
    if i == 0 or i == len(knots):
        k = knots[0] if i == 0 else knots[-1]
        return k['mean'], k['sigma'], 'calibrated_edge_clamp'
    left, right = knots[i-1], knots[i]
    gap = right['rigidity']-left['rigidity']
    if gap > calibration['maxInterpolationGap']+TOL:
        half = item['supportBinHalfWidth']
        if rigidity <= left['rigidity']+half+TOL:
            return left['mean'], left['sigma'], 'calibrated_gap_edge_clamp'
        if rigidity >= right['rigidity']-half-TOL:
            return right['mean'], right['sigma'], 'calibrated_gap_edge_clamp'
        return (*legacy, 'legacy_unsupported_gap')
    fraction = (rigidity-left['rigidity'])/gap
    return tuple(left[k]+fraction*(right[k]-left[k]) for k in ('mean', 'sigma'))+('calibrated_interpolation',)


def current_active(runtime, species, rigidity):
    mean, sigma, _ = current_parameters(runtime, species, rigidity)
    return bool(rigidity > 0 and rigidity < runtime['cuts']['maxPOverQ'] and mean > 0 and sigma > 0)


def gl3_nodes(lo, hi):
    if not (math.isfinite(lo) and math.isfinite(hi) and hi > lo > 0):
        raise ValueError('Bad positive rigidity bin.')
    mid, half = (hi+lo)/2., (hi-lo)/2.
    offset = math.sqrt(3./5.)*half
    return ((mid-offset, 5./18.), (mid, 8./18.), (mid+offset, 5./18.))
