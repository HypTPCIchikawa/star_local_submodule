import math

def compare(expected,actual,policy):
    if set(expected)!=set(actual):raise RuntimeError('Histogram names mismatch')
    waived=[];allowed=set(policy['histograms'])
    for name,(schema,values) in expected.items():
        other_schema,other=actual[name]
        if schema!=other_schema:raise RuntimeError('Histogram schema mismatch: '+name)
        if any(not math.isfinite(x) for x in values+other):raise RuntimeError('Non-finite histogram: '+name)
        if values[0]!=other[0]:raise RuntimeError('Histogram Entries mismatch: '+name)
        different=any(not math.isclose(a,b,rel_tol=policy['otherHistogramsRelativeTolerance'],abs_tol=policy['otherHistogramsAbsoluteTolerance']) for a,b in zip(values[1:],other[1:]))
        if different:
            if name not in allowed:raise RuntimeError('Unapproved histogram sum difference: '+name)
            waived.append({'histogram':name,'schema':schema,'expected':values,'actual':other,'relativePercent':[(b-a)/a*100 if a else None for a,b in zip(values,other)]})
    return waived
