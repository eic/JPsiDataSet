#!/usr/bin/env python3
"""Event-by-event comparison of Lorentz-invariant quantities between lAger samples,
so that afterburned (boosted/rotated) events can be compared with generator output.

usage: compare_invariants.py LABEL=file.hepmc LABEL=file.hepmc [...]
The first sample is the one every other sample is compared to.
"""
import importlib.util
import os
import sys
import numpy as np
from scipy import stats

spec = importlib.util.spec_from_file_location(
    "compare", os.path.join(os.path.dirname(os.path.abspath(__file__)), "compare.py"))
compare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(compare)
mdot, read_hepmc = compare.mdot, compare.read_hepmc


def gamma_vectors(path):
    """Summed FSR photon 4-vector per event (zero if no FSR)."""
    out, cur, jpsi_id, gams = [], None, None, None
    with open(path) as f:
        for line in f:
            if line.startswith("E "):
                if cur is not None:
                    out.append(cur)
                cur, jpsi_id = np.zeros(4), None
            elif line.startswith("P "):
                t = line.split()
                if int(t[3]) == 443:
                    jpsi_id = int(t[1])
                elif int(t[3]) == 22 and int(t[2]) == jpsi_id:
                    cur += np.array([float(x) for x in t[4:8]])
        out.append(cur)
    return np.array(out)


def invariants(path):
    d, _ = read_hepmc(path)
    k, kp, P, Pp, V, lm, lp = (d[x] for x in ("k", "kp", "P", "Pp", "jpsi", "lm", "lp"))
    q = k - kp
    MV = np.sqrt(mdot(V, V))
    g = gamma_vectors(path)
    ll = lm + lp
    return {
        "Q2": -mdot(q, q),
        "y": mdot(P, q) / mdot(P, k),
        "W": np.sqrt(mdot(P + q, P + q)),
        "-t": -mdot(P - Pp, P - Pp),
        "M(J/psi)": MV,
        "N FSR": d["ngam"].astype(float),
        "E_gamma (J/psi frame)": mdot(V, g) / MV,
        "M(e+e-)": np.sqrt(mdot(ll, ll)),
        "decay angle (q.l-/q.V)": mdot(q, lm) / mdot(q, V),
    }


PRODUCTION = ("Q2", "y", "W", "-t", "M(J/psi)")
DECAY = ("N FSR", "E_gamma (J/psi frame)", "M(e+e-)", "decay angle (q.l-/q.V)")


def same(a, b, rtol=1e-6, atol=1e-9):
    return np.isclose(a, b, rtol=rtol, atol=atol)


def main():
    samples = [a.split("=", 1) for a in sys.argv[1:]]
    inv = {}
    for label, path in samples:
        inv[label] = invariants(path)
        print(f"read {label}: {len(inv[label]['Q2'])} events")
    ref = samples[0][0]
    for label, _ in samples[1:]:
        a, b = inv[ref], inv[label]
        n = min(len(a["Q2"]), len(b["Q2"]))
        print(f"\n== {ref} vs {label}: first {n} events, event by event (rel. tol. 1e-6)")
        prod = np.all([same(a[x][:n], b[x][:n]) for x in PRODUCTION], axis=0)
        dec = np.all([same(a[x][:n], b[x][:n]) for x in DECAY], axis=0)
        for x in PRODUCTION + DECAY:
            m = same(a[x][:n], b[x][:n])
            print(f"  {x:26s} identical in {100 * m.mean():6.2f}% of events")
        print(f"  {'all production invariants':26s} identical in {100 * prod.mean():6.2f}% of events")
        print(f"  {'all decay invariants':26s} identical in {100 * dec.mean():6.2f}% of events")
        print("  distributions (KS p-value, mean A / mean B):")
        for x in PRODUCTION + DECAY:
            p = stats.ks_2samp(a[x][:n], b[x][:n]).pvalue
            print(f"    {x:26s} p={p:6.3f}   {np.mean(a[x][:n]):.5g} / {np.mean(b[x][:n]):.5g}")


if __name__ == "__main__":
    main()
