#!/usr/bin/env python3
"""Compare lAger HepMC3 ASCII samples (e p -> e' J/psi p', J/psi -> e+ e- (gamma)).

usage: compare.py OUTDIR LABEL=file.hepmc [LABEL=file.hepmc ...]
The first sample is the reference.
"""
import sys
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def mdot(a, b):
    return a[..., 3] * b[..., 3] - np.sum(a[..., :3] * b[..., :3], axis=-1)


def read_hepmc(path):
    """Return dict of per-event 4-vectors (px,py,pz,E) and the raw particle rows."""
    out = {k: [] for k in ("k", "kp", "P", "Pp", "jpsi", "lm", "lp", "ngam", "egam")}
    raw = []

    def flush(ev):
        if ev is None:
            return
        parts = ev
        beam_e = next(p for p in parts.values() if p["st"] == 4 and p["pid"] == 11)
        beam_p = next(p for p in parts.values() if p["st"] == 4 and p["pid"] == 2212)
        jpsi = next(p for p in parts.values() if p["pid"] == 443)
        scat = next(p for p in parts.values() if p["pid"] == 11 and p["st"] == 1 and p["par"] == beam_e["id"])
        rec = next(p for p in parts.values() if p["pid"] == 2212 and p["st"] == 1)
        lm = next(p for p in parts.values() if p["pid"] == 11 and p["par"] == jpsi["id"])
        lp = next(p for p in parts.values() if p["pid"] == -11 and p["par"] == jpsi["id"])
        gams = [p for p in parts.values() if p["pid"] == 22 and p["par"] == jpsi["id"]]
        for key, p in (("k", beam_e), ("kp", scat), ("P", beam_p), ("Pp", rec), ("jpsi", jpsi), ("lm", lm), ("lp", lp)):
            out[key].append(p["v"])
        out["ngam"].append(len(gams))
        out["egam"].append(sum(g["v"][3] for g in gams))
        raw.append(np.concatenate([p["v"] for p in sorted(parts.values(), key=lambda p: p["id"])]))

    ev = None
    vtx_parent = {}
    with open(path) as f:
        for line in f:
            c = line[:2]
            if c == "E ":
                flush(ev)
                ev, vtx_parent = {}, {}
            elif c == "V ":
                tok = line.split()
                ins = tok[3].strip("[]").split(",")
                vtx_parent[int(tok[1])] = int(ins[0])
            elif c == "P ":
                tok = line.split()
                pid_, par = int(tok[1]), int(tok[2])
                if par < 0:  # produced at a vertex: use its first incoming particle
                    par = vtx_parent.get(par, par)
                ev[pid_] = dict(id=pid_, par=par, pid=int(tok[3]), st=int(tok[9]),
                                v=np.array([float(x) for x in tok[4:8]]))
        flush(ev)
    d = {k: np.array(v) for k, v in out.items()}
    return d, raw


def observables(d):
    k, kp, P, Pp, V = d["k"], d["kp"], d["P"], d["Pp"], d["jpsi"]
    q = k - kp
    o = {}
    o["Q2 [GeV2]"] = -mdot(q, q)
    o["y"] = mdot(P, q) / mdot(P, k)
    o["W [GeV]"] = np.sqrt(mdot(P + q, P + q))
    dP = P - Pp
    o["-t [GeV2]"] = -mdot(dP, dP)
    o["x_B"] = o["Q2 [GeV2]"] / (2 * mdot(P, q))

    def pt(v): return np.hypot(v[:, 0], v[:, 1])
    def eta(v):
        p = np.linalg.norm(v[:, :3], axis=1)
        return 0.5 * np.log((p + v[:, 2]) / (p - v[:, 2]))
    def mom(v): return np.linalg.norm(v[:, :3], axis=1)

    o["J/psi pT [GeV]"] = pt(V)
    o["J/psi eta"] = eta(V)
    o["J/psi p [GeV]"] = mom(V)
    o["scat. e- eta"] = eta(kp)
    o["scat. e- p [GeV]"] = mom(kp)
    o["recoil p pT [GeV]"] = pt(Pp)
    o["recoil p eta"] = eta(Pp)
    o["decay e- eta"] = eta(d["lm"])
    o["decay e+ eta"] = eta(d["lp"])
    o["decay e- p [GeV]"] = mom(d["lm"])
    ll = d["lm"] + d["lp"]
    o["M(e+e-) [GeV]"] = np.sqrt(mdot(ll, ll))
    o["FSR photon E [GeV]"] = d["egam"]
    o["N FSR photons"] = d["ngam"].astype(float)
    # four-momentum conservation check (beam in = everything final)
    fin = kp + Pp + d["lm"] + d["lp"]
    o["|sum p_in - p_out| (no FSR) [GeV]"] = np.linalg.norm((k + P - fin)[:, :3], axis=1)
    return o


LOGX = {"Q2 [GeV2]", "y", "-t [GeV2]", "x_B", "FSR photon E [GeV]"}


def main():
    outdir = sys.argv[1]
    samples = [a.split("=", 1) for a in sys.argv[2:]]
    data = {}
    for label, path in samples:
        d, raw = read_hepmc(path)
        data[label] = (observables(d), raw)
        print(f"read {label}: {len(raw)} events from {path}")

    labels = [s[0] for s in samples]
    ref = labels[0]
    ref_obs = data[ref][0]

    # exact-match check of every other sample against the reference
    print("\n== event-by-event comparison vs", ref)
    for lab in labels[1:]:
        a, b = data[ref][1], data[lab][1]
        n = min(len(a), len(b))
        same_shape = sum(1 for i in range(n) if a[i].shape == b[i].shape)
        ident = sum(1 for i in range(n) if a[i].shape == b[i].shape and np.allclose(a[i], b[i], rtol=1e-9, atol=1e-12))
        first_diff = next((i for i in range(n) if not (a[i].shape == b[i].shape and np.allclose(a[i], b[i], rtol=1e-9, atol=1e-12))), None)
        print(f"  {lab}: {ident}/{n} events identical (rel 1e-9); same topology {same_shape}/{n}; first differing event: {first_diff}")

    print("\n== distribution comparison (KS two-sample p-value; means)")
    hdr = f"  {'observable':38s}" + "".join(f"{l:>22s}" for l in labels) + "".join(f"{'KS p '+l:>18s}" for l in labels[1:])
    print(hdr)
    rows = []
    for name in ref_obs:
        line = f"  {name:38s}"
        for lab in labels:
            x = data[lab][0][name]
            line += f"{np.mean(x):>12.5g} ±{np.std(x)/np.sqrt(len(x)):<9.2g}"
        ps = []
        for lab in labels[1:]:
            p = stats.ks_2samp(ref_obs[name], data[lab][0][name]).pvalue
            ps.append(p)
            line += f"{p:>18.3f}"
        rows.append((name, ps))
        print(line)

    with PdfPages(f"{outdir}/comparison.pdf") as pdf:
        names = list(ref_obs)
        for i0 in range(0, len(names), 6):
            fig, axes = plt.subplots(2, 3, figsize=(15, 8))
            for ax, name in zip(axes.flat, names[i0:i0 + 6]):
                allx = np.concatenate([data[l][0][name] for l in labels])
                lo, hi = np.percentile(allx, [0.1, 99.9])
                if name in LOGX and lo > 0:
                    bins = np.logspace(np.log10(lo), np.log10(hi), 50); ax.set_xscale("log")
                elif name == "N FSR photons":
                    bins = np.arange(-0.5, allx.max() + 1.5)
                else:
                    bins = np.linspace(lo, hi, 50)
                for lab in labels:
                    x = data[lab][0][name]
                    h, e = np.histogram(x, bins=bins)
                    w = 1.0 / len(x)
                    c = 0.5 * (e[1:] + e[:-1])
                    if lab == ref:
                        ax.hist(x, bins=bins, weights=np.full(len(x), w), histtype="stepfilled", alpha=0.3, label=lab)
                    else:
                        ax.errorbar(c, h * w, yerr=np.sqrt(h) * w, fmt=".", ms=3, label=lab)
                ax.set_xlabel(name); ax.set_ylabel("fraction / bin")
                if name in ("-t [GeV2]", "FSR photon E [GeV]", "M(e+e-) [GeV]", "|sum p_in - p_out| (no FSR) [GeV]"):
                    ax.set_yscale("log")
                ax.legend(fontsize=7)
            fig.tight_layout(); pdf.savefig(fig); fig.savefig(f"{outdir}/comparison_p{i0 // 6 + 1}.png", dpi=90); plt.close(fig)
    print(f"\nplots: {outdir}/comparison.pdf")


if __name__ == "__main__":
    main()
