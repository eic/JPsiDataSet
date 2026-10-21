"""Per-event FSR photon multiplicity correlation between lAger samples.

usage: fsr_corr.py REFERENCE.hepmc SEED1.hepmc SEED2.hepmc
"""
import sys, numpy as np
def ngam(path):
    out=[]; n=None
    for line in open(path):
        if line.startswith("E "):
            if n is not None: out.append(n)
            n=int(line.split()[3])
    out.append(n); return np.array(out)-9
a,b,c=(ngam(p) for p in sys.argv[1:4])
m=min(len(b),len(c))
print("per-event FSR multiplicity sequence, seed1 vs seed2: identical in %.1f%% of events, corr=%.3f" % (100*np.mean(b[:m]==c[:m]), np.corrcoef(b[:m],c[:m])[0,1]))
m=min(len(a),len(b))
print("reference vs seed1: identical in %.1f%% of events, corr=%.3f" % (100*np.mean(a[:m]==b[:m]), np.corrcoef(a[:m],b[:m])[0,1]))
