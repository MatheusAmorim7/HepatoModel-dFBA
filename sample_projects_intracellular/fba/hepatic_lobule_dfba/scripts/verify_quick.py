#!/usr/bin/env python3
"""Compares the first steps of runs/quick_* with results/reference/quick_reference.json (metrics computed from the original published runs). Tolerance: relative 1e-9."""
import gzip, json, sys
import numpy as np
from common import *
RX = reaction_ids(); NR = len(RX); col = np.array([RX.index(EX[n]) for n in NAMES]); kj = [RX.index(k) for k in KEY]; cells, vox, seg_h, seg_v = cells_info(); unit = 0.12 / 60
def metrics(gz, n):
    f = gzip.open(gz); upt = np.zeros((n, 25)); pepck = np.zeros((n, 8))
    for i in range(n):
        a = np.frombuffer(f.read(3000 * NR * 8), dtype="<f8").reshape(3000, NR); upt[i] = -a[:, col].sum(axis=0) * unit
        pepck[i] = [a[seg_h == s][:, RX.index("r0123")].mean() for s in range(8)]
    return dict(uptake_last=upt[-1].tolist(), uptake_sum=upt.sum(axis=0).tolist(), PEPCK_last_by_segment=pepck[-1].tolist())
if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--make-reference":      # usage: verify_quick.py --make-reference <quick_label>=<path to original fluxes.f64.gz> ...
        ref = {}
        for kv in sys.argv[2:]: lab, path = kv.split("="); ref[lab] = metrics(path, 20)
        json.dump(ref, open(RES / "reference/quick_reference.json", "w"), indent=1); print("reference written", list(ref)); sys.exit(0)
    ref = json.load(open(RES / "reference/quick_reference.json")); ok = True
    for lab, r in ref.items():
        f = RUNS / lab / "fluxes.f64.gz"
        if not f.exists(): print("MISSING", lab); ok = False; continue
        m = metrics(f, 20); err = max(float(np.max(np.abs(np.array(m[k]) - np.array(r[k])) / (np.abs(np.array(r[k])) + 1e-12))) for k in r)
        print(f"{lab}: max relative difference vs reference = {err:.2e}", "OK" if err < 1e-9 else "FAIL"); ok &= err < 1e-9
    sys.exit(0 if ok else 1)
