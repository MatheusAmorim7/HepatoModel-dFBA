#!/usr/bin/env python3
"""Stream runs/<arm> (fluxes/states/vessel .f64.gz) -> results/analysis/ANALYSIS_<arm>.{json,npz}. Usage: python3 scripts/analyze.py arm [arm ...]"""
import csv, gzip, json, sys
import numpy as np
from common import *
RX = reaction_ids(); NR = len(RX); kj = [RX.index(k) for k in KEY]; col = np.array([RX.index(EX[n]) for n in NAMES]); unit = 0.12 / 60
cells, vox, seg_h, seg_v = cells_info(); A = arms()
def analyze(label):
    RUN = RUNS / label; dt = A[label]["dt"]; out = dict(label=label, arm=A[label])
    if (RUN / "rejection.json").exists():
        out.update(status="REJECTED_ATOMICALLY", rejection=json.load(open(RUN / "rejection.json"))); v = out["rejection"]["failed_voxel"]; out["failed_voxel_segment"] = int(seg_v[v])
        json.dump(out, open(ANA / f"ANALYSIS_{label}.json", "w"), indent=1); return out
    S = json.load(open(RUN / "summary.json")); n = S["accepted_steps"]; out.update(status="completed", accepted_steps=n)
    led = list(csv.DictReader((RUN / "ledger.csv").open())); out["max_abs_ledger_error_pmol"] = max(abs(float(r["ledger_error_pmol"])) for r in led)
    segflux = np.zeros((n, len(KEY), 8)); upt = np.zeros((n, 25)); ves = np.zeros((n, 25, 8)); keep = sorted(set([0, 49, 149, n - 1]) & set(range(n))); maps = {}; hostlast = None
    fs = gzip.open(RUN / "fluxes.f64.gz"); ss = gzip.open(RUN / "states.f64.gz"); vs = gzip.open(RUN / "vessel.f64.gz")
    for i in range(n):
        a = np.frombuffer(fs.read(3000 * NR * 8), dtype="<f8").reshape(3000, NR); st = np.frombuffer(ss.read(3 * 2500 * 25 * 8), dtype="<f8").reshape(3, 2500, 25); vz = np.frombuffer(vs.read(25 * 48 * 8), dtype="<f8").reshape(25, 48)
        upt[i] = -a[:, col].sum(axis=0) * unit
        for s in range(8): segflux[i, :, s] = a[seg_h == s][:, kj].mean(axis=0)
        ves[i] = vz.reshape(25, 6, 8).mean(axis=1)
        if i in keep: maps[i] = st[2] / 0.008
        if i == n - 1: hostlast = a[:, kj].copy()
    np.savez_compressed(ANA / f"ANALYSIS_{label}.npz", segflux=segflux, upt=upt, ves=ves, hostlast=hostlast, seg_h=seg_h, seg_v=seg_v, keep=np.array(keep), maps=np.array([maps[k] for k in keep]), KEY=np.array(KEY), NAMES=np.array(NAMES), dt=dt)
    json.dump(out, open(ANA / f"ANALYSIS_{label}.json", "w"), indent=1); return out
if __name__ == "__main__":
    for lab in sys.argv[1:]:
        o = analyze(lab); print(lab, o["status"], o.get("accepted_steps"), o.get("max_abs_ledger_error_pmol"))
