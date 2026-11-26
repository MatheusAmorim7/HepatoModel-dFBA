#!/usr/bin/env python3
"""Scores the registered predictions of the three families from results/analysis -> results/SCORES.json (see docs/PROTOCOLS for the predictions). Runs only the parts whose analyses exist."""
import json
import numpy as np
from common import *
cells, vox, seg_h, seg_v = cells_info(); cnt = np.bincount(seg_h, minlength=8); IN = inlet_mM(); A = arms(); Q = 16.125; BETA = 107.45239901452385
def load(a):
    f = ANA / f"ANALYSIS_{a}.npz"; return np.load(f) if f.exists() else None
D = {a: load(a) for a in A}; D = {a: d for a, d in D.items() if d is not None}; S = {}
def tot(d, f, sl=slice(-50, None)): return float(d["upt"][sl, NAMES.index(f)].mean())
conv = lambda p: p / (3000 * 0.002) / UNIT                      # domain pmol/min -> umol/min of whole liver
def per_seg(d, key, i=-1): return d["segflux"][i][KEY.index(key)]
ref = D.get("Z0_uniform")
if ref is not None:
    z = {}
    for a in ("Z0_uniform", "Z3_s3_k15", "Z10_s10_k15"):
        if a in D:
            d = D[a]; z[a] = dict(PEPCK_by_segment=[round(float(x), 4) for x in per_seg(d, "r0123")], PEPCK_ratio_seg0_seg7=round(float(per_seg(d, "r0123")[0] / per_seg(d, "r0123")[7]), 3), monotone_decreasing=bool(np.all(np.diff(per_seg(d, "r0123")) < 0)),
                                  lactate_blood_PC_vs_PP_pct=round(float(100 * (d["ves"][-1][NAMES.index("lactate")][7] / d["ves"][-1][NAMES.index("lactate")][0] - 1)), 2),
                                  organ_totals_late_pmol_min={f: round(tot(d, f), 4) for f in ("lactate", "glycerol", "oxygen")}, glucose_export_mean=round(float((d["hostlast"][:, KEY.index("EX_HC00040_s")]).mean()), 4))
    S["zonation_prescribed"] = z
for a in A:
    if a.startswith("E") and a in D:
        d = D[a]; s = A[a]["s"]; k = A[a]["kappa"]
        m_seg = np.array([ (1 + (1 / s - 1) * ((j + .5) / 8)) for j in range(8)]) / np.mean(1 + (1 / s - 1) * ((seg_h + .5) / 8)); J = {"r0123": 1.482, "r0487": 0.856}
        r = dict(glucose_export_mean_per_host=round(float((per_seg(d, "EX_HC00040_s") * cnt).sum() / cnt.sum()), 4), PEPCK_ratio_seg0_seg7=round(float(per_seg(d, "r0123")[0] / per_seg(d, "r0123")[7]), 3), capacity_ratio_seg0_seg7=round(float(m_seg[0] / m_seg[7]), 3),
                 PEPCK_over_cap=[round(float(x), 3) for x in per_seg(d, "r0123") / (k * m_seg * J["r0123"])], FBP1_over_cap=[round(float(x), 3) for x in per_seg(d, "r0487") / (k * m_seg * J["r0487"])],
                 lactate_extraction=round(tot(d, "lactate") / (Q * IN["lactate"]), 4))
        if ref is not None: r["vs_Z0_organ_totals_pct"] = {f: round(100 * (tot(d, f) / tot(ref, f) - 1), 2) for f in ("lactate", "glycerol", "oxygen")}
        S.setdefault("emergent", {})[a] = r
for a in A:
    if a.startswith("A") and a in D:
        d = D[a]; o2 = tot(d, "oxygen"); v = d["ves"][-1][NAMES.index("oxygen")] / IN["oxygen"]
        r = dict(O2_mL_min=round(conv(o2) * 22.4 / 1000, 2), O2_extraction=round(o2 / (Q * IN["oxygen"] * (1 + BETA)), 3), palmitate_umol_min=round(conv(tot(d, "palmitate")), 1), lactate_extraction=round(tot(d, "lactate") / (Q * IN["lactate"]), 3),
                 blood_O2_content_rel_inlet_by_segment=[round(float(x), 3) for x in v], O2_content_drop_seg0_to_7_pct=round(100 * float(1 - v[7] / v[0]), 1), O2_per_host_by_segment=[round(float(x), 3) for x in -per_seg(d, "EX_HC00017_s")])
        S.setdefault("atp_calibration", {})[a] = r
for a, d in A.items():
    f = ANA / f"ANALYSIS_{a}.json"
    if f.exists(): j = json.load(open(f)); S.setdefault("status", {})[a] = dict(status=j["status"], steps=j.get("accepted_steps"), max_ledger_error_pmol=j.get("max_abs_ledger_error_pmol"), failed_voxel_segment=j.get("failed_voxel_segment"))
json.dump(S, open(RES / "SCORES.json", "w"), indent=1); print(json.dumps(S.get("status", {}), indent=1))
