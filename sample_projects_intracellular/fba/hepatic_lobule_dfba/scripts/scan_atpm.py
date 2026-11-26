"""Static scan of ATP maintenance (CULTURE_ATPM) in the post-absorptive reference scenario (package v1, f=1): feasibility, L1 solution, O2 uptake, fuels. New files only.
Target arithmetic in TARGETS (liver O2 uptake) is documented in PROTOCOL_ATP.md."""
import csv, json, logging
from pathlib import Path
import numpy as np, cobra
from cobra.util.array import create_stoichiometric_matrix
from scipy import sparse
from scipy.optimize import linprog
logging.disable(50)
HERE = Path(__file__).resolve().parent; PKG = HERE.parent / "config"; RESDIR = HERE.parent / "results"
M = cobra.io.read_sbml_model(str(PKG / "postabsorptive_host_v1.xml")); ids = M.reactions.list_attr("id"); n = len(ids)
for r in csv.DictReader((PKG / "scenario_reference_f1.csv").open()):
    M.reactions.get_by_id(r["reaction_id"]).bounds = (float(r["lower"]), float(r["upper"]))
S = sparse.csr_matrix(create_stoichiometric_matrix(M)); UNIT = 0.0031194                  # units per umol/min of whole liver
jATP = ids.index("CULTURE_ATPM"); lo0 = np.array([r.lower_bound for r in M.reactions]); hi0 = np.array([r.upper_bound for r in M.reactions])
def l1(atpm, free=False):
    lo, hi = lo0.copy(), hi0.copy(); lo[jATP] = hi[jATP] = atpm
    # variables: v (n), a (n) >= |v|
    A_eq = sparse.hstack([S, sparse.csr_matrix((S.shape[0], n))]); b_eq = np.zeros(S.shape[0])
    I = sparse.identity(n, format="csr"); A_ub = sparse.vstack([sparse.hstack([I, -I]), sparse.hstack([-I, -I])]); b_ub = np.zeros(2 * n)
    c = np.concatenate([np.zeros(n), np.ones(n)]); bounds = [(lo[j], hi[j]) for j in range(n)] + [(0, None)] * n
    return linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method="highs")
def maxatp():
    lo, hi = lo0.copy(), hi0.copy(); lo[jATP] = 0; hi[jATP] = 1000; c = np.zeros(n); c[jATP] = -1
    r = linprog(c, A_eq=S, b_eq=np.zeros(S.shape[0]), bounds=list(zip(lo, hi)), method="highs"); return -r.fun if r.status == 0 else None
out = {"max_feasible_ATPM_in_scenario": maxatp(), "scan": {}}
show = {"O2": "EX_HC00017_s", "palmitate": "EX_HC00226_s", "lactate": "EX_HC00177_s", "glycerol": "EX_HC00114_s", "alanine": "EX_HC00048_s", "glucose": "EX_HC00040_s", "urea": "EX_HC00089_s", "glycogen": "GLYCOGEN_SUPPLY", "ammonia": "EX_HC00024_s", "CO2": "EX_HC00021_s", "glutamine": "EX_HC00067_s"}
for a in (5, 8, 10, 12, 15, 20, 25, 30, 40):
    r = l1(a)
    if r.status != 0: out["scan"][a] = {"status": "infeasible"}; continue
    v = r.x[:n]; d = {k: round(float(v[ids.index(x)]), 4) for k, x in show.items()}
    d["O2_mL_min_liver"] = round(-v[ids.index("EX_HC00017_s")] / UNIT * 22.4 / 1000, 2); d["status"] = "optimal"; d["L1"] = round(float(r.fun), 3)
    d["at_ceiling"] = [i for i in ("EX_HC00226_s", "EX_HC00114_s", "EX_HC00177_s", "EX_HC00017_s") if abs(v[ids.index(i)] - lo0[ids.index(i)]) < 1e-6 and lo0[ids.index(i)] < 0]
    out["scan"][a] = d
json.dump(out, open(RESDIR / "SCAN_ATPM.json", "w"), indent=1)
print("max feasible ATPM", out["max_feasible_ATPM_in_scenario"])
for a, d in out["scan"].items(): print(a, d)
