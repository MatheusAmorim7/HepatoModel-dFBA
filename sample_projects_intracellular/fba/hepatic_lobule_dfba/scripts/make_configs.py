#!/usr/bin/env python3
"""Generates the run configuration (PhysiCell-style control XML) of every arm from config/control_base_phaseC.xml.gz.
Families: Z (prescribed zonation), E (emergent glucose output), A (ATP calibration); 'quick' = short smoke test. Usage: python3 scripts/make_configs.py [--steps N] (writes config/arms/<arm>/control.xml and config/arms.json)."""
import argparse, gzip, hashlib, json, io
from pathlib import Path
import numpy as np, pandas as pd
import xml.etree.ElementTree as ET
HERE = Path(__file__).resolve().parent; PRJ = HERE.parent; CFG = PRJ / "config"
import sys; sys.path.insert(0, str(HERE)); import vascular_ref as V
ap = argparse.ArgumentParser(); ap.add_argument("--steps-quick", type=int, default=20); ap.add_argument("--default", action="store_true", help="only write config/PhysiCell_settings.xml (the calibrated lobule, arm A3) for the PhysiCell-style ./project run"); args = ap.parse_args()
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
cells = pd.read_csv(CFG / "cells.csv"); assert len(cells) == 3000
vox = (np.floor((cells.y + 500) / 20).astype(int) * 50 + np.floor((cells.x + 500) / 20).astype(int)).values; seg = V.node_map()[vox] % 8; x = (seg + .5) / 8
JREF = {"r0123": 1.482, "r0487": 0.856, "r0396": 2.0232727272727}; NOM = 1000.0; LAMBDA = 200.0
PAL449 = (124 + 60 * 18 / 16 + 150 * 18 / 16 + 79 * 18 / 16) / 1000
NAMES = ["oxygen", "phosphate", "CO2", "glucose", "lactate", "water", "serine", "glutamate", "glutamine", "ammonia", "aspartate", "alanine", "pyruvate", "bicarbonate", "urea", "asparagine", "tyrosine", "phenylalanine", "acetoacetate", "glycine", "palmitate", "arginine", "ornithine", "citrulline", "glycerol"]
BASE = gzip.open(CFG / "control_base_phaseC.xml.gz").read()
def load(): return ET.parse(io.BytesIO(BASE)).getroot()
def case_template():
    fc = load().find("finite_resource_control"); return [c for c in fc.findall("case") if c.get("label") == "cL3_q086_dt01"][0]
def weights(s): E = 1 + (1 / s - 1) * x; return E / E.mean()
def common(root, label, host, case_over, atpm=5):
    fc = root.find("finite_resource_control"); fc.set("model", str(CFG / host)); fc.set("atpm_expected", repr(atpm)); out = PRJ / "runs"; fc.set("output", str(out)); root.find("save/folder").text = str(out)
    for tag in ("initial_conditions", "cell_rules"):
        for e in root.iter("folder"):
            if e.text and ("PhysiCelldFBA/config" in e.text or "PhysiCellLobule" in e.text): e.text = str(CFG)
    for c in list(fc.findall("case")): fc.remove(c)
    c = case_template(); c.set("label", label)
    for k, v in case_over.items(): c.set(k, v)
    fc.append(c); return fc
def write(root, label):
    d = CFG / "arms" / label; d.mkdir(parents=True, exist_ok=True); ET.ElementTree(root).write(d / "control.xml", encoding="utf-8", xml_declaration=True); return sha(d / "control.xml")
def profiles(fc, kappa, m, task=False, glucose_free=False):
    cm = fc.find("capacity_maps"); [cm.remove(p) for p in list(cm)]
    for c in range(3000):
        at = {"cell_id": str(c), "r0054": "1.0", "r0353": "1.0"}
        for r in ("r0396", "r0123", "r0487"): at[r] = repr(1.0 if kappa is None else float(kappa * m[c] * JREF[r] / NOM))
        if task: at["task_glucose"] = repr(float(m[c]))
        if glucose_free: at["glucose_ub"] = "10.0"; at["glycogen_cap"] = repr(float(m[c]))
        ET.SubElement(cm, "profile", at)
arms = {}
for name, (s, kappa, task) in {"Z0_uniform": (1, None, False), "Z3_s3_k15": (3, 1.5, True), "Z10_s10_k15": (10, 1.5, True), "Z3c_caps_only": (3, 1.5, False)}.items():
    if args.default: continue
    root = load(); fc = common(root, name, "host_f1_atpm5.xml", {}); profiles(fc, kappa, weights(s), task=task)
    arms[name] = dict(family="Z", atpm=5, dt=0.1, s=s, kappa=kappa, sha256=write(root, name))
for name, (s, kappa) in {"E0_uniform_k10": (1, 1.0), "E3_s3_k10": (3, 1.0), "E0_uniform_k15": (1, 1.5), "E3_s3_k15": (3, 1.5)}.items():
    if args.default: continue
    root = load(); fc = common(root, name, "host_f1_atpm5.xml", {"glucose_reward": repr(LAMBDA)}); profiles(fc, kappa, weights(s), glucose_free=True)
    arms[name] = dict(family="E", atpm=5, dt=0.1, s=s, kappa=kappa, sha256=write(root, name))
jp = NAMES.index("palmitate")
for name, (s, kappa, task, pal) in {"A1_uniform_pal124": (1, None, False, 0.124), "A2_uniform_pal449": (1, None, False, PAL449), "A3_s3_pal449": (3, 1.5, True, PAL449)}.items():
    if args.default and name != "A3_s3_pal449": continue
    root = load(); fc = common(root, name, "host_atpm21.6.xml", {"dt": "0.02", "steps": "300", "dt_mechanics": "0.1", "dt_phenotype": "2.0"}, atpm=21.6); profiles(fc, kappa, weights(s), task=task)
    [e for e in fc.find("vascular").findall("inlet") if e.get("field") == "palmitate"][0].set("mM", repr(pal))
    for row in fc.find("initial_stocks"): v = row.text.split(); v[jp] = repr(pal * 0.008); row.text = " ".join(v)
    if args.default and name == "A3_s3_pal449":
        fc.set("output", str(PRJ / "output/raw")); root.find("save/folder").text = str(PRJ / "output")
        ET.ElementTree(root).write(CFG / "PhysiCell_settings.xml", encoding="utf-8", xml_declaration=True); print("wrote", CFG / "PhysiCell_settings.xml"); sys.exit(0)
    arms[name] = dict(family="A", atpm=21.6, dt=0.02, s=s, kappa=kappa, palmitate_mM=pal, sha256=write(root, name))
for name, base, atp in (("quick_Z0", "Z0_uniform", 5), ("quick_E3", "E3_s3_k10", 5), ("quick_A2", "A2_uniform_pal449", 21.6)):
    root = ET.parse(CFG / "arms" / base / "control.xml").getroot(); fc = root.find("finite_resource_control"); c = fc.find("case"); c.set("label", name); c.set("steps", str(args.steps_quick))
    arms[name] = dict(family="quick", atpm=atp, dt=float(c.get("dt")), base=base, steps=args.steps_quick, sha256=write(root, name))
json.dump(dict(inlet_note="inlet concentrations are stored in each control.xml (vascular/inlet)", arms=arms), open(CFG / "arms.json", "w"), indent=1); print("arms:", list(arms))
