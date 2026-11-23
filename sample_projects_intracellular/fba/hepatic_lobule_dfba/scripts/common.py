"""Shared constants/paths for the post-processing scripts."""
import json
from pathlib import Path
import numpy as np, pandas as pd
HERE = Path(__file__).resolve().parent; PRJ = HERE.parent; CFG = PRJ / "config"; RUNS = PRJ / "runs"; RES = PRJ / "results"; ANA = RES / "analysis"; FIG = RES / "figures"
for d in (ANA, FIG): d.mkdir(parents=True, exist_ok=True)
import sys; sys.path.insert(0, str(HERE)); import vascular_ref as V
NAMES = ["oxygen", "phosphate", "CO2", "glucose", "lactate", "water", "serine", "glutamate", "glutamine", "ammonia", "aspartate", "alanine", "pyruvate", "bicarbonate", "urea", "asparagine", "tyrosine", "phenylalanine", "acetoacetate", "glycine", "palmitate", "arginine", "ornithine", "citrulline", "glycerol"]
EX = {"oxygen": "EX_HC00017_s", "phosphate": "EX_HC00019_s", "CO2": "EX_HC00021_s", "glucose": "EX_HC00040_s", "lactate": "EX_HC00177_s", "water": "EX_HC00011_s", "serine": "EX_HC00068_s", "glutamate": "EX_HC00034_s", "glutamine": "EX_HC00067_s", "ammonia": "EX_HC00024_s",
      "aspartate": "EX_HC00055_s", "alanine": "EX_HC00048_s", "pyruvate": "EX_HC00032_s", "bicarbonate": "EX_HC00253_s", "urea": "EX_HC00089_s", "asparagine": "EX_HC00148_s", "tyrosine": "EX_HC00085_s", "phenylalanine": "EX_HC00081_s", "acetoacetate": "EX_HC00159_s",
      "glycine": "EX_HC00045_s", "palmitate": "EX_HC00226_s", "arginine": "EX_HC00065_s", "ornithine": "EX_HC00079_s", "citrulline": "EX_HC00277_s", "glycerol": "EX_HC00114_s"}
KEY = ["r0123", "r0487", "r0396", "GLYCOGEN_SUPPLY", "EX_HC00040_s", "EX_HC00177_s", "EX_HC00114_s", "EX_HC00017_s", "EX_HC00089_s", "EX_HC00226_s", "EX_HC00048_s", "EX_HC00024_s"]
UNIT = 0.0031194          # model flux units per (umol/min of the whole liver): 6e7 cells... see docs/MODEL_NOTES.md
def reaction_ids():
    import cobra, logging; logging.disable(50)
    return cobra.io.read_sbml_model(str(CFG / "postabsorptive_host_v1.xml")).reactions.list_attr("id")
def cells_info():
    cells = pd.read_csv(CFG / "cells.csv"); vox = (np.floor((cells.y + 500) / 20).astype(int) * 50 + np.floor((cells.x + 500) / 20).astype(int)).values
    return cells, vox, V.node_map()[vox] % 8, V.node_map() % 8
def arms(): return json.load(open(CFG / "arms.json"))["arms"]
def inlet_mM():
    import xml.etree.ElementTree as ET
    fc = ET.parse(CFG / "arms/Z0_uniform/control.xml").getroot().find("finite_resource_control"); return {e.get("field"): float(e.get("mM")) for e in fc.find("vascular").findall("inlet")}
