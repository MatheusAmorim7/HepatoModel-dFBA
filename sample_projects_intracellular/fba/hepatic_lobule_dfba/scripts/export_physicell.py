"""Exports a native run (runs/<label>: fluxes/states .f64.gz) to PhysiCell-standard snapshots (MultiCellDS XML + MATLAB v4 .mat) so PhysiCell Studio's Plot tab can open them.
Substrates = the 25 resource fields (mM, tissue concentration after the step = stock/0.008 nL). Cells = the 3000 hosts (positions from hosts.csv) with cell_type = zone (cells.csv) and custom rows = per-host fluxes
(model units; negative = uptake for exchanges) and segment (0 periportal .. 7 central). VISUALISATION ONLY: not an independent result; per-host values inside a voxel are a degenerate split (use voxel/segment aggregates)."""
import gzip, json, re, sys
from pathlib import Path
import numpy as np, pandas as pd, cobra, logging, scipy.io as sio
logging.disable(50)
HERE = Path(__file__).resolve().parent; PRJ = HERE.parent; CFG = PRJ / "config"; T = CFG / "studio_template"
sys.path.insert(0, str(HERE)); import vascular_ref as V
NAMES = ["oxygen", "phosphate", "CO2", "glucose", "lactate", "water", "serine", "glutamate", "glutamine", "ammonia", "aspartate", "alanine", "pyruvate", "bicarbonate", "urea", "asparagine", "tyrosine", "phenylalanine", "acetoacetate", "glycine", "palmitate", "arginine", "ornithine", "citrulline", "glycerol"]
CUSTOM = [("PEPCK_r0123", "r0123"), ("FBP1_r0487", "r0487"), ("G6Pase_r0396", "r0396"), ("glycogenolysis", "GLYCOGEN_SUPPLY"), ("glucose_export", "EX_HC00040_s"), ("lactate_flux", "EX_HC00177_s"), ("glycerol_flux", "EX_HC00114_s"),
          ("O2_flux", "EX_HC00017_s"), ("urea_export", "EX_HC00089_s"), ("palmitate_flux", "EX_HC00226_s"), ("alanine_flux", "EX_HC00048_s"), ("NH3_flux", "EX_HC00024_s")]
RX = cobra.io.read_sbml_model(str(CFG / "postabsorptive_host_v1.xml")).reactions.list_attr("id"); NR = len(RX)
cells = pd.read_csv(CFG / "cells.csv"); zone = cells.type.map({"zone_1": 0, "zone_2": 1, "zone_3": 2}).values
vox = (np.floor((cells.y + 500) / 20).astype(int) * 50 + np.floor((cells.x + 500) / 20).astype(int)).values; seg = V.node_map()[vox] % 8
tmpl = (T / "output00000003.xml").read_text(); mesh = T / "initial_mesh0.mat"
xs = np.tile(-490 + 20 * np.arange(50), 50).astype(float); ys = np.repeat(-490 + 20 * np.arange(50), 50).astype(float)
def labels_xml():
    L = [("ID", 0, 1, "none"), ("position", 1, 3, "microns"), ("total_volume", 4, 1, "cubic microns"), ("cell_type", 5, 1, "none"), ("cycle_model", 6, 1, "none"), ("current_phase", 7, 1, "none"), ("elapsed_time_in_phase", 8, 1, "min"), ("segment_PP0_to_PC7", 9, 1, "none")]
    for i, (n, _) in enumerate(CUSTOM): L.append((n, 10 + i, 1, "model flux units"))
    return "<labels>\n" + "\n".join(f'\t\t\t\t\t\t\t<label index="{i}" size="{s}" units="{u}">{n}</label>' for n, i, s, u in L) + "\n\t\t\t\t\t\t</labels>", 10 + len(CUSTOM)
LAB, NROWS = labels_xml()
VARS = "<variables>\n" + "\n".join(f'\t\t\t\t<variable name="{n}" units="mM" ID="{i}">\n\t\t\t\t\t<physical_parameter_set>\n\t\t\t\t\t\t<conditions />\n\t\t\t\t\t\t<diffusion_coefficient units="micron^2/min">0.000000</diffusion_coefficient>\n\t\t\t\t\t\t<decay_rate units="1/min">0.000000</decay_rate>\n\t\t\t\t\t</physical_parameter_set>\n\t\t\t\t</variable>' for i, n in enumerate(NAMES)) + "\n\t\t\t</variables>"
TYPES = '<cell_types>\n' + "\n".join(f'\t\t\t\t\t\t\t<type ID="{i}" type="{i}">zone_{i + 1}</type>' for i in range(3)) + '\n\t\t\t\t\t\t</cell_types>'
def export(src, label, dt, out, n_snap=None):
    S = json.load(open(src / "summary.json")); n = S["accepted_steps"]; idx = sorted(set([0, n // 20, n // 10, n // 4, n // 2, n - 1])); out.mkdir(parents=True, exist_ok=True)
    import shutil; shutil.copy(mesh, out / "initial_mesh0.mat"); fs = gzip.open(src / "fluxes.f64.gz"); ss = gzip.open(src / "states.f64.gz"); k = 0; cj = [RX.index(r) for _, r in CUSTOM]
    for i in range(n):
        a = np.frombuffer(fs.read(3000 * NR * 8), dtype="<f8").reshape(3000, NR); st = np.frombuffer(ss.read(3 * 2500 * 25 * 8), dtype="<f8").reshape(3, 2500, 25)
        if i not in idx: continue
        name = f"output{k:08d}"; conc = (st[2] / 0.008).T
        sio.savemat(out / f"{name}_microenvironment0.mat", {"multiscale_microenvironment": np.vstack([xs, ys, np.zeros(2500), np.full(2500, 8000.0), conc])}, format="4")
        rows = np.zeros((NROWS, 3000)); rows[0] = np.arange(3000); rows[1] = cells.x.values; rows[2] = cells.y.values; rows[3] = cells.z.values; rows[4] = 2494.0; rows[5] = zone; rows[9] = seg
        for q, j in enumerate(cj): rows[10 + q] = a[:, j]
        sio.savemat(out / f"{name}_cells.mat", {"cells": rows}, format="4")
        x = tmpl.replace("output00000003", name); x = re.sub(r"<current_time units=\"min\">[^<]*<", f'<current_time units="min">{(i + 1) * dt:.6f}<', x)
        x = re.sub(r"<variables>.*?</variables>", lambda m: VARS, x, flags=re.S); x = re.sub(r"<labels>.*?</labels>", lambda m: LAB, x, flags=re.S); x = re.sub(r"<cell_types>.*?</cell_types>", lambda m: TYPES, x, flags=re.S)
        x = re.sub(r"\s*<neighbor_graph.*?</neighbor_graph>|\s*<attached_cells_graph.*?</attached_cells_graph>|\s*<spring_attached_cells_graph.*?</spring_attached_cells_graph>", "", x, flags=re.S)
        (out / f"{name}.xml").write_text(x); k += 1
    (out / "LEIA-ME.txt").write_text(f"Exportado de {src.name} (dt={dt} min). Snapshots nos passos {idx} (tempo = (passo+1)*dt). Apenas visualizacao; ver export_physicell.py.\n"); return idx
if __name__ == "__main__":
    import argparse, xml.etree.ElementTree as ET
    ap = argparse.ArgumentParser(description="Export native runs to PhysiCell-standard output (MultiCellDS XML + .mat) for PhysiCell Studio.")
    ap.add_argument("arms", nargs="*", help="arm names under runs/ (default: every finished non-quick arm)")
    ap.add_argument("--run", help="a run directory (raw output of ./project, e.g. output/raw)"); ap.add_argument("--control", help="control XML of that run (gives the time step)"); ap.add_argument("--out", help="output directory for --run mode")
    a = ap.parse_args()
    if a.run:
        case = ET.parse(a.control).getroot().find("finite_resource_control").find("case"); dt = float(case.get("dt"))
        print(case.get("label"), export(Path(a.run), case.get("label"), dt, Path(a.out or "output"))); sys.exit(0)
    A = json.load(open(CFG / "arms.json"))["arms"]; wanted = a.arms or [x for x in A if (PRJ / "runs" / x / "summary.json").exists() and not x.startswith("quick")]
    for label in wanted: print(label, export(PRJ / "runs" / label, label, A[label]["dt"], PRJ / "studio_output" / label))
