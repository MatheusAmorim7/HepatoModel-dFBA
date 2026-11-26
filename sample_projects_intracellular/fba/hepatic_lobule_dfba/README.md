# Hepatic lobule dFBA — zonation, vascular transport and ATP calibration (PhysiCell sample project)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/MatheusAmorim7/HepatoModel-dFBA/blob/lobule-dfba-reproducible/sample_projects_intracellular/fba/hepatic_lobule_dfba/notebooks/run_lobule_dfba.ipynb)

**Fastest way to run it:** click the badge → in Colab choose *Runtime → Run all* (about 15 min for build + quick test, nothing to install). Other ways (Docker, no-git bundle, Ubuntu) are in *Easiest ways to run* below.


A PhysiCell + dFBA model of a liver lobule: 3000 hepatocyte agents on a 50×50 voxel domain, each with a 225-reaction metabolic network derived from **HepatoNet1** (post-absorptive reference scenario), serial portal→central blood transport with finite stocks, enzyme zonation, and an ATP-maintenance calibration against hepatic O₂ consumption. It is the platform for the thesis question *how much serine can the liver supply to liver metastases of pancreatic cancer* (the tumour is **not** in this project yet).

> **Status: exploratory.** Nothing here is a validated model. Zonation profiles, ATP maintenance and several parameters are declared hypotheses; independent review is pending. Read `docs/` before citing any number.

## What is in the box
| Path | Content |
|---|---|
| `main.cpp`, `src/` | The simulation program (native PhysiCell driver with an atomic per-voxel joint LP selector, serial vascular transport, finite stocks, LP cache). |
| `config/` | Host networks (`host_f1_atpm5.xml`, `host_atpm21.6.xml`, package v1 = 225 reactions + scenario), cell positions (`cells.csv`), base configuration (`control_base_phaseC.xml.gz`). |
| `scripts/` | `make_configs.py` (generates all arms), `run_arm.sh`, `analyze.py`, `score.py`, `make_figures.py`, `export_physicell.py` (PhysiCell-standard output for PhysiCell Studio), `scan_atpm.py`, `verify_quick.py`. |
| `results/reference/` | Figures, scores and the ATPM scan **from the original runs** (what we report). Your runs write to `results/` (git-ignored). |
| `docs/` | Protocols with the predictions registered *before* each run, and result reports (Portuguese). |

## Easiest ways to run (pick one)
| You have | Do this | Status |
|---|---|---|
| **A browser** (no install) | Open the notebook in **Google Colab**: [`notebooks/run_lobule_dfba.ipynb`](https://colab.research.google.com/github/MatheusAmorim7/HepatoModel-dFBA/blob/lobule-dfba-reproducible/sample_projects_intracellular/fba/hepatic_lobule_dfba/notebooks/run_lobule_dfba.ipynb) and run all cells. It clones, installs, builds, runs the quick test and shows the figures. | notebook tested headless on Ubuntu; **not yet tested on Colab itself** |
| **Jupyter on Ubuntu / WSL2** | Same notebook (`jupyter lab`, open `notebooks/run_lobule_dfba.ipynb`). | tested headless on Ubuntu |
| **Windows or macOS** | **Docker** (see below). Native Windows is not supported; on Windows you can also install WSL2 + Ubuntu and follow the Ubuntu route. | Dockerfile **not tested** (no Docker where it was written) |
| **No git** | Download `lobule-dfba-standalone.tar.gz` (2 MB; `make bundle` creates it), `tar xzf`, `cd hepatomodel-dfba/sample_projects_intracellular/fba/hepatic_lobule_dfba`, `make deps && make quick`. | tested (build + quick) |
| **Ubuntu + git** | `git clone --depth 1 --branch lobule-dfba-reproducible https://github.com/MatheusAmorim7/HepatoModel-dFBA.git` (**always use `--depth 1`**: the full history is ~1 GB), then `cd` into the project and `make deps && make quick`. | tested |

**Docker:** `docker build -t lobule-dfba .` (from the repository root), then `docker run --rm -it -p 8888:8888 lobule-dfba` and open the link it prints; the notebook is in the file list. Headless: `docker run --rm -it lobule-dfba make quick`.

> If something fails, `make quick` is the first thing to run: it must print `OK` three times. Please report the **last 30 lines of its output** and your OS/Python versions.

## Requirements
Ubuntu/WSL (tested with g++ 13), Python ≥ 3.10.
```bash
sudo apt install g++ make coinor-libclp-dev coinor-libcoinutils-dev libsbml5-dev \
                 liblapack-dev libxml2-dev libbz2-dev zlib1g-dev
pip install numpy scipy pandas cobra matplotlib        # use a virtualenv; pass PY=... to make if needed
```
About 4 GB RAM (each run ≈1.2 GB; three can run in parallel), ~1.5 GB disk per arm of raw output.

## Quick start (PhysiCell style)
From the repository root of this project:
```bash
cd sample_projects_intracellular/fba/hepatic_lobule_dfba
make                  # builds ./hepatic_lobule_dfba (first build compiles PhysiCell core, a few minutes)
make quick            # 3 short runs (~5-10 min) and a bit-level check against the reference
```
`make quick` must end with `OK` for the three runs (relative difference ≤ 1e-9 against the published runs). If it does, your build reproduces ours.

## Reproduce the results
Each arm takes ~25 min (a family runs three arms in parallel):
```bash
make reproduce-zonation     # Z0 uniform, Z3 (s=3), Z10 (s=10), Z3c (caps only; expected to be REJECTED at step 0)
make reproduce-emergent     # E0/E3 arms: glucose output emerges from enzyme capacity
make reproduce-atp          # static ATPM scan + A1/A2/A3 (ATPM 21.6, liver O2 ≈ 43 mL/min)
make run ARM=Z3_s3_k15      # a single arm (list: config/arms.json)
```
Outputs: `runs/<arm>/` (raw binary: `fluxes.f64.gz`, `states.f64.gz`, `vessel.f64.gz`, `ledger.csv`, `summary.json`), `results/analysis/ANALYSIS_<arm>.{json,npz}`, `results/SCORES.json`, `results/figures/*.png`. Compare with `results/reference/`.

## Look at it in PhysiCell Studio
```bash
make studio            # writes studio_output/<arm>/output0000000N.xml (+ .mat) in standard PhysiCell format
```
In Studio open the **Plot** tab and point the output directory at one of the folders (see `docs/STUDIO_PT.md`). Do not open `config/PhysiCell_settings.xml` of the repository root in Studio (it crashes on the `dfba` intracellular block). The export is visualisation only.

## How a run works (short)
1. `scripts/make_configs.py` writes `config/arms/<arm>/control.xml`: a PhysiCell settings file plus a `finite_resource_control` block (25 resource fields, blood inlet concentrations, per-host capacity profiles, one `case` with time step, steps and flow).
2. `./hepatic_lobule_dfba <control.xml> <label> <outdir> replay warm gzip` advances the metabolism in steps of `dt` min: per voxel one joint LP (L1 objective, shared finite stock); a step is accepted only if every LP is feasible (**atomic transaction, no repair, no clipping**; exit code 42 = rejected, see `rejection.json`). Blood transport is serial and implicit, with a mass-balance certificate per field.
3. Units: one model flux unit = 0.12 pmol/h per host; 1 µmol/min of the whole liver = 0.0031194 units (107e6 hepatocytes/g × 21.4 g/kg × 70 kg); `docs/` has the derivations.

## Arms
| Family | Arms | Question |
|---|---|---|
| Z (prescribed zonation) | `Z0_uniform`, `Z3_s3_k15`, `Z10_s10_k15`, `Z3c_caps_only` | Does zonating glucose output and PEPCK/FBP1/G6Pase capacity change organ totals or gradients? |
| E (emergent) | `E0_uniform_k10`, `E3_s3_k10`, `E0_uniform_k15`, `E3_s3_k15` | If only enzyme capacity is given, does the flux zonation emerge? |
| A (ATP calibration) | `A1_uniform_pal124`, `A2_uniform_pal449`, `A3_s3_pal449` | With ATPM calibrated to liver O₂ uptake, does the O₂ gradient come out right? |

## Headline findings (details and failures in `docs/`)
* Organ totals are identical with and without zonation: zonation redistributes where, not how much.
* The flux gradient follows the enzyme-capacity gradient; supply modifies it by ≲10%.
* One calibrated parameter (ATPM 5 → 21.6; liver O₂ 16.7 → 43 mL/min) brings the sinusoidal O₂ content drop from ~4% to ~20%, close to literature PO₂ ranges that were **not** used in the fit; O₂ extraction (0.26 vs ~0.37) and palmitate uptake remain low.
* Registered predictions that failed are listed in each report; none was tuned away.

## Provenance and limits
HepatoNet1 (Gille et al., Mol Syst Biol 2010) is the source of the network; the compact host and the package-v1 additions are documented in `config/PACKAGE.json` and `config/NATIVE_REQUIREMENTS.md`. Literature values used as targets come from web-search summaries and must be checked against the originals. Raw runs are not versioned (regenerate them with the targets above).

## Citation / license
Built on PhysiCell (BSD-3, see repository root licenses). If you use this project, please cite PhysiCell (Ghaffarizadeh et al., PLoS Comput Biol 2018) and HepatoNet1.
