# Model notes

## Network
* Source: HepatoNet1 (Gille et al., Mol Syst Biol 2010; SBML L2V3). The simulations use a compact curated host (179 reactions) extended to **package v1 (225 reactions, 217 metabolites)**: urea-cycle/N exchanges, a palmitate β-oxidation block, glycerol utilisation and a glycogen-mobilisation reaction that is **ours** (`GLYCOGEN_SUPPLY`, 11 glucose residues per cycle). See `config/PACKAGE.json`, `config/NATIVE_REQUIREMENTS.md`.
* Reference scenario (post-absorptive): glucose export 2.14 units (9.8 µmol/kg/min), urea 0.19, gluconeogenic share ≥ 40%, lactate uptake ceiling 6, palmitate ≤ 0.4, glycerol ≤ 0.23, ATP maintenance 5 (`host_f1_atpm5.xml`) or 21.6 after calibration (`host_atpm21.6.xml`, O₂ uptake ceiling −20).

## Units and scale
* 1 model flux unit = 0.12 pmol/h per host. 3000 hosts represent 1.6×10¹¹ hepatocytes (107×10⁶ cells/g × 21.4 g/kg × 70 kg), so 1 µmol/min of the whole liver = 0.0031194 units per host.
* Voxel = 20×20×20 µm (0.008 nL; mM = pmol/nL). Domain plasma flow = 18.75 nL/min per L/min of liver plasma flow (0.86 L/min → 16.125 nL/min).
* Blood: 6 sectors × 8 segments serial nodes (segment 0 periportal … 7 central), blood volume 2 nL, linear O₂ buffer β = 107.45 (haemoglobin; a linearisation, not Hb kinetics), tissue–blood conductance 3.504 nL/min per voxel.

## What “atomic” means
Each metabolic step solves one joint LP per voxel (L1 objective, shared finite stock). If any LP is infeasible the whole step is rejected (exit code 42, `rejection.json`): no repair, no clipping, no ATP relaxation, no dropped hosts.

## Time step
The step must keep (consumption per step) < (stock available to the voxel): dt ≲ 0.1 min at ATPM 5 with the local blood pool, 0.02 min at ATPM 21.6.

## Reproducibility check
`make quick` compares the first 20 steps of three runs with the published runs (relative difference ≤ 1e-9; we obtained 0.0 from a fresh clone with only the files in this branch).

## Not in this project
The tumour focus, serine export capacity (earlier static analyses), cell cycle/death/motion, independent validation. See `docs/` for the status of each.
