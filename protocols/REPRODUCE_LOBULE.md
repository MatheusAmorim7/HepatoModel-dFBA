# Reproduce the lobular simulation

This update contains the current `config`, `custom_modules`, dFBA sources, `main.cpp` and `Makefile`. Keep the other PhysiCell core files from the existing repository. No simulation output or Python notebook is required to run this lobular configuration.

## Ubuntu / WSL

```bash
sudo apt update
sudo apt install build-essential pkg-config coinor-libclp-dev libsbml5-dev libxml2-dev libbz2-dev zlib1g-dev liblapack-dev libblas-dev
bash protocols/build_ubuntu.sh
./project config/PhysiCell_settings.xml
```

For a short numerical check after building, run `python3 protocols/smoke_lobule.py`. It uses a temporary output directory and checks optimal solver status, glucose/O2 uptake and memory logging without overwriting the normal simulation output.

The build script uses system libraries rather than requiring precompiled files in `addons/dFBA/ext`. It compiles in a temporary checkout so older object files cannot silently be reused. It places the new `project` executable in the repository root.

The current configuration loads `config/cells.csv` (3,000 hepatocytes) and creates the supply/drainage agents in the custom code. Alternative CSVs are supplied for other agent counts. Change the cells filename deliberately before changing scale. The simulation writes its own output files to the configured output folder.

`main.cpp` includes Linux memory reporting (`/proc/self/status`) and writes `memory_usage.csv`. The code uses the configuration parameters for distributed O2/glucose replenishment and includes prescribed lactate/CO2 exchange rules after FBA. These fields represent effective extracellular concentrations. The current objective maximises glucose transport (`r1032`), not a validated hepatocyte growth objective.

The tumour serine/HPRT assays are a separate Python package, not part of this executable. Libraries and third-party models retain their existing terms and citations.
