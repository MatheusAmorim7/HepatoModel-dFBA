# legacy/ — previous single-driver lobule model

Code that is **not used by the current model** (`sample_projects_intracellular/fba/hepatic_lobule_dfba`) and was moved here from the repository root to keep the active code small:

* `main.cpp` — the standard PhysiCell main loop used with this project's sinusoidal-supply lobule (it called `custom_microenvironment_function` every diffusion step and logged memory use).
* `custom_modules/custom.cpp`, `custom.h` — the full custom module: sinusoidal oxygen/glucose supply (`sinusoidal_distributed_supply`, `dynamic_oxygen_supply`, `dynamic_glucose_supply`), zone initialisation, peripheral sources, `setup_tissue`, coloring function and the cell-definition defaults. The slim version that the current model needs stays in `custom_modules/`.

The legacy driver still reads `./config/PhysiCell_settings.xml` and `./config/cells.csv` (unchanged). To build it, from the repository root: `make -f legacy/Makefile`, then `./legacy_project`.

* `custom_modules/exchange_audit.h` — helper used by the legacy `main.cpp` to write the exchange inventory.

Removed earlier (see git history, commit "Remove legacy harness code ..."): the shared-pool event-control harness and the finite-oxygen feed. They were only used by experiments whose sealed results live in `notebooks/`.

Status: kept for reference; **not verified against the current results**. The current model does not execute any of this code.
