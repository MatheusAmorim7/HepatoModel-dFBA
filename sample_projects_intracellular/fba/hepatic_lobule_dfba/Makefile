# Hepatic lobule dFBA sample project (PhysiCell style).
#   make              build ./hepatic_lobule_dfba against the repository root (core/, BioFVM/, addons/dFBA, custom_modules/)
#   make configs      generate every run configuration (config/arms/*/control.xml) from config/control_base_phaseC.xml.gz
#   make quick        ~5 min smoke test (3 short runs) + comparison with the shipped reference
#   make run ARM=Z3_s3_k15   one arm (about 7-10 min alone; arms: see config/arms.json)
#   make reproduce-zonation | reproduce-emergent | reproduce-atp   all arms of a family, then analysis/figures
#   make analyze figures studio   post-processing (results/), Studio-viewable export (studio_output/)
# Requirements (Ubuntu): g++ with OpenMP, coinor-libclp-dev coinor-libcoinutils-dev libsbml5-dev liblapack-dev libxml2-dev libbz2-dev zlib1g-dev, python3 with numpy scipy pandas cobra matplotlib.
ROOT ?= ../../..
PROGRAM_NAME := hepatic_lobule_dfba
CC := g++
ifdef PHYSICELL_CPP
	CC := $(PHYSICELL_CPP)
endif
PY ?= python3
ARCH := native
CFLAGS := -march=$(ARCH) -O3 -fomit-frame-pointer -mfpmath=both -fopenmp -m64 -std=c++11 -DADDON_PHYSIDFBA
CFLAGS_LINK := $(shell echo $(CFLAGS) | sed -e "s/-fopenmp//g")
INC := -I$(ROOT)/addons/dFBA/src -I$(ROOT) -I.
LIBS := -lClp -lClpSolver -lCoinUtils -llapack -lsbml-static -lxml2 -lbz2 -lz
WRAP := -Wl,--wrap=_ZN6BioFVM11Basic_Agent29simulate_secretion_and_uptakeEPNS_16MicroenvironmentEd -Wl,--wrap=_ZN6BioFVM16Microenvironment24simulate_diffusion_decayEd -Wl,--wrap=_Z32custom_microenvironment_functiond

BioFVM_OBJECTS := BioFVM_vector.o BioFVM_mesh.o BioFVM_microenvironment.o BioFVM_solvers.o BioFVM_matlab.o BioFVM_utilities.o BioFVM_basic_agent.o BioFVM_MultiCellDS.o BioFVM_agent_container.o
CORE_OBJECTS := pugixml.o PhysiCell_phenotype.o PhysiCell_cell_container.o PhysiCell_standard_models.o PhysiCell_cell.o PhysiCell_custom.o PhysiCell_utilities.o PhysiCell_constants.o PhysiCell_basic_signaling.o PhysiCell_signal_behavior.o PhysiCell_rules.o
MODULE_OBJECTS := PhysiCell_SVG.o PhysiCell_pathology.o PhysiCell_MultiCellDS.o PhysiCell_various_outputs.o PhysiCell_pugixml.o PhysiCell_settings.o PhysiCell_geometry.o
DFBA_OBJECTS := dfba_Metabolite.o dfba_Reaction.o dfba_Solution.o dfba_Model.o dfba_intracellular.o
ALL_OBJECTS := $(BioFVM_OBJECTS) $(CORE_OBJECTS) $(MODULE_OBJECTS) custom.o $(DFBA_OBJECTS)
OBJ_PATHS := $(addprefix $(ROOT)/,$(ALL_OBJECTS))
SRCS := main.cpp $(wildcard src/*.h src/*/*.h src/*/*.cpp)

all: $(PROGRAM_NAME)

# object files are built by the repository root Makefile (same flags as the rest of the project); order-only so a plain "make" does not relink needlessly
.PHONY: objects
objects:
	$(MAKE) -C $(ROOT) $(ALL_OBJECTS)

$(PROGRAM_NAME): $(SRCS) | objects
	$(CC) $(CFLAGS_LINK) $(CFLAGS) $(INC) -o $@ $(OBJ_PATHS) main.cpp -L$(ROOT)/addons/dFBA/ext/libsbml/lib $(LIBS) $(WRAP)
	@echo "Executable name is $(PROGRAM_NAME)"

configs:
	$(PY) scripts/make_configs.py

# usage: make run ARM=Z3_s3_k15
run: $(PROGRAM_NAME)
	@test -n "$(ARM)" || (echo "usage: make run ARM=<arm> (see config/arms.json)"; exit 1)
	./scripts/run_arm.sh $(ARM)

quick: $(PROGRAM_NAME) configs
	./scripts/run_arm.sh quick_Z0 && ./scripts/run_arm.sh quick_E3 && ./scripts/run_arm.sh quick_A2
	$(PY) scripts/verify_quick.py

reproduce-zonation: $(PROGRAM_NAME) configs
	./scripts/run_family.sh Z0_uniform Z3_s3_k15 Z10_s10_k15 Z3c_caps_only
	$(PY) scripts/analyze.py Z0_uniform Z3_s3_k15 Z10_s10_k15 Z3c_caps_only && $(PY) scripts/score.py && $(PY) scripts/make_figures.py

reproduce-emergent: $(PROGRAM_NAME) configs
	./scripts/run_family.sh E0_uniform_k10 E3_s3_k10 E3_s3_k15 E0_uniform_k15
	$(PY) scripts/analyze.py E0_uniform_k10 E3_s3_k10 E3_s3_k15 E0_uniform_k15 && $(PY) scripts/score.py && $(PY) scripts/make_figures.py

reproduce-atp: $(PROGRAM_NAME) configs
	$(PY) scripts/scan_atpm.py
	./scripts/run_family.sh A1_uniform_pal124 A2_uniform_pal449 A3_s3_pal449
	$(PY) scripts/analyze.py A1_uniform_pal124 A2_uniform_pal449 A3_s3_pal449 && $(PY) scripts/score.py && $(PY) scripts/make_figures.py

analyze:
	$(PY) scripts/analyze.py $(ARMS) && $(PY) scripts/score.py
figures:
	$(PY) scripts/make_figures.py
studio:
	$(PY) scripts/export_physicell.py

clean:
	rm -f $(PROGRAM_NAME)
	rm -rf runs config/arms
data-cleanup:
	rm -rf runs studio_output
# Install system packages and Python modules (Ubuntu/Debian/WSL/Colab; needs root or sudo for apt).
deps:
	./scripts/setup_ubuntu.sh

# Standalone bundle: only what is needed to build and run (no git history, no notebooks). Unpack and run "make quick" inside the project.
BUNDLE ?= lobule-dfba-standalone.tar.gz
bundle:
	tar -C $(ROOT) -czf $(BUNDLE) --transform 's,^,hepatomodel-dfba/,' \
	  --exclude='*.o' --exclude='sample_projects_intracellular/fba/hepatic_lobule_dfba/runs' --exclude='sample_projects_intracellular/fba/hepatic_lobule_dfba/config/arms' \
	  --exclude='sample_projects_intracellular/fba/hepatic_lobule_dfba/studio_output' --exclude='sample_projects_intracellular/fba/hepatic_lobule_dfba/hepatic_lobule_dfba' --exclude='__pycache__' \
	  Makefile VERSION.txt licenses BioFVM core modules addons/dFBA/src addons/dFBA/ext custom_modules/custom.cpp custom_modules/custom.h \
	  sample_projects_intracellular/fba/hepatic_lobule_dfba
	@echo "wrote $(BUNDLE)"; ls -la $(BUNDLE)

.PHONY: all deps bundle configs run quick reproduce-zonation reproduce-emergent reproduce-atp analyze figures studio clean data-cleanup
