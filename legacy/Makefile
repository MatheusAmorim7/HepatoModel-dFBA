# Legacy single-driver lobule model (sinusoidal oxygen/glucose supply, PhysiCell main loop).
# Kept for reference; it is NOT used by the current model (sample_projects_intracellular/fba/hepatic_lobule_dfba).
# Build and run from the repository root:
#   make -f legacy/Makefile            # builds ./legacy_project (needs the shared objects: make -C . first builds them)
#   ./legacy_project                   # reads ./config/PhysiCell_settings.xml and ./config/cells.csv
PROGRAM_NAME := legacy_project
CC := g++
ifdef PHYSICELL_CPP
	CC := $(PHYSICELL_CPP)
endif
CFLAGS := -march=native -O3 -fomit-frame-pointer -mfpmath=both -fopenmp -m64 -std=c++11 -DADDON_PHYSIDFBA
INC := -I./addons/dFBA/src -I. -I./custom_modules
LIBS := -lClp -lClpSolver -lCoinUtils -llapack -lsbml-static -lxml2 -lbz2 -lz
BioFVM_OBJECTS := BioFVM_vector.o BioFVM_mesh.o BioFVM_microenvironment.o BioFVM_solvers.o BioFVM_matlab.o BioFVM_utilities.o BioFVM_basic_agent.o BioFVM_MultiCellDS.o BioFVM_agent_container.o
CORE_OBJECTS := pugixml.o PhysiCell_phenotype.o PhysiCell_cell_container.o PhysiCell_standard_models.o PhysiCell_cell.o PhysiCell_custom.o PhysiCell_utilities.o PhysiCell_constants.o PhysiCell_basic_signaling.o PhysiCell_signal_behavior.o PhysiCell_rules.o
MODULE_OBJECTS := PhysiCell_SVG.o PhysiCell_pathology.o PhysiCell_MultiCellDS.o PhysiCell_various_outputs.o PhysiCell_pugixml.o PhysiCell_settings.o PhysiCell_geometry.o
DFBA_OBJECTS := dfba_Metabolite.o dfba_Reaction.o dfba_Solution.o dfba_Model.o dfba_intracellular.o
OBJECTS := $(BioFVM_OBJECTS) $(CORE_OBJECTS) $(MODULE_OBJECTS) $(DFBA_OBJECTS)

$(PROGRAM_NAME): legacy/main.cpp legacy/custom_modules/custom.cpp legacy/custom_modules/custom.h
	$(MAKE) $(OBJECTS)
	$(CC) $(CFLAGS) $(INC) -o $@ $(OBJECTS) legacy/main.cpp legacy/custom_modules/custom.cpp $(LIBS)
