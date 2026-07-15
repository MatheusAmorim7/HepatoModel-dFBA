VERSION := $(shell grep . VERSION.txt | cut -f1 -d:)
PROGRAM_NAME := project

CC := g++

ifdef PHYSICELL_CPP
	CC := $(PHYSICELL_CPP)
endif

ifndef STATIC_OPENMP
	STATIC_OPENMP = -fopenmp
endif

ARCH := native

# ---------------------------------------------------------------
# dFBA dependencies — paths relative to project root
# ---------------------------------------------------------------
DFBA_SRC     := ./addons/dFBA/src
DFBA_EXT     := ./addons/dFBA/ext
DFBA_INC     := -I$(DFBA_SRC) -I$(DFBA_EXT)/libsbml/include -I$(DFBA_EXT)/coin-or/include
DFBA_LD      := -L./addons/dFBA/ext/libsbml/lib
DFBA_RPATH   := -Wl,-rpath,./addons/dFBA/ext/libsbml/lib
DFBA_LIBS    := -lClp -lClpSolver -lCoinUtils -llapack -lsbml-static -lxml2 -lbz2 -lz

# ---------------------------------------------------------------
# Compiler flags
# ---------------------------------------------------------------
CFLAGS := -march=$(ARCH) -O3 -fomit-frame-pointer -mfpmath=both -fopenmp -m64 -std=c++11 -DADDON_PHYSIDFBA

ifeq ($(OS),Windows_NT)
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Darwin)
		UNAME_P := $(shell uname -p)
		var := $(shell which $(CC) | xargs file)
		ifeq ($(lastword $(var)),arm64)
			CFLAGS := -march=$(ARCH) -O3 -fomit-frame-pointer -fopenmp -m64 -std=c++11 -DADDON_PHYSIDFBA
		endif
	endif
endif

CFLAGS_LINK     := $(shell echo $(CFLAGS) | sed -e "s/-fopenmp//g")
COMPILE_COMMAND := $(CC) $(CFLAGS) $(DFBA_INC)
LINK_COMMAND    := $(CC) $(CFLAGS_LINK) $(DFBA_INC)

# ---------------------------------------------------------------
# PhysiCell core objects
# ---------------------------------------------------------------
BioFVM_OBJECTS := BioFVM_vector.o BioFVM_mesh.o BioFVM_microenvironment.o BioFVM_solvers.o \
BioFVM_matlab.o BioFVM_utilities.o BioFVM_basic_agent.o BioFVM_MultiCellDS.o \
BioFVM_agent_container.o

PhysiCell_core_OBJECTS := PhysiCell_phenotype.o PhysiCell_cell_container.o \
PhysiCell_standard_models.o PhysiCell_cell.o PhysiCell_custom.o PhysiCell_utilities.o \
PhysiCell_constants.o PhysiCell_basic_signaling.o PhysiCell_signal_behavior.o PhysiCell_rules.o

PhysiCell_module_OBJECTS := PhysiCell_SVG.o PhysiCell_pathology.o PhysiCell_MultiCellDS.o \
PhysiCell_various_outputs.o PhysiCell_pugixml.o PhysiCell_settings.o PhysiCell_geometry.o

# ---------------------------------------------------------------
# dFBA addon objects
# ---------------------------------------------------------------
dFBA_OBJECTS := dfba_Metabolite.o dfba_Reaction.o dfba_Solution.o dfba_Model.o \
dfba_intracellular.o

# ---------------------------------------------------------------
# Custom module objects
# ---------------------------------------------------------------
PhysiCell_custom_module_OBJECTS := custom.o

pugixml_OBJECTS := pugixml.o

PhysiCell_OBJECTS := $(BioFVM_OBJECTS) $(pugixml_OBJECTS) $(PhysiCell_core_OBJECTS) \
$(PhysiCell_module_OBJECTS)

ALL_OBJECTS := $(PhysiCell_OBJECTS) $(PhysiCell_custom_module_OBJECTS) $(dFBA_OBJECTS)

# ---------------------------------------------------------------
# Main build target
# ---------------------------------------------------------------
all: main.cpp $(ALL_OBJECTS)
	$(LINK_COMMAND) -o $(PROGRAM_NAME) $(ALL_OBJECTS) main.cpp \
	$(DFBA_LD) $(DFBA_RPATH) $(DFBA_LIBS) -fopenmp
	make name

name:
	@echo ""
	@echo "Executable name is $(PROGRAM_NAME)"
	@echo ""

# ---------------------------------------------------------------
# dFBA addon sources
# ---------------------------------------------------------------
dfba_intracellular.o: $(DFBA_SRC)/dfba_intracellular.cpp
	$(COMPILE_COMMAND) -c $< -o $@

dfba_Model.o: $(DFBA_SRC)/dfba_Model.cpp
	$(COMPILE_COMMAND) -c $< -o $@

dfba_Reaction.o: $(DFBA_SRC)/dfba_Reaction.cpp
	$(COMPILE_COMMAND) -c $< -o $@

dfba_Metabolite.o: $(DFBA_SRC)/dfba_Metabolite.cpp
	$(COMPILE_COMMAND) -c $< -o $@

dfba_Solution.o: $(DFBA_SRC)/dfba_Solution.cpp
	$(COMPILE_COMMAND) -c $< -o $@

# ---------------------------------------------------------------
# PhysiCell core components
# ---------------------------------------------------------------
PhysiCell_phenotype.o: ./core/PhysiCell_phenotype.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_phenotype.cpp

PhysiCell_digital_cell_line.o: ./core/PhysiCell_digital_cell_line.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_digital_cell_line.cpp

PhysiCell_cell.o: ./core/PhysiCell_cell.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_cell.cpp

PhysiCell_cell_container.o: ./core/PhysiCell_cell_container.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_cell_container.cpp

PhysiCell_standard_models.o: ./core/PhysiCell_standard_models.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_standard_models.cpp

PhysiCell_utilities.o: ./core/PhysiCell_utilities.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_utilities.cpp

PhysiCell_custom.o: ./core/PhysiCell_custom.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_custom.cpp

PhysiCell_constants.o: ./core/PhysiCell_constants.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_constants.cpp

PhysiCell_signal_behavior.o: ./core/PhysiCell_signal_behavior.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_signal_behavior.cpp

PhysiCell_rules.o: ./core/PhysiCell_rules.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_rules.cpp

PhysiCell_basic_signaling.o: ./core/PhysiCell_basic_signaling.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_basic_signaling.cpp

# ---------------------------------------------------------------
# BioFVM core components
# ---------------------------------------------------------------
BioFVM_vector.o: ./BioFVM/BioFVM_vector.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_vector.cpp

BioFVM_agent_container.o: ./BioFVM/BioFVM_agent_container.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_agent_container.cpp

BioFVM_mesh.o: ./BioFVM/BioFVM_mesh.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_mesh.cpp

BioFVM_microenvironment.o: ./BioFVM/BioFVM_microenvironment.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_microenvironment.cpp

BioFVM_solvers.o: ./BioFVM/BioFVM_solvers.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_solvers.cpp

BioFVM_utilities.o: ./BioFVM/BioFVM_utilities.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_utilities.cpp

BioFVM_basic_agent.o: ./BioFVM/BioFVM_basic_agent.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_basic_agent.cpp

BioFVM_matlab.o: ./BioFVM/BioFVM_matlab.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_matlab.cpp

BioFVM_MultiCellDS.o: ./BioFVM/BioFVM_MultiCellDS.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_MultiCellDS.cpp

pugixml.o: ./BioFVM/pugixml.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/pugixml.cpp

# ---------------------------------------------------------------
# Standard PhysiCell modules
# ---------------------------------------------------------------
PhysiCell_SVG.o: ./modules/PhysiCell_SVG.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_SVG.cpp

PhysiCell_pathology.o: ./modules/PhysiCell_pathology.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_pathology.cpp

PhysiCell_MultiCellDS.o: ./modules/PhysiCell_MultiCellDS.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_MultiCellDS.cpp

PhysiCell_various_outputs.o: ./modules/PhysiCell_various_outputs.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_various_outputs.cpp

PhysiCell_pugixml.o: ./modules/PhysiCell_pugixml.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_pugixml.cpp

PhysiCell_settings.o: ./modules/PhysiCell_settings.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_settings.cpp

PhysiCell_geometry.o: ./modules/PhysiCell_geometry.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_geometry.cpp

# ---------------------------------------------------------------
# Custom modules
# ---------------------------------------------------------------
custom.o: ./custom_modules/custom.cpp
	$(COMPILE_COMMAND) -c ./custom_modules/custom.cpp

# ---------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------
clean:
	rm -f *.o
	rm -f $(PROGRAM_NAME)*

data-cleanup:
	rm -rf ./output
	mkdir ./output
	touch ./output/empty.txt

# ---------------------------------------------------------------
# Archival
# ---------------------------------------------------------------
checkpoint:
	zip -r $$(date +%b_%d_%Y_%H%M).zip Makefile *.cpp *.h config/*.xml custom_modules/*

zip:
	zip -r latest.zip Makefile* *.cpp *.h BioFVM/* config/* core/* custom_modules/* \
	matlab/* modules/* sample_projects/* addons/dFBA/src/* addons/dFBA/ext/*
