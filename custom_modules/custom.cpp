/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"
#include "../BioFVM/BioFVM.h"

using namespace BioFVM;

// Cell-definition defaults used by the hepatic lobule dFBA native model.
// The legacy sinusoidal-supply driver (setup_tissue, oxygen/glucose supply,
// zones) lives in legacy/custom_modules/custom.cpp.

// =============================================================================
// fixed_volume_update_function: substitui standard_volume_update_function.
//
// O dFBA chama multiply_by_ratio() em update_dfba_outputs(), que aumenta o
// volume diretamente. Depois, o PhysiCell chama volume_update_function() para
// atualizar sub-volumes e verificar divisão.
//
// Aqui nós resetamos o volume para o valor fixo de referência ANTES que o
// PhysiCell verifique se a célula deve dividir.
// Isso preserva todos os fluxos dFBA sem alterar o solver.
//
// ORDEM REAL DE EXECUÇÃO (por célula, a cada dt):
//   1. phenotype_function()
//   2. custom_function()
//   3. intracellular->update()     <- dFBA: multiply_by_ratio aumenta volume
//   4. volume_update_function()    <- AQUI: resetamos para TARGET_VOLUME
//   5. check_for_division()        <- volume já está no alvo -> sem divisão
// =============================================================================
static const double HEPATOCYTE_TARGET_VOLUME = 2494.0; // um³ — reference_volume do XML

void fixed_volume_update_function( Cell* pCell, Phenotype& phenotype, double dt )
{
    // Células source_sink não têm dFBA — deixa o comportamento padrão
    if( pCell->phenotype.intracellular == nullptr )
    {
        standard_volume_update_function( pCell, phenotype, dt );
        return;
    }

    // Reseta volume total para o valor fixo de referência
    double target = HEPATOCYTE_TARGET_VOLUME;

    phenotype.volume.total             = target;
    phenotype.volume.nuclear           = 0.1 * target;
    phenotype.volume.nuclear_fluid     = phenotype.volume.nuclear * phenotype.volume.fluid_fraction;
    phenotype.volume.nuclear_solid     = phenotype.volume.nuclear - phenotype.volume.nuclear_fluid;
    phenotype.volume.cytoplasmic       = target - phenotype.volume.nuclear;
    phenotype.volume.cytoplasmic_fluid = phenotype.volume.cytoplasmic * phenotype.volume.fluid_fraction;
    phenotype.volume.cytoplasmic_solid = phenotype.volume.cytoplasmic - phenotype.volume.cytoplasmic_fluid;
    phenotype.volume.solid             = phenotype.volume.nuclear_solid + phenotype.volume.cytoplasmic_solid;
    phenotype.volume.fluid             = phenotype.volume.nuclear_fluid + phenotype.volume.cytoplasmic_fluid;

    pCell->set_total_volume( target );
    phenotype.geometry.update( pCell, phenotype, dt );

    // Zera transições de ciclo por garantia (evita divisão residual)
    int n_phases = (int) phenotype.cycle.model().phases.size();
    for( int i = 0; i < n_phases; i++ )
        for( int j = 0; j < n_phases; j++ )
            phenotype.cycle.data.transition_rate(i, j) = 0.0;
}


void create_cell_types( void )
{
	// set the random seed 
	if (parameters.ints.find_index("random_seed") != -1)
	{
		SeedRandom(parameters.ints("random_seed"));
	}
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = fixed_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/*
       Cell rule definitions 
	*/

	setup_cell_rules(); 

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

// set up the BioFVM microenvironment
void setup_microenvironment()
{
    initialize_microenvironment();
    return;
}

static void zero_all_cycle_transitions( Phenotype& phenotype )
{
    // phenotype.cycle.model().phases é o std::vector<Phase> do modelo de ciclo.
    // É a forma correta de obter o número de fases no PhysiCell.
    int n_phases = (int) phenotype.cycle.model().phases.size();
    for( int i = 0; i < n_phases; i++ )
        for( int j = 0; j < n_phases; j++ )
            phenotype.cycle.data.transition_rate(i, j) = 0.0;
}

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{
    // phenotype_function roda ANTES do update intracellular (dFBA).
    // NÃO bloqueamos crescimento aqui — o dFBA sobrescreveria depois.
    // O bloqueio real fica em custom_function(), que roda APÓS o dFBA.
    return;
}

// =============================================================================
// custom_function: roda DEPOIS do update intracellular (dFBA) a cada passo.
// É aqui que sobrescrevemos o growth_rate que o dFBA escreveu no ciclo.
//
// ORDEM DE EXECUÇÃO NO PhysiCell (por célula, a cada dt):
//   1. phenotype_function()          ← antes do dFBA  (não usamos para ciclo)
//   2. intracellular->update()       ← dFBA roda, escreve growth_rate no ciclo
//   3. custom_function()             ← AQUI zeramos a transição após o dFBA
//   4. standard_volume_update()      ← volume cresce com base na transição
//
// CICLO "Live" (code=5): tem APENAS 1 fase (Phase 0: Live → Live).
//   A transição 0→0 É a taxa de proliferação neste modelo.
//   growth_rate=1 do dFBA vira transition_rate(0,0)=1 → célula dobra volume.
//   Zeramos transition_rate(0,0)=0 aqui para bloquear completamente.
// =============================================================================
void custom_function( Cell* pCell, Phenotype& phenotype, double dt )
{
    // Volume e ciclo são controlados em fixed_volume_update_function().
    // Esta função fica disponível para regras adicionais (quimiotaxia, etc.).
    return;
}

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; }
