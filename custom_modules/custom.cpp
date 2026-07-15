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
#include <array> 
#include <cmath> 
#include <vector>

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

void initialize_oxygen_zones()
{
    int oxygen_index = microenvironment.find_density_index("oxygen");

    // centro do domínio
    double x_center = 0.5 * (microenvironment.mesh.bounding_box[0] + microenvironment.mesh.bounding_box[3]);
    double y_center = 0.5 * (microenvironment.mesh.bounding_box[1] + microenvironment.mesh.bounding_box[4]);

    // parâmetros do gradiente
    double radius_max = 330.0; // alcance da zona 1 (periferia)
    double o2_max = 0.13;      // oxigênio máximo (zona 1)
    double o2_min = 0.02;      // oxigênio mínimo (zona 3, centro)

    for (int n = 0; n < microenvironment.number_of_voxels(); n++)
    {
        double x = microenvironment.mesh.voxels[n].center[0];
        double y = microenvironment.mesh.voxels[n].center[1];

        // distância ao centro (2D)
        double dist = sqrt((x - x_center)*(x - x_center) + (y - y_center)*(y - y_center));

        double frac = std::min(dist / radius_max, 1.0);
        microenvironment(n)[oxygen_index] = o2_min + (o2_max - o2_min) * frac;
    }
}

void initialize_glucose_zones()
{
    int glucose_index = microenvironment.find_density_index("glucose");

    double x_center = 0.5 * (microenvironment.mesh.bounding_box[0] + microenvironment.mesh.bounding_box[3]);
    double y_center = 0.5 * (microenvironment.mesh.bounding_box[1] + microenvironment.mesh.bounding_box[4]);

    double radius_max = 330.0;
    double glc_max = 5.0;   // mM portal (periferia)
    double glc_min = 1.0;   // mM veia central

    for (int n = 0; n < microenvironment.number_of_voxels(); n++)
    {
        double x = microenvironment.mesh.voxels[n].center[0];
        double y = microenvironment.mesh.voxels[n].center[1];
        double dist = sqrt((x - x_center)*(x - x_center) + (y - y_center)*(y - y_center));

        double frac = std::min(dist / radius_max, 1.0);
        microenvironment(n)[glucose_index] = glc_min + (glc_max - glc_min) * frac;
    }
}

// --- dentro do setup_microenvironment() ---
void setup_microenvironment()
{
    // inicializa microambiente
    initialize_microenvironment();

    initialize_oxygen_zones();

    initialize_glucose_zones();

    return;
}

void dynamic_oxygen_supply(double dt)
{
    int oxygen_index = microenvironment.find_density_index("oxygen");

    double limite  = 0.08;
    double max_o2  = 0.13;
    double raio_influencia = 500.0;

    double x_center = 0.5 * (microenvironment.mesh.bounding_box[0] + microenvironment.mesh.bounding_box[3]);
    double y_center = 0.5 * (microenvironment.mesh.bounding_box[1] + microenvironment.mesh.bounding_box[4]);
    double z_center = 0.5 * (microenvironment.mesh.bounding_box[2] + microenvironment.mesh.bounding_box[5]);

    // escala a injeção com o número de células
    int n_cells = (*all_cells).size();
    double fator_escala = std::max(1.0, n_cells / 9.0);

    std::vector<std::array<double,3>> fontes;
    double radius_max = 450.0;
    int n_sources = 6;
    for(int i = 0; i < n_sources; i++)
    {
        double theta = i * 2.0 * M_PI / n_sources;
        fontes.push_back({x_center + radius_max * cos(theta),
                          y_center + radius_max * sin(theta),
                          z_center});
    }

    for (int n = 0; n < microenvironment.number_of_voxels(); n++)
    {
        double voxel_x = microenvironment.mesh.voxels[n].center[0];
        double voxel_y = microenvironment.mesh.voxels[n].center[1];
        double voxel_z = microenvironment.mesh.voxels[n].center[2];

        double& o2 = microenvironment(n)[oxygen_index];

        bool na_borda = (fabs(voxel_x - x_center) > 400.0 || fabs(voxel_y - y_center) > 400.0);
        if (na_borda) { o2 = max_o2; continue; }

        if (o2 >= limite) continue;

        double incremento = 0.0;
        for(auto& f : fontes)
        {
            double dx = voxel_x - f[0];
            double dy = voxel_y - f[1];
            double dz = voxel_z - f[2];
            double dist = sqrt(dx*dx + dy*dy + dz*dz);
            if(dist <= raio_influencia)
                incremento += (1.0 - dist / raio_influencia);
        }

        if(incremento > 0)
        {
            double taxa = fator_escala * 2.0 * (max_o2 - o2) * incremento * dt;
            o2 = std::min(o2 + taxa, max_o2);
        }
    }
}

void setup_peripheral_sources()
{
    int oxygen_index = microenvironment.find_density_index("oxygen");
    int glucose_index = microenvironment.find_density_index("glucose");

    Cell_Definition* pSourceDef = find_cell_definition("source_sink");

    double x_center = 0.5 * (microenvironment.mesh.bounding_box[0] + microenvironment.mesh.bounding_box[3]);
    double y_center = 0.5 * (microenvironment.mesh.bounding_box[1] + microenvironment.mesh.bounding_box[4]);
    double z_center = 0.5 * (microenvironment.mesh.bounding_box[2] + microenvironment.mesh.bounding_box[5]);

    double radius_max = 450.0; // periferia
    int n_sources = 6;

    for(int i = 0; i < n_sources; i++)
    {
        double theta = i * 2.0 * M_PI / n_sources;

        double x = x_center + radius_max * cos(theta);
        double y = y_center + radius_max * sin(theta);

        Cell* source = create_cell(*pSourceDef);
        source->assign_position({x, y, z_center});

        source->phenotype.intracellular = nullptr;

        source->phenotype.secretion.uptake_rates[oxygen_index] = 0.0;    // não consome
        source->phenotype.secretion.secretion_rates[oxygen_index] = 0.2; // secreta oxigênio
        source->phenotype.secretion.saturation_densities[oxygen_index] = 0.13; // máximo

        source->phenotype.secretion.uptake_rates[glucose_index] = 0.0;
        source->phenotype.secretion.secretion_rates[glucose_index] = 0.2;
        source->phenotype.secretion.saturation_densities[glucose_index] = 5.0; // mM
    }
}

void setup_tissue(void) 
{ 
    int oxygen_substrate = microenvironment.find_density_index("oxygen");
    int glucose_index = microenvironment.find_density_index("glucose");
    int co2_index = microenvironment.find_density_index("CO2");

    Cell_Definition* pSourceDef = find_cell_definition("source_sink");

    double x_center = 0.5 * (microenvironment.mesh.bounding_box[0] + microenvironment.mesh.bounding_box[3]);
    double y_center = 0.5 * (microenvironment.mesh.bounding_box[1] + microenvironment.mesh.bounding_box[4]);
    double z_center = 0.5 * (microenvironment.mesh.bounding_box[2] + microenvironment.mesh.bounding_box[5]);

    // --- Veia central: super-sink ---
    Cell* sink = create_cell(*pSourceDef);
    sink->assign_position({x_center, y_center, z_center});

    sink->phenotype.intracellular = nullptr;

    sink->phenotype.secretion.uptake_rates[oxygen_substrate] = 0.05; //Diminuir para maior hipoxia central(Faz demorar consumo de O2 no centro)
    sink->phenotype.secretion.secretion_rates[oxygen_substrate] = 0.0;
    sink->phenotype.secretion.saturation_densities[oxygen_substrate] = 0.0;

    // adicionar glicose no sink
    sink->phenotype.secretion.uptake_rates[glucose_index]  = 0.0;  // /min
    sink->phenotype.secretion.secretion_rates[glucose_index]  = 0.0;
    sink->phenotype.secretion.saturation_densities[glucose_index]  = 0.0;

    // Dreno de CO2: representa a drenagem venosa do CO2 metabólico produzido no lóbulo.
    // O sink não tem dFBA (intracellular = nullptr), então o consumo aqui é puramente
    // fenomenológico (uptake de primeira ordem do PhysiCell), não uma reação de troca do LP.
    sink->phenotype.secretion.uptake_rates[co2_index] = 0.0005;   // /min -- ponto de partida, calibrar conforme item abaixo
    sink->phenotype.secretion.secretion_rates[co2_index] = 0.0;
    sink->phenotype.secretion.saturation_densities[co2_index] = 0.0;

    setup_peripheral_sources();

}

void dynamic_glucose_supply(double dt)
{
    int glucose_index = microenvironment.find_density_index("glucose");

    double limite   = 2.0;
    double max_glc  = 5.0;
    double raio_influencia = 500.0;

    double x_center = 0.5 * (microenvironment.mesh.bounding_box[0] + microenvironment.mesh.bounding_box[3]);
    double y_center = 0.5 * (microenvironment.mesh.bounding_box[1] + microenvironment.mesh.bounding_box[4]);
    double z_center = 0.5 * (microenvironment.mesh.bounding_box[2] + microenvironment.mesh.bounding_box[5]);

    int n_cells = (*all_cells).size();
    double fator_escala = std::max(1.0, n_cells / 9.0);

    std::vector<std::array<double,3>> fontes;
    double radius_max = 450.0;
    int n_sources = 6;
    for(int i = 0; i < n_sources; i++)
    {
        double theta = i * 2.0 * M_PI / n_sources;
        fontes.push_back({x_center + radius_max * cos(theta),
                          y_center + radius_max * sin(theta),
                          z_center});
    }

    for(int n = 0; n < microenvironment.number_of_voxels(); n++)
    {
        double voxel_x = microenvironment.mesh.voxels[n].center[0];
        double voxel_y = microenvironment.mesh.voxels[n].center[1];
        double voxel_z = microenvironment.mesh.voxels[n].center[2];

        double& glc = microenvironment(n)[glucose_index];

        bool na_borda = (fabs(voxel_x - x_center) > 400.0 || fabs(voxel_y - y_center) > 400.0);
        if (na_borda) { glc = max_glc; continue; }

        if (glc >= limite) continue;

        double incremento = 0.0;
        for(auto& f : fontes)
        {
            double dx = voxel_x - f[0];
            double dy = voxel_y - f[1];
            double dz = voxel_z - f[2];
            double dist = sqrt(dx*dx + dy*dy + dz*dz);
            if(dist <= raio_influencia)
                incremento += (1.0 - dist / raio_influencia);
        }

        if(incremento > 0)
        {
            double taxa = fator_escala * 2.0 * (max_glc - glc) * incremento * dt;
            glc = std::min(glc + taxa, max_glc);
        }
    }
}

void sinusoidal_distributed_supply(double dt)
{
    int oxygen_index = microenvironment.find_density_index("oxygen");
    int glucose_index = microenvironment.find_density_index("glucose");

    if( oxygen_index < 0 && glucose_index < 0 )
    { return; }

    double x_center = 0.5 * (microenvironment.mesh.bounding_box[0] + microenvironment.mesh.bounding_box[3]);
    double y_center = 0.5 * (microenvironment.mesh.bounding_box[1] + microenvironment.mesh.bounding_box[4]);

    // Level 1 sinusoidal perfusion approximation:
    // weak radial relaxation toward portal-to-central blood-side targets.
    const double radius_max = 450.0;

    const double o2_portal_target = 0.12;   // mM, periportal/sinusoidal inlet side
    const double o2_central_target = 0.035; // mM, pericentral/sinusoidal outlet side
    const double glc_portal_target = 5.0;   // mM
    const double glc_central_target = 3.5;  // mM

    const double k_o2 = 0.02;  // 1/min
    const double k_glc = 0.01; // 1/min

    for( int n = 0; n < microenvironment.number_of_voxels(); n++ )
    {
        double voxel_x = microenvironment.mesh.voxels[n].center[0];
        double voxel_y = microenvironment.mesh.voxels[n].center[1];

        double dx = voxel_x - x_center;
        double dy = voxel_y - y_center;
        double distance_from_center = sqrt( dx*dx + dy*dy );
        double portal_fraction = std::min( distance_from_center / radius_max, 1.0 );

        if( oxygen_index >= 0 )
        {
            double oxygen_target = o2_central_target
                + ( o2_portal_target - o2_central_target ) * portal_fraction;
            double& oxygen = microenvironment(n)[oxygen_index];
            oxygen += k_o2 * ( oxygen_target - oxygen ) * dt;
            oxygen = std::max( 0.0, std::min( oxygen, o2_portal_target ) );
        }

        if( glucose_index >= 0 )
        {
            double glucose_target = glc_central_target
                + ( glc_portal_target - glc_central_target ) * portal_fraction;
            double& glucose = microenvironment(n)[glucose_index];
            glucose += k_glc * ( glucose_target - glucose ) * dt;
            glucose = std::max( 0.0, std::min( glucose, glc_portal_target ) );
        }
    }
}

void custom_microenvironment_function(double dt)
{
    sinusoidal_distributed_supply(dt);
    return;
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

// =============================================================================
// CONTROLE DE PROLIFERAÇÃO VIA dFBA — OPÇÃO 3
//
// LÓGICA GERAL:
//   Hepatócitos adultos são quiescentes. Não queremos proliferação livre.
//   Usamos o fluxo da reação objetivo de cada zona como "sinal metabólico":
//
//   ZONA 1/2 — hepatonet.xml
//     Objetivo: R_r1032 (minimize) = Glucose(s) → Glucose(c)
//     Sentido biológico: a célula IMPORTA glicose do espaço sinusoidal.
//     O fluxo otimizado é NEGATIVO (fluxo reverso = liberação de glicose),
//     ou seja, quanto mais negativo, mais glicose a célula está EXPORTANDO.
//     Usamos |fluxo| como proxy de atividade metabólica.
//
//   ZONA 3 — hepatonet_zone3.xml
//     Objetivo: R_r1389 (minimize, lb=1) = Glycogenin-G11 → Glycogenin-G4G7
//     Sentido biológico: síntese de glicogênio.
//     O fluxo é POSITIVO (armazenamento ativo). Quanto maior, mais síntese.
//
// THRESHOLD:
//   Se |fluxo_objetivo| < threshold  →  célula quiescente (sem divisão)
//   Se |fluxo_objetivo| >= threshold →  permite divisão, mas com taxa baixa
//
//   O threshold evita que ruído numérico do solver dispare divisão.
//   Recomendado: 0.01 mmol/gDW/hr para começar, ajuste conforme outputs.
//
// ONDE FICA:
//   - phenotype_function(): lê o fluxo do intracellular e seta a taxa de
//     transição de ciclo. Roda a cada dt para CADA célula.
//   - custom_function(): deixamos vazio (pode usar para outras regras).
// =============================================================================

// =============================================================================
// CONFIGURAÇÕES DE QUIESCÊNCIA
// =============================================================================

// Threshold mínimo de fluxo para considerar a célula metabolicamente ativa.
// Usado apenas para logging/debug — crescimento está SEMPRE bloqueado abaixo.
static const double FLUX_THRESHOLD_Z12 = 0.01;  // mmol/gDW/hr — Zona 1/2
static const double FLUX_THRESHOLD_Z3  = 0.5;   // mmol/gDW/hr — Zona 3

// Taxa de crescimento máxima permitida.
// 0.0 = hepatócitos completamente quiescentes (divisão bloqueada).
// ~7e-4 = divisão em ~24h se quiser ativar proliferação lenta no futuro.
static const double MAX_GROWTH_RATE_ACTIVE = 0.0;

// =============================================================================
// FUNÇÃO AUXILIAR: zera TODAS as taxas de transição do ciclo celular.
// Necessário porque o dFBA (via roadrunner) escreve growth_rate no ciclo
// ANTES do phenotype_function rodar — precisamos sobrescrever depois.
// =============================================================================
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