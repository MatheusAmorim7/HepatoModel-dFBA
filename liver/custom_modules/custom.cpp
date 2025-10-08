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
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
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
    double o2_max = 38.0;      // oxigênio máximo (zona 1)
    double o2_min = 5.0;      // oxigênio mínimo (zona 3, centro)

    for (int n = 0; n < microenvironment.number_of_voxels(); n++)
    {
        double x = microenvironment.mesh.voxels[n].center[0];
        double y = microenvironment.mesh.voxels[n].center[1];

        // distância ao centro (2D)
        double dist = sqrt((x - x_center)*(x - x_center) + (y - y_center)*(y - y_center));

        // gradiente linear de O2 do centro para a periferia
        double o2 = o2_max - (o2_max - o2_min) * (radius_max - dist) / radius_max;

        // limitar para não ficar abaixo de o2_min
        if (o2 < o2_min) o2 = o2_min;

        microenvironment(n)[oxygen_index] = o2;
    }
}

// --- dentro do setup_microenvironment() ---
void setup_microenvironment()
{
    // inicializa microambiente
    initialize_microenvironment();

    initialize_oxygen_zones();

    return;
}

void dynamic_oxygen_supply(double dt)
{
    int oxygen_index = microenvironment.find_density_index("oxygen");

    // Parâmetros de ajuste
    double limite = 20.0;     // valor crítico para ativar injeção
    double max_o2 = 38.0;     // valor alvo máximo
    double raio_influencia = 100.0; // distância máxima da fonte para injetar O2

    // Coordenadas das fontes periféricas
    std::vector<std::array<double,3>> fontes;
    double x_center = 0.5 * (microenvironment.mesh.bounding_box[0] + microenvironment.mesh.bounding_box[3]);
    double y_center = 0.5 * (microenvironment.mesh.bounding_box[1] + microenvironment.mesh.bounding_box[4]);
    double z_center = 0.5 * (microenvironment.mesh.bounding_box[2] + microenvironment.mesh.bounding_box[5]);
    double radius_max = 350.0;
    int n_sources = 6;

    for(int i = 0; i < n_sources; i++)
    {
        double theta = i * 2.0 * M_PI / n_sources;
        double x = x_center + radius_max * cos(theta);
        double y = y_center + radius_max * sin(theta);
        fontes.push_back({x, y, z_center});
    }

    for (int n = 0; n < microenvironment.number_of_voxels(); n++)
    {
        double& o2 = microenvironment(n)[oxygen_index];

        if(o2 >= limite) continue; // só aumenta voxels abaixo do limite

        double voxel_x = microenvironment.mesh.voxels[n].center[0];
        double voxel_y = microenvironment.mesh.voxels[n].center[1];
        double voxel_z = microenvironment.mesh.voxels[n].center[2];

        // Calcula contribuição de cada fonte
        double incremento = 0.0;
        for(auto &f : fontes)
        {
            double dx = voxel_x - f[0];
            double dy = voxel_y - f[1];
            double dz = voxel_z - f[2];
            double dist = sqrt(dx*dx + dy*dy + dz*dz);

            if(dist <= raio_influencia)
            {
                double peso = 1.0 - dist/raio_influencia; // menor O2 quanto mais longe da fonte
                incremento += peso;
            }
        }

        // Normaliza incremento e aplica ao voxel
        if(incremento > 0)
        {
            double deficit = max_o2 - o2;
            double taxa_injecao = 0.5 * deficit * incremento * dt; // ajustável
            o2 += taxa_injecao;
            if(o2 > max_o2) o2 = max_o2;
        }
    }
}

void setup_peripheral_sources()
{
    int oxygen_index = microenvironment.find_density_index("oxygen");

    double x_center = 0.5 * (microenvironment.mesh.bounding_box[0] + microenvironment.mesh.bounding_box[3]);
    double y_center = 0.5 * (microenvironment.mesh.bounding_box[1] + microenvironment.mesh.bounding_box[4]);
    double z_center = 0.5 * (microenvironment.mesh.bounding_box[2] + microenvironment.mesh.bounding_box[5]);

    double radius_max = 350.0; // periferia
    int n_sources = 6;

    for(int i = 0; i < n_sources; i++)
    {
        double theta = i * 2.0 * M_PI / n_sources;

        double x = x_center + radius_max * cos(theta);
        double y = y_center + radius_max * sin(theta);

        Cell* source = create_cell();
        source->assign_position({x, y, z_center});

        source->phenotype.secretion.uptake_rates[oxygen_index] = 0.0;    // não consome
        source->phenotype.secretion.secretion_rates[oxygen_index] = 1.0; // secreta oxigênio
        source->phenotype.secretion.saturation_densities[oxygen_index] = 38.0; // máximo
    }
}

void setup_tissue(void) 
{ 
    int oxygen_substrate = microenvironment.find_density_index("oxygen");

    double x_center = 0.5 * (microenvironment.mesh.bounding_box[0] + microenvironment.mesh.bounding_box[3]);
    double y_center = 0.5 * (microenvironment.mesh.bounding_box[1] + microenvironment.mesh.bounding_box[4]);
    double z_center = 0.5 * (microenvironment.mesh.bounding_box[2] + microenvironment.mesh.bounding_box[5]);

    // --- Veia central: super-sink ---
    Cell* sink = create_cell();
    sink->assign_position({x_center, y_center, z_center});

    sink->phenotype.secretion.uptake_rates[oxygen_substrate] = 5.0; //Diminuir para maior hipoxia central(Faz demorar consumo de O2 no centro)
    sink->phenotype.secretion.secretion_rates[oxygen_substrate] = 0.0;
    sink->phenotype.secretion.saturation_densities[oxygen_substrate] = 0.0;

    setup_peripheral_sources();

}


void custom_microenvironment_function(double dt)
{
    dynamic_oxygen_supply(dt);
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 