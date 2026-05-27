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
    double o2_max = 0.05;      // oxigênio máximo (zona 1)
    double o2_min = 0.007;      // oxigênio mínimo (zona 3, centro)

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

    double limite  = 0.03;
    double max_o2  = 0.05;
    double raio_influencia = 500.0;

    double x_center = 0.5 * (microenvironment.mesh.bounding_box[0] + microenvironment.mesh.bounding_box[3]);
    double y_center = 0.5 * (microenvironment.mesh.bounding_box[1] + microenvironment.mesh.bounding_box[4]);
    double z_center = 0.5 * (microenvironment.mesh.bounding_box[2] + microenvironment.mesh.bounding_box[5]);

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

        // borda fixa no máximo (vasculatura portal)
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
            double taxa = 2.0 * (max_o2 - o2) * incremento * dt;
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
        source->phenotype.secretion.secretion_rates[oxygen_index] = 1.0; // secreta oxigênio
        source->phenotype.secretion.saturation_densities[oxygen_index] = 0.05; // máximo

        source->phenotype.secretion.uptake_rates[glucose_index] = 0.0;
        source->phenotype.secretion.secretion_rates[glucose_index] = 1.0;
        source->phenotype.secretion.saturation_densities[glucose_index] = 5.0; // mM
    }
}

void setup_tissue(void) 
{ 
    int oxygen_substrate = microenvironment.find_density_index("oxygen");

    Cell_Definition* pSourceDef = find_cell_definition("source_sink");

    double x_center = 0.5 * (microenvironment.mesh.bounding_box[0] + microenvironment.mesh.bounding_box[3]);
    double y_center = 0.5 * (microenvironment.mesh.bounding_box[1] + microenvironment.mesh.bounding_box[4]);
    double z_center = 0.5 * (microenvironment.mesh.bounding_box[2] + microenvironment.mesh.bounding_box[5]);

    // --- Veia central: super-sink ---
    Cell* sink = create_cell(*pSourceDef);
    sink->assign_position({x_center, y_center, z_center});

    sink->phenotype.intracellular = nullptr;

    sink->phenotype.secretion.uptake_rates[oxygen_substrate] = 0.5; //Diminuir para maior hipoxia central(Faz demorar consumo de O2 no centro)
    sink->phenotype.secretion.secretion_rates[oxygen_substrate] = 0.0;
    sink->phenotype.secretion.saturation_densities[oxygen_substrate] = 0.0;

    setup_peripheral_sources();

}

void dynamic_glucose_supply(double dt)
{
    int glucose_index = microenvironment.find_density_index("glucose");

    double limite   = 3.0;
    double max_glc  = 5.0;
    double raio_influencia = 500.0;  // igual ao O2

    double x_center = 0.5 * (microenvironment.mesh.bounding_box[0] + microenvironment.mesh.bounding_box[3]);
    double y_center = 0.5 * (microenvironment.mesh.bounding_box[1] + microenvironment.mesh.bounding_box[4]);
    double z_center = 0.5 * (microenvironment.mesh.bounding_box[2] + microenvironment.mesh.bounding_box[5]);

    std::vector<std::array<double,3>> fontes;
    int n_sources = 6;
    double radius_max = 450.0;
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

        // borda fixa no máximo (portal)
        bool na_borda = (fabs(voxel_x - x_center) > 400.0 || fabs(voxel_y - y_center) > 400.0);
        if (na_borda) { glc = max_glc; continue; }

        if (glc >= limite) continue;

        double incremento = 0.0;
        for(auto& f : fontes)
        {
            double dist = sqrt(pow(voxel_x-f[0],2) + pow(voxel_y-f[1],2) + pow(voxel_z-f[2],2));
            if(dist <= raio_influencia)
                incremento += (1.0 - dist / raio_influencia);
        }

        if(incremento > 0)
        {
            double taxa = 2.0 * (max_glc - glc) * incremento * dt;
            glc = std::min(glc + taxa, max_glc);
        }
    }
}

void custom_microenvironment_function(double dt)
{
    dynamic_oxygen_supply(dt);
    dynamic_glucose_supply(dt);
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 