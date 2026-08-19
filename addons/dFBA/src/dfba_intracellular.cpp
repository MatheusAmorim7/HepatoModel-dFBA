#include "dfba_Reaction.h"
#include "dfba_Model.h"
#include "dfba_intracellular.h"

#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <mutex>


namespace PhysiCelldFBA {

dFBAIntracellular::dFBAIntracellular() : Intracellular()
{
    intracellular_type = "dfba";
    objective_reaction = "";
    sbml_filename = "";
    cell_density = 0.0;
    reference_volume = 0.0;
    max_growth_rate = 0.0;
    current_growth_rate = 0.0;
    last_solution_status = "not_run";
    next_dfba_run = 0.0;
    use_metabolic_death = true;
    death_type = "";
    death_trigger_flux = "";
    death_flux_threshold = 0.0;
    death_rate_increase = 0.0;
    substrate_exchanges.clear();
    flag_for_death = false;
    dfba_time_step = PhysiCell::diffusion_dt; // Default value
}

dFBAIntracellular::dFBAIntracellular(pugi::xml_node& node)
{
    intracellular_type = "dfba";
    objective_reaction = "";
    sbml_filename = "";
    cell_density = 0.0;
    reference_volume = 0.0;
    max_growth_rate = 0.0;
    current_growth_rate = 0.0;
    last_solution_status = "not_run";
    next_dfba_run = 0.0;
    use_metabolic_death = true;
    death_type = "";
    death_trigger_flux = "";
    death_flux_threshold = 0.0;
    death_rate_increase = 0.0;
    substrate_exchanges.clear();
    flag_for_death = false;
    sbml_model.clear();
    is_initialized = false;
    dfba_time_step = PhysiCell::diffusion_dt; // Default value
    this->initialize_intracellular_from_pugixml(node);
}

dFBAIntracellular::dFBAIntracellular(const dFBAIntracellular& copy) : Intracellular() {
    intracellular_type = copy.intracellular_type;
    objective_reaction = copy.objective_reaction;
    sbml_filename = copy.sbml_filename;
    cell_density = copy.cell_density;
    reference_volume = copy.reference_volume;
    max_growth_rate = copy.max_growth_rate;
    current_growth_rate = copy.current_growth_rate;
    last_solution_status = copy.last_solution_status;
    next_dfba_run = copy.next_dfba_run;
    dfba_time_step = copy.dfba_time_step; // Copy the time step

    // Copy sbml_model
    sbml_model = copy.sbml_model;

    // Copy substrate_exchanges (map is copied by value, which is fine here)
    substrate_exchanges = copy.substrate_exchanges;

    // Copy the initialization flag
    is_initialized = copy.is_initialized;
    use_metabolic_death = copy.use_metabolic_death;
	death_type = copy.death_type;
	death_trigger_flux = copy.death_trigger_flux;
	death_flux_threshold = copy.death_flux_threshold;
	death_rate_increase = copy.death_rate_increase;
	flag_for_death = copy.flag_for_death;
}


int dFBAIntracellular::parse_transport_model(pugi::xml_node& node)
{
    ExchangeFluxData exchange_flux;
    KineticParam Km;
    KineticParam Vmax;
    
    pugi::xml_node node_exchange = node.child( "exchange" );
	int num_exchanges = 0;
    while( node_exchange )
	{
		string density_name = node_exchange.attribute( "substrate" ).value(); 
        int density_index = microenvironment.find_density_index( density_name ); 

        if (density_index == -1)
        {
            std::cout << "Error: attempted to set secretion/uptake/export for \""; 
            std::cout << density_name << "\", which was not found in the microenvironment." << std::endl;
            std::cout << "Please double-check your substrate name in the config file." << std::endl;
            std::cout << std::endl;
            exit(-1); 
        }

        std::string actual_name = microenvironment.density_names[ density_index ]; 
			
        // error check 
        if( std::strcmp( density_name.c_str() , actual_name.c_str() ) != 0 )
        {
            std::cout << "Error: attempted to set secretion/uptake/export for \""; 
            std::cout << density_name << "\", which was not found in the microenvironment." << std::endl;
            std::cout << "Please double-check your substrate name in the config file." << std::endl;
            std::cout << std::endl;
            exit(-1); 
        }
        
        pugi::xml_node node_fba_flux = node_exchange.child( "fba_flux" ); 
		if( node_fba_flux )
		{  
            exchange_flux.fba_flux_id = PhysiCell::xml_get_my_string_value(node_fba_flux);
        }
        else {
            std::cout << "Error: attempted get fba_flux node for "; 
            std::cout << exchange_flux.density_name << "\", but not found." << std::endl;
            std::cout << "Please double-check your exchange nodes in the config file." << std::endl;
            std::cout << std::endl; 
            exit(-1); 
        }

        pugi::xml_node node_Km = node_exchange.child( "Km" ); 
		if( node_Km )
		{
            Km.name = "Km";
            Km.units = node_Km.attribute("units").value();
            Km.value = PhysiCell::xml_get_my_double_value(node_Km);
            assert(Km.value > 0);
            
        }
        else {
            std::cout << "Error: attempted get Km node for "; 
            std::cout << density_name << "\", but not found." << std::endl;
            std::cout << "Please double-check your exchange nodes in the config file." << std::endl;
            std::cout << std::endl; 
            exit(-1); 
        }

        pugi::xml_node node_Vmax = node_exchange.child( "Vmax" ); 
		if( node_Vmax )
		{
            Vmax.name = "Vmax";
            Vmax.units = node_Vmax.attribute("units").value();
            Vmax.value = PhysiCell::xml_get_my_double_value(node_Vmax);
            assert(Vmax.value >= 0);
            
        }
        else {
            std::cout << "Error: attempted get Vmax node for "; 
            std::cout << density_name << "\", but not found." << std::endl;
            std::cout << "Please double-check your exchange nodes in the config file." << std::endl;
            std::cout << std::endl; 
            exit(-1); 
        }

        bool allow_uptake = true;
        pugi::xml_node node_allow_uptake = node_exchange.child( "allow_uptake" );
        if( node_allow_uptake )
        {
            std::string value = PhysiCell::xml_get_my_string_value( node_allow_uptake );
            std::transform( value.begin(), value.end(), value.begin(), ::tolower );
            allow_uptake = !( value == "false" || value == "0" || value == "no" );
        }
		
        
        double minimum_uptake_fraction = 0.0;
        pugi::xml_node node_minimum_uptake_fraction = node_exchange.child( "minimum_uptake_fraction" );
        if( node_minimum_uptake_fraction )
        {
            minimum_uptake_fraction = PhysiCell::xml_get_my_double_value( node_minimum_uptake_fraction );
            minimum_uptake_fraction = std::max( 0.0, std::min( minimum_uptake_fraction, 1.0 ) );
        }

        exchange_flux.density_name = density_name;
        exchange_flux.density_index = density_index;
        exchange_flux.Km = Km;
        exchange_flux.Vmax = Vmax;
        exchange_flux.allow_uptake = allow_uptake;
        exchange_flux.minimum_uptake_fraction = minimum_uptake_fraction;

        this->substrate_exchanges[density_name] = exchange_flux;
        num_exchanges++;
		node_exchange = node_exchange.next_sibling( "exchange" ); 
	}
    return num_exchanges;
}

void dFBAIntracellular::parse_growth_model(pugi::xml_node& parent)
{
    
    pugi::xml_node node = parent.child( "cell_density" );
	if ( node )
	{ 
        this->cell_density = PhysiCell::xml_get_my_double_value(node);
    }
    else
    {
        std::cout << "Error: attempted to read cell_density attribute." << std::endl;
        std::cout << "Please double-check your growth model in the XML setting." << std::endl;
        std::cout << std::endl; 
        exit(-1); 
    }

    node = parent.child( "reference_volume" );
	if ( node )
	{ 
        this->reference_volume = PhysiCell::xml_get_my_double_value(node);
    }
    else
    {
        std::cout << "Error: attempted to read reference_volume attribute." << std::endl;
        std::cout << "Please double-check your reference_volume in the XML setting." << std::endl;
        std::cout << "This value is used to determine the reference volume for cell division." << std::endl;
        std::cout << "Ensure consistency with the total cell volume defined in the cell phenotype." << std::endl;
        std::cout << std::endl; 
        exit(-1); 
    }


    node = parent.child( "max_growth_rate" );
	if ( node )
	{ 
        this->max_growth_rate = PhysiCell::xml_get_my_double_value(node);
        assert(this->max_growth_rate >= 0);
    }
    else
    {
        std::cout << "Error: attempted to read max_growth_rate attribute." << std::endl;
        std::cout << "Please double-check your growth model in the XML setting." << std::endl;
        std::cout << std::endl; 
        exit(-1); 
    }
    
    node = parent.child( "objective_reaction" );
	if ( node )
	{ 
        this->objective_reaction = PhysiCell::xml_get_my_string_value(node);
    }
    else
    {
        std::cout << "Error: attempted to read objective_reaction attribute." << std::endl;
        std::cout << "Please double-check your growth model in the XML setting." << std::endl;
        std::cout << std::endl; 
        exit(-1); 
    }

    
}

void dFBAIntracellular::parse_death_model(pugi::xml_node& parent){


    // If death_model is specified, parse its parameters
    this->use_metabolic_death = true;

    // Parse death_type (apoptosis or necrosis)
    pugi::xml_node node = parent.child("death_type");
    if (node){
        this->death_type = PhysiCell::xml_get_my_string_value(node);
        assert(this->death_type == "apoptosis" || this->death_type == "necrosis");
    } else {
        std::cout << "Warning: death_type not specified. Defaulting to apoptosis." << std::endl;
        this->death_type = "apoptosis";
    }

    // Parse monitored_flux (optional)
    node = parent.child("death_trigger_flux");
    if (node){
        this->death_trigger_flux = PhysiCell::xml_get_my_string_value(node);
    } else {
        std::cout << "Warning: monitored_flux not specified. Metabolic-dependent death disabled." << std::endl;
        this->death_trigger_flux = "";
        return;
    }

    // Parse flux_threshold (optional, default provided)
    node = parent.child("death_flux_threshold");
    if (node){
        this->death_flux_threshold = PhysiCell::xml_get_my_double_value(node);
    } else {
        std::cout << "Warning: flux_threshold not specified. Using default 1e-6." << std::endl;
        this->death_flux_threshold = 1e-6;
    }

    // Parse death_rate (optional, default provided)
    node = parent.child("death_rate_increase");
    if (node){
        this->death_rate_increase = PhysiCell::xml_get_my_double_value(node);
    } else {
        std::cout << "Warning: death_rate not specified. Using default 0.01." << std::endl;
        this->death_rate_increase = 0.01;
    }

    std::cout << "Metabolic-dependent death enabled: " << this->death_type << ", monitored flux: " << this->death_trigger_flux << ", threshold: " << this->death_flux_threshold << ", rate: " << this->death_rate_increase << std::endl;
}

void dFBAIntracellular::initialize_intracellular_from_pugixml(pugi::xml_node& node)
{

    // Setting all the rest to default values : Nothing should be kept from the existing intracellular object (NO INHERITANCE)
    objective_reaction = "";
    cell_density = 0.0;
    max_growth_rate = 0.0;
    current_growth_rate = 0.0;
    last_solution_status = "not_run";
    use_metabolic_death = false;
	death_type = "";
	death_trigger_flux = "";
	death_flux_threshold = 0.0;
	death_rate_increase = 0.0;
    substrate_exchanges.clear();
	flag_for_death = false;
    sbml_model.clear();
    is_initialized = false;

    std::cout << "===================================================" << std::endl;
    std::cout << "Parssing dFBA intracellular model" << std::endl;
    
    // Parse dFBA sbml_filename and time step from XML settings
    pugi::xml_node node_settings = node.child("settings");
    if (node_settings) 
    {
        pugi::xml_node node_intracellular_dt = node_settings.child("intracellular_dt");
        pugi::xml_node node_time_step = node_settings.child("time_step");
        if (node_intracellular_dt) 
        {
            dfba_time_step = PhysiCell::xml_get_my_double_value(node_intracellular_dt);
        } else if (node_time_step) 
        {
            std::cout << "[Warning] The setting 'time_step' is deprecated. Please use 'intracellular_dt' instead." << std::endl;
            dfba_time_step = PhysiCell::xml_get_my_double_value(node_time_step);
        } else 
        {
            std::cout << "[Warning] No intracellular_dt or time_step specified. Using default value: " << PhysiCell::diffusion_dt << std::endl;
            dfba_time_step = PhysiCell::diffusion_dt; // Default value
        }
        // Getting sbml file name for reading the model
        pugi::xml_node node_sbml = node_settings.child( "sbml_filename" );
        if ( node_sbml ) 
        { 
            
            this->sbml_filename = PhysiCell::xml_get_my_string_value (node_sbml);
        }
        else
        {
            std::cout << "Error: attempted to read sbml_filename path but not found." << std::endl;
            std::cout << PhysiCell::xml_get_my_string_value(node_sbml) << std::endl;
            std::cout << "Please double-check your <intracellular> node in the XML setting." << std::endl;
            std::cout << std::endl; 
            exit(-1); 
        }
    } else 
    {
        std::cout << "Error: No settings block found in the XML" << std::endl;
        std::cout << "Please double-check your intracellular node in the XML setting." << std::endl;
        std::cout << "intracellular.settings node in the XML setting should include <sbml_filename> and <intracellular_dt>" << std::endl;
        std::cout << std::endl; 
        exit(-1); 
    }

    // parsing the transport model
    pugi::xml_node node_transport_model = node.child( "transport_model" );
    if ( node_transport_model )
	{ 
        int num_exchanges = parse_transport_model(node_transport_model);
        if (num_exchanges == 0){
            std::cout << "Error, dFBA model must have at least one exchange flux." << std::endl;
            std::cout << "Please double-check your exchange nodes in the XML setting file." << std::endl;
            std::cout << std::endl; 
            exit(-1);
        }
    }
    else
    {
        std::cout << "Error: attempted to parse transport_model but not found." << std::endl;
        std::cout << "Please double-check intracellular model nodes in the XML setting file." << std::endl;
        std::cout << std::endl; 
        exit(-1); 
    }


    // parsing the transport model
    pugi::xml_node node_growth_model = node.child( "growth_model" );
    if ( node_growth_model )
	{ 
        parse_growth_model(node_growth_model);
    }
    else
    {
        std::cout << "Error: attempted to parse growth_model but not found." << std::endl;
        std::cout << "Please double-check intracellular model nodes in the XML setting file." << std::endl;
        std::cout << std::endl; 
        exit(-1); 
    }

    // parsing the death model
    pugi::xml_node node_death_model = node.child( "death_model" ); 
    if ( node_death_model )
	{ 
        bool death_enabled = node_death_model.attribute("enabled").as_bool();
        if (!death_enabled){
            std::cout << "Death model disabled. Using default behavior (no metabolic-dependent death)." << std::endl;
            this->use_metabolic_death = false;
        }
        else{
        parse_death_model(node_death_model);
        }
    }
    else
    {
        std::cout << "Death model not specified. Using default behavior (no metabolic-dependent death)." << std::endl;
        this->use_metabolic_death = false;
    }


    // ==============================================================================================

    std::cout << "Loading SBML model from: " << this->sbml_filename << std::endl;
    this->sbml_model.initModel(this->sbml_filename.c_str());
    
    map<std::string, ExchangeFluxData>::iterator it;
    for(it = this->substrate_exchanges.begin(); it != this->substrate_exchanges.end(); it++)
    {
        std::string substrate_name = it->first;
        ExchangeFluxData ex_strut = it->second;
        // TODO
        // check this ex_strut.density_index is a defined density at Microenviroment
        dFBAReaction* rxn = this->sbml_model.getReaction(ex_strut.fba_flux_id);
        if( rxn == nullptr ){
            std::cout << "Error: exchange reaction not found in model: " << ex_strut.fba_flux_id << " for substrate " << substrate_name << std::endl;
            exit(-1);
        }
    }
    dFBAReaction* growth_rxn = this->sbml_model.getReaction(this->objective_reaction);
    assert( growth_rxn != nullptr );
    this->sbml_model.setReactionLowerBound(this->objective_reaction, 0.0);
    this->sbml_model.setReactionUpperBound(this->objective_reaction, this->max_growth_rate);

    if(this->use_metabolic_death && !this->death_trigger_flux.empty()){
        dFBAReaction* death_rxn = this->sbml_model.getReaction(this->death_trigger_flux);
        if(death_rxn == nullptr){
            std::cout << "ERROR: Specified death_trigger_flux (" << this->death_trigger_flux << ") does not exist in the SBML model." << std::endl;
            exit(-1);
        }
        else{
            this->sbml_model.setReactionLowerBound(this->death_trigger_flux, this->death_flux_threshold);
        }
    }

    std::cout << "Done!" << std::endl;
    std::cout << "===================================================" << std::endl;
}

void dFBAIntracellular::start()
{
    this->next_dfba_run = PhysiCell::PhysiCell_globals.current_time + dfba_time_step;
    /*
    for (const auto& reaction : this->sbml_model.getListOfReactions())
    {   
    if (reaction != nullptr) {
        std::cout << "Reaction: " << this->sbml_model.getReactionIndex(reaction->getId()) 
                  << " Reaction id: " << reaction->getId() 
                  << " Reaction string: " << reaction->getReactionString(this->sbml_model) << std::endl;
        }
    }
    */

}


void dFBAIntracellular::update_dfba_inputs( PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt ){

    // HeLa cell mass 2.3 ng
    // cell volume fL (um³)
    // mM: mmol / L = 10⁻³ mol / 1e¹⁵ um³ = 10⁻¹⁵mol / um³ = amol /  um³
    // cell_density --> 1.04 g/ml = 1.04 ug/nL = 1.04 ng / pL = 1.04 pg / fL = pg / um³
    // cell.mass.total (ng) = cell.volume (um³) * cell.density (pg / um³) = pg
    // cell_density 1.04 pg/um³
    // cell volume ~ 2500 (um³) * ~ 1.04 pg / (um³) 
    // gDW cell (cell.volume.total * cell_density) = mass.total 
    // 2500 (um³) * 1.04 (pg/um³) =~ 2600 (pg) ~= 2.3 ng
    // cell.mass.solid = cell.mass.total * (1-fluid_frac) = 780 pg
 
    // Cell Volume
    double current_volume = phenotype.volume.total;
    // unitless ~30%
    double solid_fraction = 1 - phenotype.volume.fluid_fraction;     
    // um³ = um³
    double solid_volume = current_volume * solid_fraction;   
    // um³ * pg / um³ = pg
    double cell_dry_weight = solid_volume * this->cell_density ; 
    // pg * 1e-12 = g
    double mass_scaling = cell_dry_weight * 1e-12;
    
    // r = (3V / 4π))^1/3 (um)
    double radius = cbrt((3.0 * current_volume) / (4.0 * PI)); 
    // um * um = um²
    double cell_surface = 4 * PI * pow(radius, 2);

    // Setp 1 - Update exchange fluxes lower
    //  Standard FBA fluxes units:  mmol /  gDW cell / h
    //  Km: mM = mmol/L
    //  Vmax: mmol/g DW cell/min

    int voxel_index = pCell->get_current_voxel_index();
    double dV = microenvironment.voxels(voxel_index).volume;
    double local_cell_count = 1.0;
    if( microenvironment.agent_container && voxel_index >= 0 && voxel_index < (int) microenvironment.agent_container->agent_grid.size() )
    {
        local_cell_count = std::max(1.0, (double) microenvironment.agent_container->agent_grid[voxel_index].size());
    }
    const double substrate_availability_safety = 0.05;

    std::vector<double> density_vector = pCell->nearest_density_vector(); 
    for( auto& clamp_it : this->substrate_exchanges )
    {
        int clamp_index = clamp_it.second.density_index;
        if( density_vector[clamp_index] < 0.0 )
        {
            density_vector[clamp_index] = 0.0;
        }
    }

    // Only declare the iterator ONCE per function
    map<std::string, ExchangeFluxData>::iterator it;
    for(it = this->substrate_exchanges.begin(); it != this->substrate_exchanges.end(); it++)
    {
        
        std::string substrate_name = it->first;
        ExchangeFluxData ex_strut = it->second;

        // geting the amount of substrate
        double substrate_conc = density_vector[ex_strut.density_index];
        substrate_conc = max(substrate_conc, 0.0); // mM 
        // scaling Vmax based on cell volume
        double Vmax = ex_strut.Vmax.value; // mmol / gDWcell / hours
        double Km   = ex_strut.Km.value; // mM 
        
        // using irreversible Michaelis Menten kinetics to estimate the flux bound; should be calculated from density
        double uptake_rate = (Vmax * substrate_conc) / (Km + substrate_conc); //  mmol / gDWcell / hours shoudl 
        

        // Here we estimate the amount BioFVM can safely remove over this dFBA step.
        // Each cell in the same voxel only receives a share of the local substrate pool.
        double total_substrate = substrate_conc * dV / liter_micron_cubes_conversion; // mmol
        double available_substrate_per_cell = substrate_availability_safety * total_substrate / local_cell_count;
        double dfba_dt_hours = std::max(this->dfba_time_step * hours_to_minutes, 1e-12);
        double availability_limited_uptake_rate = 0.0;
        if( mass_scaling > 0.0 )
        {
            availability_limited_uptake_rate = available_substrate_per_cell / (mass_scaling * dfba_dt_hours); // mmol/gDW/h
        }
        if( ex_strut.allow_uptake )
        {
            uptake_rate = std::min(uptake_rate, availability_limited_uptake_rate);
        }
        if( substrate_conc <= 1e-12 )
        {
            uptake_rate = 0.0;
        }

        double total_uptake = uptake_rate * dfba_dt_hours * mass_scaling; // mmol

        std::ostringstream substrate_log;
        substrate_log << "Substrate: " << substrate_name << " -- Total substrate in voxel: " << total_substrate << " mmol," << "Total concentration: " << substrate_conc << " mM, Total uptake requested: " << total_uptake << " mmol. Uptake rate: " << uptake_rate << " mmol/gDW/h\n";

        // Change sign to use as lower bound of the exchange flux.
        // Species marked export-only must not become extracellular nutrients.
        double exchange_flux_lb = ex_strut.allow_uptake ? -1 * uptake_rate : 0.0;
        bool force_minimum_uptake = ex_strut.allow_uptake && ex_strut.minimum_uptake_fraction > 0.0 && uptake_rate > 1e-12;
        double exchange_flux_ub = force_minimum_uptake ? -1.0 * ex_strut.minimum_uptake_fraction * uptake_rate : 1000.0;
        if( exchange_flux_lb > exchange_flux_ub )
        {
            exchange_flux_lb = exchange_flux_ub;
        }
        substrate_log << "Substrate: " << substrate_name << " Density: " << substrate_conc << " Vmax: " << Vmax << " Km: " << Km << " Max Uptake rate (mmol/gDW/h): " << uptake_rate <<  " Exchange flux LB: " << exchange_flux_lb << " Exchange flux UB: " << exchange_flux_ub << "\n";
        // Uncomment for detailed per-cell/per-substrate diagnostics.
        // {
        //     static std::mutex substrate_log_mutex;
        //     std::lock_guard<std::mutex> lock(substrate_log_mutex);
        //     std::cout << substrate_log.str() << std::flush;
        // }
        // Updateing the lower bound of the corresponding exchange flux
        this->sbml_model.setReactionLowerBound(ex_strut.fba_flux_id, exchange_flux_lb);
        this->sbml_model.setReactionUpperBound(ex_strut.fba_flux_id, exchange_flux_ub);

    }
}

void dFBAIntracellular::update(){
    // Only run dFBA if current_time >= next_dfba_run
    dFBASolution solution = this->sbml_model.optimize();
    this->last_solution_status = solution.status;

    // Detailed per-cell dFBA logging is disabled by default; use fba_exchanges.csv
    // for flux lower bounds and optimized fluxes.

    //next_dfba_run = PhysiCell::PhysiCell_globals.current_time + dfba_time_step;
    if (solution.status == "infeasible"){
        //std::cout << "I'm dead from the metabolic point of view" << std::endl;
        this->flag_for_death = true;
        this->current_growth_rate = 0.0;
    }
    else if(solution.status == "unknown"){
        std::cout << "ERROR: Unknown status for the FBA problem!" << std::endl;
        exit(1);
    }
    else{
        this->current_growth_rate = solution.getObjectiveValue();
        this->flag_for_death = false;
    }
}

void dFBAIntracellular::save_fluxes_to_csv(PhysiCell::Cell* pCell, double current_time, std::string output_folder)
{
    static std::mutex save_fluxes_mutex;
    std::lock_guard<std::mutex> lock(save_fluxes_mutex);
    static bool header_summary = false;
    std::string path_summary = output_folder + "/fba_summary.csv";
    std::ofstream f_summary;

    if (!header_summary) {
        f_summary.open(path_summary, std::ios::out);
        f_summary << "time,cell_id,cell_type,x,y,dfba_status,growth_rate,obj_flux";
        for (auto& it : this->substrate_exchanges)
            f_summary << "," << it.second.density_name << "_flux";
        f_summary << "\n";
        header_summary = true;
    } else {
        f_summary.open(path_summary, std::ios::app);
    }

    f_summary << current_time 
              << "," << pCell->ID 
              << "," << pCell->type_name
              << "," << pCell->position[0] 
              << "," << pCell->position[1]
              << "," << this->last_solution_status
              << "," << this->current_growth_rate;

    dFBAReaction* obj = this->sbml_model.getReaction(this->objective_reaction);
    f_summary << "," << (obj ? obj->getFluxValue() : 0.0);


    std::vector<double> density_vector = pCell->nearest_density_vector();
    for( auto& clamp_it : this->substrate_exchanges )
    {
        int clamp_index = clamp_it.second.density_index;
        if( density_vector[clamp_index] < 0.0 )
        {
            density_vector[clamp_index] = 0.0;
        }
    }
    double glucose_flux_for_co2 = 0.0;
    auto glucose_exchange_for_co2 = this->substrate_exchanges.find("glucose");
    if( glucose_exchange_for_co2 != this->substrate_exchanges.end() )
    {
        dFBAReaction* glucose_rxn_for_co2 = this->sbml_model.getReaction(glucose_exchange_for_co2->second.fba_flux_id);
        if( glucose_rxn_for_co2 )
        {
            glucose_flux_for_co2 = glucose_rxn_for_co2->getFluxValue();
        }
    }
    double oxygen_flux_for_co2 = 0.0;
    double oxygen_conc_for_lactate = 0.0;
    auto oxygen_exchange_for_lactate = this->substrate_exchanges.find("oxygen");
    if( oxygen_exchange_for_lactate != this->substrate_exchanges.end() )
    {
        oxygen_conc_for_lactate = std::max(density_vector[oxygen_exchange_for_lactate->second.density_index], 0.0);
        dFBAReaction* oxygen_rxn_for_co2 = this->sbml_model.getReaction(oxygen_exchange_for_lactate->second.fba_flux_id);
        if( oxygen_rxn_for_co2 )
        {
            oxygen_flux_for_co2 = oxygen_rxn_for_co2->getFluxValue();
        }
    }
    double lactate_conc_for_exchange = 0.0;
    auto lactate_exchange_for_exchange = this->substrate_exchanges.find("lactate");
    if( lactate_exchange_for_exchange != this->substrate_exchanges.end() )
    {
        lactate_conc_for_exchange = std::max(density_vector[lactate_exchange_for_exchange->second.density_index], 0.0);
    }
    const double lactate_hypoxia_threshold = 0.06; // mM, experimental hypoxic lactate test threshold
    const double lactate_anoxia_threshold = 0.02; // mM
    double hypoxia_fraction_for_lactate = 0.0;
    if( oxygen_conc_for_lactate < lactate_hypoxia_threshold )
    {
        hypoxia_fraction_for_lactate = (lactate_hypoxia_threshold - oxygen_conc_for_lactate) / (lactate_hypoxia_threshold - lactate_anoxia_threshold);
        hypoxia_fraction_for_lactate = std::min(1.0, std::max(0.0, hypoxia_fraction_for_lactate));
    }
    const double lactate_oxidation_o2_threshold = 0.06; // mM
    double oxygen_fraction_for_lactate_uptake = 0.0;
    if( oxygen_conc_for_lactate > lactate_hypoxia_threshold )
    {
        oxygen_fraction_for_lactate_uptake = (oxygen_conc_for_lactate - lactate_hypoxia_threshold) / (lactate_oxidation_o2_threshold - lactate_hypoxia_threshold);
        oxygen_fraction_for_lactate_uptake = std::min(1.0, std::max(0.0, oxygen_fraction_for_lactate_uptake));
    }
    double lactate_zone_uptake_factor = 0.25;
    if( pCell->type_name == "zone_1" )
    {
        lactate_zone_uptake_factor = 1.0;
    }
    else if( pCell->type_name == "zone_2" )
    {
        lactate_zone_uptake_factor = 0.5;
    }
    else if( pCell->type_name == "zone_3" )
    {
        lactate_zone_uptake_factor = 0.1;
    }
    for (auto& it : this->substrate_exchanges) {
        dFBAReaction* rxn = this->sbml_model.getReaction(it.second.fba_flux_id);
        double flux_value = rxn ? rxn->getFluxValue() : 0.0;
        if( it.second.density_name == "CO2" && glucose_flux_for_co2 < 0.0 )
        {
            // Oxidative CO2 cannot exceed either glucose-derived carbon or O2-supported oxidation.
            double co2_from_glucose = 6.0 * (-glucose_flux_for_co2);
            double co2_from_oxygen = oxygen_flux_for_co2 < 0.0 ? -oxygen_flux_for_co2 : co2_from_glucose;
            double max_oxidative_co2 = std::max(0.0, std::min(co2_from_glucose, co2_from_oxygen));
            if( flux_value <= 1e-12 )
            {
                flux_value = max_oxidative_co2;
            }
            else
            {
                flux_value = std::min(flux_value, max_oxidative_co2);
            }
        }
        if( it.second.density_name == "lactate" )
        {
            double lactate_release_zone_factor = 0.1;
            if( pCell->type_name == "zone_1" )
            {
                lactate_release_zone_factor = 0.0;
            }
            else if( pCell->type_name == "zone_2" )
            {
                lactate_release_zone_factor = 0.25;
            }
            else if( pCell->type_name == "zone_3" )
            {
                lactate_release_zone_factor = 1.0;
            }

            double max_glycolytic_lactate = glucose_flux_for_co2 < 0.0 ? 2.0 * (-glucose_flux_for_co2) : 0.0;
            double hypoxic_lactate_release = max_glycolytic_lactate
                * hypoxia_fraction_for_lactate
                * lactate_release_zone_factor;

            if( flux_value > 0.0 )
            {
                flux_value = std::min(flux_value, hypoxic_lactate_release);
            }
            else if( std::abs(flux_value) <= 1e-12 && hypoxic_lactate_release > 0.0 )
            {
                flux_value = hypoxic_lactate_release;
            }
            else if( std::abs(flux_value) <= 1e-12 && lactate_conc_for_exchange > 1e-12 && oxygen_fraction_for_lactate_uptake > 0.0 )
            {
                double lactate_Km = std::max(static_cast<double>(it.second.Km.value), 1e-12);
                double lactate_Vmax = std::max(static_cast<double>(it.second.Vmax.value), 1.0);
                flux_value = -lactate_Vmax * lactate_conc_for_exchange / (lactate_Km + lactate_conc_for_exchange) * oxygen_fraction_for_lactate_uptake * lactate_zone_uptake_factor;
            }
        }
        f_summary << "," << flux_value;
    }
    f_summary << "\n";
    f_summary.close();

    static bool header_exchanges = false;
    std::string path_exchanges = output_folder + "/fba_exchanges.csv";
    std::ofstream f_exchanges;

    if (!header_exchanges) {
        f_exchanges.open(path_exchanges, std::ios::out);
        f_exchanges << "time,cell_id,cell_type,dfba_status,substrate,fba_reaction,flux_lb,flux_value,concentration_mM\n";
        header_exchanges = true;
    } else {
        f_exchanges.open(path_exchanges, std::ios::app);
    }

    for (auto& it : this->substrate_exchanges) {
        ExchangeFluxData ex = it.second;
        dFBAReaction* rxn = this->sbml_model.getReaction(ex.fba_flux_id);
        double conc = std::max(density_vector[ex.density_index], 0.0);
        double flux_value = rxn ? rxn->getFluxValue() : 0.0;
        if( ex.density_name == "CO2" && glucose_flux_for_co2 < 0.0 )
        {
            double co2_from_glucose = 6.0 * (-glucose_flux_for_co2);
            double co2_from_oxygen = oxygen_flux_for_co2 < 0.0 ? -oxygen_flux_for_co2 : co2_from_glucose;
            double max_oxidative_co2 = std::max(0.0, std::min(co2_from_glucose, co2_from_oxygen));
            if( flux_value <= 1e-12 )
            {
                flux_value = max_oxidative_co2;
            }
            else
            {
                flux_value = std::min(flux_value, max_oxidative_co2);
            }
        }
        if( ex.density_name == "lactate" )
        {
            double lactate_release_zone_factor = 0.1;
            if( pCell->type_name == "zone_1" )
            {
                lactate_release_zone_factor = 0.0;
            }
            else if( pCell->type_name == "zone_2" )
            {
                lactate_release_zone_factor = 0.25;
            }
            else if( pCell->type_name == "zone_3" )
            {
                lactate_release_zone_factor = 1.0;
            }

            double max_glycolytic_lactate = glucose_flux_for_co2 < 0.0 ? 2.0 * (-glucose_flux_for_co2) : 0.0;
            double hypoxic_lactate_release = max_glycolytic_lactate
                * hypoxia_fraction_for_lactate
                * lactate_release_zone_factor;

            if( flux_value > 0.0 )
            {
                flux_value = std::min(flux_value, hypoxic_lactate_release);
            }
            else if( std::abs(flux_value) <= 1e-12 && hypoxic_lactate_release > 0.0 )
            {
                flux_value = hypoxic_lactate_release;
            }
            else if( std::abs(flux_value) <= 1e-12 && lactate_conc_for_exchange > 1e-12 && oxygen_fraction_for_lactate_uptake > 0.0 )
            {
                double lactate_Km = std::max(static_cast<double>(ex.Km.value), 1e-12);
                double lactate_Vmax = std::max(static_cast<double>(ex.Vmax.value), 1.0);
                flux_value = -lactate_Vmax * lactate_conc_for_exchange / (lactate_Km + lactate_conc_for_exchange) * oxygen_fraction_for_lactate_uptake * lactate_zone_uptake_factor;
            }
        }
        f_exchanges << current_time
                    << "," << pCell->ID
                    << "," << pCell->type_name
                    << "," << this->last_solution_status
                    << "," << ex.density_name
                    << "," << ex.fba_flux_id
                    << "," << (rxn ? rxn->getLowerBound() : 0.0)
                    << "," << flux_value
                    << "," << conc
                    << "\n";
    }
    f_exchanges.close();
}

void dFBAIntracellular::update_dfba_outputs(PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt )
{
    // Metabolic-dependent death check
    //std::cout << "Cell ready to die? --> " << this->flag_for_death << std::endl;
    //std::cout << "Death parameters: " << this->use_metabolic_death << " " << this->death_trigger_flux << " " << this->death_flux_threshold << " " << this->death_rate_increase << std::endl;
    static int nApoptosis = phenotype.death.find_death_model_index(PhysiCell::PhysiCell_constants::apoptosis_death_model );
    static int nNecrosis = phenotype.death.find_death_model_index(PhysiCell::PhysiCell_constants::necrosis_death_model );

    if(this->flag_for_death && this->use_metabolic_death)
    {
        if(this->death_type == "apoptosis")
        {
            pCell->phenotype.death.rates[nApoptosis] += this->death_rate_increase;
        }
        else if(this->death_type == "necrosis")
        {
            pCell->phenotype.death.rates[nNecrosis] += this->death_rate_increase;
        }
    } 
    else if(this->use_metabolic_death && !this->flag_for_death)
    {
        if(this->death_type == "apoptosis")
        {
            pCell->phenotype.death.rates[nApoptosis] = 0.0;
        }
        else if(this->death_type == "necrosis")
        {
            pCell->phenotype.death.rates[nNecrosis] = 0.0;
        }
    }

    double fba_growth_rate = this->current_growth_rate;
    
    
    if (fba_growth_rate > fba_epsilon_tolerance) {
        // SCALE GROWTH RATE BY this->dfba_time_step
        double growth_rate = fba_growth_rate * this->dfba_time_step * hours_to_minutes; // growth_rate 1/h * dt (h)

        // exact solution (fallback to linear for very small x to avoid exp overhead)
        double factor = (std::abs(growth_rate) < 1e-6) ? (1.0 + growth_rate) : std::exp(growth_rate);
        double safe_factor = std::max(factor, 1e-12);
        //double volume_increase_ratio = 1.0 + growth_rate;
        phenotype.volume.multiply_by_ratio( safe_factor );
        pCell->set_total_volume( phenotype.volume.total );
        phenotype.geometry.update(pCell, phenotype, this->dfba_time_step);
    }

    
    // Cell Volume
    double current_volume = phenotype.volume.total;
    // unitless ~30%
    double solid_fraction = 1.0 - phenotype.volume.fluid_fraction;     
    // um³ = um³
    double solid_volume = current_volume * solid_fraction;   
    // um³ * pg / um³ = pg
    double cell_dry_weight = solid_volume * this->cell_density;
    
    // pg * 1e-12 = g
    // pCell->custom_data["cell_dry_weight"] = cell_dry_weight;
    cell_dry_weight *= 1e-12;


    
	
    float current_growth_rate = this->get_growth_rate();
    // std::cout << "Current growth rate: " << current_growth_rate << std::endl;
    // pCell->custom_data["growth_rate"] = current_growth_rate;
    // float R_biomass_reaction = this->get_flux_value(this->objective_reaction);
    // float biomass_ub = this->sbml_model.getReactionUpperBound("R_biomass_reaction");
    // float biomass_lb = this->sbml_model.getReactionLowerBound("R_biomass_reaction");
    // std::cout << "Biomass reaction flux value: " << R_biomass_reaction << endl;
	// pCell->custom_data["R_biomass_reaction"] = R_biomass_reaction;
    
    


    std::vector<double> density_vector = pCell->nearest_density_vector(); 
    for( auto& clamp_it : this->substrate_exchanges )
    {
        int clamp_index = clamp_it.second.density_index;
        if( density_vector[clamp_index] < 0.0 )
        {
            density_vector[clamp_index] = 0.0;
        }
    }
    int voxel_index = pCell->get_current_voxel_index();
    double dV = microenvironment.voxels(voxel_index).volume;
    double local_cell_count = 1.0;
    if( microenvironment.agent_container && voxel_index >= 0 && voxel_index < (int) microenvironment.agent_container->agent_grid.size() )
    {
        local_cell_count = std::max(1.0, (double) microenvironment.agent_container->agent_grid[voxel_index].size());
    }
    const double substrate_availability_safety = 0.05;
    // For each substrate exchange, scale net_export_rates by this->dfba_time_step
    double glucose_flux_for_co2_output = 0.0;
    auto glucose_exchange_for_co2_output = this->substrate_exchanges.find("glucose");
    if( glucose_exchange_for_co2_output != this->substrate_exchanges.end() )
    {
        dFBAReaction* glucose_rxn_for_co2_output = this->sbml_model.getReaction(glucose_exchange_for_co2_output->second.fba_flux_id);
        if( glucose_rxn_for_co2_output )
        {
            glucose_flux_for_co2_output = glucose_rxn_for_co2_output->getFluxValue();
        }
    }
    double oxygen_flux_for_co2_output = 0.0;
    double oxygen_conc_for_lactate_output = 0.0;
    auto oxygen_exchange_for_lactate_output = this->substrate_exchanges.find("oxygen");
    if( oxygen_exchange_for_lactate_output != this->substrate_exchanges.end() )
    {
        oxygen_conc_for_lactate_output = std::max(density_vector[oxygen_exchange_for_lactate_output->second.density_index], 0.0);
        dFBAReaction* oxygen_rxn_for_co2_output = this->sbml_model.getReaction(oxygen_exchange_for_lactate_output->second.fba_flux_id);
        if( oxygen_rxn_for_co2_output )
        {
            oxygen_flux_for_co2_output = oxygen_rxn_for_co2_output->getFluxValue();
        }
    }
    double lactate_conc_for_exchange_output = 0.0;
    auto lactate_exchange_for_exchange_output = this->substrate_exchanges.find("lactate");
    if( lactate_exchange_for_exchange_output != this->substrate_exchanges.end() )
    {
        lactate_conc_for_exchange_output = std::max(density_vector[lactate_exchange_for_exchange_output->second.density_index], 0.0);
    }
    const double lactate_hypoxia_threshold_output = 0.06; // mM, experimental hypoxic lactate test threshold
    const double lactate_anoxia_threshold_output = 0.02; // mM
    double hypoxia_fraction_for_lactate_output = 0.0;
    if( oxygen_conc_for_lactate_output < lactate_hypoxia_threshold_output )
    {
        hypoxia_fraction_for_lactate_output = (lactate_hypoxia_threshold_output - oxygen_conc_for_lactate_output) / (lactate_hypoxia_threshold_output - lactate_anoxia_threshold_output);
        hypoxia_fraction_for_lactate_output = std::min(1.0, std::max(0.0, hypoxia_fraction_for_lactate_output));
    }
    const double lactate_oxidation_o2_threshold_output = 0.06; // mM
    double oxygen_fraction_for_lactate_uptake_output = 0.0;
    if( oxygen_conc_for_lactate_output > lactate_hypoxia_threshold_output )
    {
        oxygen_fraction_for_lactate_uptake_output = (oxygen_conc_for_lactate_output - lactate_hypoxia_threshold_output) / (lactate_oxidation_o2_threshold_output - lactate_hypoxia_threshold_output);
        oxygen_fraction_for_lactate_uptake_output = std::min(1.0, std::max(0.0, oxygen_fraction_for_lactate_uptake_output));
    }
    double lactate_zone_uptake_factor_output = 0.25;
    if( pCell->type_name == "zone_1" )
    {
        lactate_zone_uptake_factor_output = 1.0;
    }
    else if( pCell->type_name == "zone_2" )
    {
        lactate_zone_uptake_factor_output = 0.5;
    }
    else if( pCell->type_name == "zone_3" )
    {
        lactate_zone_uptake_factor_output = 0.1;
    }
    map<std::string, ExchangeFluxData>::iterator it;
    for(it = this->substrate_exchanges.begin(); it != this->substrate_exchanges.end(); it++)
    {
        // Retrive the exchange flux and its corresponding density
        std::string substrate_name = it->first;
        ExchangeFluxData ex_strut = it->second;

        double Vmax = ex_strut.Vmax.value;
        double Km   = ex_strut.Km.value;
        
        int density_index = ex_strut.density_index;
        std::string density_name = ex_strut.density_name;
        std::string fba_flux_id = ex_strut.fba_flux_id;

        double substrate_conc = density_vector[ex_strut.density_index]; // mM = mmol/L
        double total_substrate = substrate_conc * (dV / 1e15); // mmol
        
        // pCell->custom_data[density_name + "_conc"] = substrate_conc;
		// pCell->custom_data[density_name + "_total"] = total_substrate;


        dFBAReaction* exchange_flux = this->sbml_model.getReaction(fba_flux_id);
        double flux_value =  exchange_flux->getFluxValue(); // mmol/gDW/h
        if( !ex_strut.allow_uptake && flux_value < 0.0 )
        {
            flux_value = 0.0;
        }
        if( density_name == "CO2" && glucose_flux_for_co2_output < 0.0 )
        {
            // Oxidative CO2 cannot exceed either glucose-derived carbon or O2-supported oxidation.
            double co2_from_glucose = 6.0 * (-glucose_flux_for_co2_output);
            double co2_from_oxygen = oxygen_flux_for_co2_output < 0.0 ? -oxygen_flux_for_co2_output : co2_from_glucose;
            double max_oxidative_co2 = std::max(0.0, std::min(co2_from_glucose, co2_from_oxygen));
            if( flux_value <= 1e-12 )
            {
                flux_value = max_oxidative_co2;
            }
            else
            {
                flux_value = std::min(flux_value, max_oxidative_co2);
            }
        }
        if( density_name == "lactate" )
        {
            double lactate_release_zone_factor_output = 0.1;
            if( pCell->type_name == "zone_1" )
            {
                lactate_release_zone_factor_output = 0.0;
            }
            else if( pCell->type_name == "zone_2" )
            {
                lactate_release_zone_factor_output = 0.25;
            }
            else if( pCell->type_name == "zone_3" )
            {
                lactate_release_zone_factor_output = 1.0;
            }

            double max_glycolytic_lactate_output = glucose_flux_for_co2_output < 0.0 ? 2.0 * (-glucose_flux_for_co2_output) : 0.0;
            double hypoxic_lactate_release_output = max_glycolytic_lactate_output
                * hypoxia_fraction_for_lactate_output
                * lactate_release_zone_factor_output;

            if( flux_value > 0.0 )
            {
                flux_value = std::min(flux_value, hypoxic_lactate_release_output);
            }
            else if( std::abs(flux_value) <= 1e-12 && hypoxic_lactate_release_output > 0.0 )
            {
                flux_value = hypoxic_lactate_release_output;
            }
            else if( std::abs(flux_value) <= 1e-12 && lactate_conc_for_exchange_output > 1e-12 && oxygen_fraction_for_lactate_uptake_output > 0.0 )
            {
                double lactate_Km = std::max(static_cast<double>(ex_strut.Km.value), 1e-12);
                double lactate_Vmax = std::max(static_cast<double>(ex_strut.Vmax.value), 1.0);
                flux_value = -lactate_Vmax * lactate_conc_for_exchange_output / (lactate_Km + lactate_conc_for_exchange_output) * oxygen_fraction_for_lactate_uptake_output * lactate_zone_uptake_factor_output;
            }
        }
        // std::cout << "Exchange flux : " << fba_flux_id << " Flux value (mmol/gDW/h): " << flux_value << std::endl;
        
        // pCell->custom_data[fba_flux_id] = flux_value;

        // Rescaling FBA exchanges flux into net_export_rates
        // Net export rates are expressed in substance/time
        // flux_value: mmol/gDW/h --> mmol/min
        // flux_value: mmol/gDW/h; cell_dry_weight was already converted from pg to g above.
        // net_export_rate (mmol/min) = flux_value / 60 * cell_dry_weight_gDW.
        //double net_export_rate = flux_value * cell_dry_weight * hours_to_minutes; // mmol/min
        double cell_dry_weight_g = cell_dry_weight;
        double net_export_rate_mmol_per_min = flux_value * cell_dry_weight_g * hours_to_minutes; // mmol/min
        
 
        /*
        double substrate_conc = density_vector[ex_strut.density_index]; // mM = mmol/L
        std::cout << " Cell Type: " << pCell->type_name << std::endl;
        std::cout << "\tSubstrate: " << substrate_name << std::endl;
        std::cout << "\tconcentration: " << substrate_conc << std::endl;
        std::cout << "\tKm: " << Km << std::endl;
        std::cout << "\tVmax: " << Vmax << std::endl;
        std::cout << "\tTotal substrate: " << total_substrate << std::endl;
        std::cout << "\tConsumption: " << substrate_consumption << std::endl;
        std::cout << "\tNet export rate: " << net_export_rate << std::endl;
        std::cout << "\tFBA Flux " << flux_value << std::endl;
        std::cout << "\tFBA growth rate " << fba_growth_rate << std::endl;
        if (substrate_name == "lactate" && flux_value > 1.0){
       }
        */

        // correct scaling that takes into account the volume units (liter to um3) in BioFVM 
        net_export_rate_mmol_per_min *= liter_micron_cubes_conversion; // BioFVM units: mM*um3/min

        if( net_export_rate_mmol_per_min < 0.0 )
        {
            double safe_substrate_conc = std::max(substrate_conc, 0.0);
            double max_removal_rate = substrate_availability_safety * safe_substrate_conc * dV / std::max(this->dfba_time_step, 1e-12) / local_cell_count;
            net_export_rate_mmol_per_min = std::max(net_export_rate_mmol_per_min, -max_removal_rate);
        }

        // std::cout << "New net export rate for substrate " << substrate_name << ": " << net_export_rate_mmol_per_min << " mmol/min" << std::endl;
        phenotype.secretion.net_export_rates[density_index] = net_export_rate_mmol_per_min;

        // phenotype.secretion.net_export_rates[density_index] = net_export_rate;
        pCell->set_internal_uptake_constants(this->dfba_time_step); // Update internal uptake constants based on the new volume and rates
        if (default_microenvironment_options.track_internalized_substrates_in_each_agent)
        {
            phenotype.molecular.internalized_total_substrates[density_index] = 0;
        }
    }

    save_fluxes_to_csv(pCell, PhysiCell::PhysiCell_globals.current_time, "./output");
    return;
}

double dFBAIntracellular::get_flux_value(std::string reaction_name)
{
        
    auto it = this->sbml_model.getReaction(reaction_name);
    if (it != nullptr){
        return it->getFluxValue();
    }
    else{
        std::cout << "ERROR: Reaction " << reaction_name << " not found in the SBML model." << std::endl;
        exit(-1);
    }
}

void dFBAIntracellular::print_model(PhysiCell::Cell* pCell, double current_time, std::string output_folder){
    // Create output filename with cell ID and timestamp
    std::stringstream filename;
    int cell_id = pCell->ID;
    filename << output_folder << "/fba_model_cell_" << cell_id 
             << "_t_" << std::fixed << std::setprecision(2) << current_time << ".txt";
    
    std::ofstream outfile(filename.str());
    if (!outfile.is_open()) {
        std::cerr << "ERROR: Could not open file for writing: " << filename.str() << std::endl;
        return;
    }

    // Header with cell information
    outfile << "=================================================================" << std::endl;
    outfile << "FBA MODEL SUMMARY - Cell ID: " << cell_id << std::endl;
    outfile << "Time: " << current_time << " min" << std::endl;
    outfile << "SBML Model: " << this->sbml_filename << std::endl;
    outfile << "=================================================================" << std::endl;
    outfile << std::endl;

    // Growth parameters
    outfile << "--- GROWTH PARAMETERS ---" << std::endl;
    outfile << "Objective Reaction: " << this->objective_reaction << std::endl;
    outfile << "Current Growth Rate: " << this->current_growth_rate << " 1/h" << std::endl;
    outfile << "Cell Density: " << this->cell_density << " pg/um³" << std::endl;
    outfile << "Reference Volume: " << this->reference_volume << " um³" << std::endl;
    outfile << "dFBA Time Step: " << this->dfba_time_step << " min" << std::endl;
    outfile << std::endl;


    // Substrate exchanges (transport model)
    std::vector<double> density_vector = pCell->nearest_density_vector(); 
    
    outfile << "--- SUBSTRATE EXCHANGES ---" << std::endl;
    outfile << "Substrate Name\tCurrent concentration\tFBA Flux ID\tKm (mM)\tVmax (mmol/gDW/h)" << std::endl;
    for (const auto& exchange : this->substrate_exchanges) {
        const ExchangeFluxData& ex = exchange.second;
        double substrate_conc = density_vector[ex.density_index]; // mM = mmol/L
        outfile << ex.density_name << "\t" 
                << substrate_conc << "\t"
                << ex.fba_flux_id << "\t"
                << ex.Km.value << "\t"
                << ex.Vmax.value << std::endl;
    }
    outfile << std::endl;

    // Metabolites
    outfile << "--- METABOLITES (" << this->sbml_model.getNumMetabolites() << " total) ---" << std::endl;
    outfile << "ID\tName" << std::endl;
    for (const auto& m : this->sbml_model.getListOfMetabolites()) {
        outfile << m->getId() << "\t" << m->getName() << std::endl;
    }
    outfile << std::endl;

    // Reactions with detailed information
    outfile << "--- REACTIONS (" << this->sbml_model.getNumReactions() << " total) ---" << std::endl;
    outfile << "Reaction ID\tReaction Name\tLower Bound\tUpper Bound\tFlux Value\tReversible\tObjective Coeff\tReaction String" << std::endl;
    
    for (const auto& r : this->sbml_model.getListOfReactions()) {
        outfile << r->getId() << "\t"
                << r->getName() << "\t"
                << r->getLowerBound() << "\t"
                << r->getUpperBound() << "\t"
                << r->getFluxValue() << "\t"
                << (r->reversible() ? "Yes" : "No") << "\t"
                << r->getObjectiveCoefficient() << "\t"
                << r->getReactionString(this->sbml_model) << std::endl;
    }
    outfile << std::endl;

    // Boundary/Exchange reactions (for easier identification)
    std::vector<dFBAReaction*> boundary_rxns = this->sbml_model.getListOfBoundaryReactions();
    if (!boundary_rxns.empty()) {
        outfile << "--- BOUNDARY/EXCHANGE REACTIONS (" << boundary_rxns.size() << " total) ---" << std::endl;
        outfile << "Reaction ID\tLower Bound\tUpper Bound\tFlux Value" << std::endl;
        for (const auto& r : boundary_rxns) {
            outfile << r->getId() << "\t"
                    << r->getLowerBound() << "\t"
                    << r->getUpperBound() << "\t"
                    << r->getFluxValue() << std::endl;
        }
        outfile << std::endl;
    }

    outfile << "=================================================================" << std::endl;
    outfile.close();
    
    std::cout << "FBA model summary saved to: " << filename.str() << std::endl;
}

void dFBAIntracellular::save_dFBA(std::string path, std::string index) 
{
	
}

}













