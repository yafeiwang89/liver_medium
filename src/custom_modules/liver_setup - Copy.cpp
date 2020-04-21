/*
#############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the ver-  #
# sion number, such as below:                                               #
#                                                                           #
# We implemented and solved the model using PhysiCell (Version 0.0.5) [1].  #
#                                                                           #
# [1] A Ghaffarizadeh, SH Friedman, and P Macklin, PhysiCell: an open       #
#    source physics-based simulator for multicellular systemssimulator,     #
#    J. Comput. Biol., 2016 (submitted).                                    #
#                                                                           #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite       #
#     BioFVM as below:                                                      #
#                                                                           #
# We implemented and solved the model using PhysiCell (Version 0.0.5) [1],  #
# with BioFVM [2] to solve the transport equations.                         #
#                                                                           #
# [1] A Ghaffarizadeh, SH Friedman, and P Macklin, PhysiCell: an open       #
#    source physics-based multicellular simulator, J. Comput. Biol., 2016   #
#   (submitted).                                                            #
#                                                                           #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient     #
#    parallelized diffusive transport solver for 3-D biological simulations,#
#    Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730 #
#                                                                           #
#############################################################################
#                                                                           #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)   #
#                                                                           #
# Copyright (c) 2015-2016, Paul Macklin and the PhysiCell Project           #
# All rights reserved.                                                      #
#                                                                           #
# Redistribution and use in source and binary forms, with or without        #
# modification, are permitted provided that the following conditions are    #
# met:                                                                      #
#                                                                           #
# 1. Redistributions of source code must retain the above copyright notice, #
# this list of conditions and the following disclaimer.                     #
#                                                                           #
# 2. Redistributions in binary form must reproduce the above copyright      #
# notice, this list of conditions and the following disclaimer in the       #
# documentation and/or other materials provided with the distribution.      #
#                                                                           #
# 3. Neither the name of the copyright holder nor the names of its          #
# contributors may be used to endorse or promote products derived from this #
# software without specific prior written permission.                       #
#                                                                           #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       #
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED #
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           #
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER #
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  #
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       #
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        #
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    #
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      #
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        #
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              #
#                                                                           #
#############################################################################
*/

#include "./liver_setup.h"

std::vector< std::vector<double> > central_vein_positions;
Cell_Container* cell_container;

Cell_Definition parenchyme;

Cell_Definition portal_triad;

Cell_Definition HCT116;

Cell_Definition MCF10A;

double dt;
Custom_Cell_Data liver_custom_data;
Custom_Cell_Data empty_data;

void set_program_metadata( void )
{
	// set up metadata

	BioFVM_metadata.program.user.surname = "Macklin";
	BioFVM_metadata.program.user.given_names = "Paul";
	BioFVM_metadata.program.user.email = "Paul.Macklin@MathCancer.org";
	BioFVM_metadata.program.user.organization = "Indiana University";
	BioFVM_metadata.program.user.department = "Intelligent Systems Engineering";

	BioFVM_metadata.program.creator.given_names = "Paul";
	BioFVM_metadata.program.creator.email = "Paul.Macklin@MathCancer.org";
	BioFVM_metadata.program.creator.organization = "Indiana University";
	BioFVM_metadata.program.creator.department = "Intelligent Systems Engineering";

	BioFVM_metadata.program.program_name = "2D_liver";
	BioFVM_metadata.program.program_version = "0.0.0";
	BioFVM_metadata.program.program_URL = "http://MathCancer.org";
}

void initialize_liver_cell_definitions( void )
{
	int number_of_substrates = 1;

	double _pi_ = 3.141592653589793;
	double four_thirds_pi = 4.0*_pi_/3.0;

	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );

	std::cout << __FUNCTION__ << " " << cell_defaults.phenotype.secretion.secretion_rates << " xx " <<  std::endl;

	// switch to the live/dead model
	cell_defaults.phenotype.cycle.sync_to_cycle_model( live );

	// O2-dependent birth

	cell_defaults.functions.update_phenotype = update_cell_and_death_parameters_O2_based;

	// 2D

	cell_defaults.functions.set_orientation = up_orientation;
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true;

	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );

	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );

	int spring_constant_index = liver_custom_data.find_variable_index( "spring constant" );
	int mechanical_relaxation_rate_index = liver_custom_data.find_variable_index( "mechanical relaxation rate" );	

	int oxygen_index = microenvironment.find_density_index( "oxygen" );

	// type 0 : tumor
	// type 1 : central vein (remove!)
	// type 2 : portal triad
	// type 3 : parenchyme



	// begin defining parenchyme

	parenchyme = cell_defaults;
	parenchyme.name = "liver parenchyme";
	parenchyme.type = 3;

	// update its phenotype

	// size information

	double radius = 15.0;
	double nuclear_radius = radius;
	parenchyme.phenotype.volume.total = four_thirds_pi*pow( radius , 3 );
	parenchyme.phenotype.volume.fluid_fraction = 0.75;
	parenchyme.phenotype.volume.fluid = parenchyme.phenotype.volume.fluid_fraction * parenchyme.phenotype.volume.total;
	parenchyme.phenotype.volume.solid = parenchyme.phenotype.volume.total - parenchyme.phenotype.volume.fluid;

	parenchyme.phenotype.volume.nuclear = four_thirds_pi*pow( nuclear_radius , 3 ) ;
	parenchyme.phenotype.volume.nuclear_fluid = parenchyme.phenotype.volume.fluid_fraction * parenchyme.phenotype.volume.nuclear;
	parenchyme.phenotype.volume.nuclear_solid = parenchyme.phenotype.volume.nuclear - parenchyme.phenotype.volume.nuclear_fluid;

	parenchyme.phenotype.volume.cytoplasmic = parenchyme.phenotype.volume.total - parenchyme.phenotype.volume.nuclear;
	parenchyme.phenotype.volume.cytoplasmic_fluid = parenchyme.phenotype.volume.fluid_fraction * parenchyme.phenotype.volume.cytoplasmic;
	parenchyme.phenotype.volume.cytoplasmic_solid = parenchyme.phenotype.volume.cytoplasmic - parenchyme.phenotype.volume.cytoplasmic_fluid;
	parenchyme.phenotype.geometry.update( NULL, parenchyme.phenotype , 0 ); // self-consistency

	parenchyme.phenotype.volume.target_fluid_fraction = 0.75;
	parenchyme.phenotype.volume.target_solid_nuclear = parenchyme.phenotype.volume.nuclear_solid;
	parenchyme.phenotype.volume.target_cytoplasmic_to_nuclear_ratio = 0;
	parenchyme.phenotype.volume.cytoplasmic_to_nuclear_ratio = 0;
	parenchyme.phenotype.volume.target_solid_cytoplasmic = parenchyme.phenotype.volume.cytoplasmic_to_nuclear_ratio * parenchyme.phenotype.volume.target_solid_nuclear;

	std::cout << parenchyme.phenotype.volume.target_solid_nuclear << " " <<  parenchyme.phenotype.volume.nuclear_solid << std::endl
		<< parenchyme.phenotype.volume.target_solid_cytoplasmic << " " <<  parenchyme.phenotype.volume.cytoplasmic_solid << std::endl;

	// quirk of the live model: cells double their target volumes upon cycle entry.
	// Fix this by getting rid of the entry function. Make them inert.
	// parenchyme.phenotype.cycle.model().phases[0].entry_function = NULL;

	// no birth or apoptosis

	parenchyme.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0.0;
	parenchyme.phenotype.death.rates[apoptosis_index] = 0.0;

	// mechanics: reduce mechanics, because these tissues are partly voide

	
	parenchyme.phenotype.mechanics.cell_cell_adhesion_strength *= 0.75;                                  
	parenchyme.phenotype.mechanics.cell_cell_repulsion_strength *= 0.75;
	

	// now set upt the custom data (as needed)

	parenchyme.custom_data = liver_custom_data;

	std::cout << cell_defaults.custom_data << std::endl;
	std::cout << parenchyme.custom_data << std::endl;

	// now set up functions

	parenchyme.functions.update_phenotype = parenchyme_phenotype_update_function;
	parenchyme.functions.custom_cell_rule = liver_plastoelastic_mechanics;

	/*
	Cell_Parameters parameters;
	Custom_Cell_Data custom_data;
	Cell_Functions functions;
	Phenotype phenotype;
*/

// In the liver paper, we don't consider the portal traid and remove them, filled with parenchyme !!!!!!!!!!!!!!!!!!
/*
	// begin defining the portal triad
	
    // directly use all information of parenchyme, don't need to modify anything of portal triad 
	portal_triad = parenchyme;       
	portal_triad.name = "portal triad";
	portal_triad.type = 2;

	// update its phenotype

	// size information

	radius = 20.0;
	nuclear_radius = radius;
	portal_triad.phenotype.volume.total = four_thirds_pi*pow( radius , 3 );
	portal_triad.phenotype.volume.fluid_fraction = 0.75;
	portal_triad.phenotype.volume.fluid = portal_triad.phenotype.volume.fluid_fraction * portal_triad.phenotype.volume.total;
	portal_triad.phenotype.volume.solid = portal_triad.phenotype.volume.total - portal_triad.phenotype.volume.fluid;

	portal_triad.phenotype.volume.nuclear = four_thirds_pi*pow( nuclear_radius , 3 ) ;
	portal_triad.phenotype.volume.nuclear_fluid = portal_triad.phenotype.volume.fluid_fraction * portal_triad.phenotype.volume.nuclear;
	portal_triad.phenotype.volume.nuclear_solid = portal_triad.phenotype.volume.nuclear - portal_triad.phenotype.volume.nuclear_fluid;

	portal_triad.phenotype.volume.cytoplasmic = portal_triad.phenotype.volume.total - portal_triad.phenotype.volume.nuclear;
	portal_triad.phenotype.volume.cytoplasmic_fluid = portal_triad.phenotype.volume.fluid_fraction * portal_triad.phenotype.volume.cytoplasmic;
	portal_triad.phenotype.volume.cytoplasmic_solid = portal_triad.phenotype.volume.cytoplasmic - portal_triad.phenotype.volume.cytoplasmic_fluid;
	portal_triad.phenotype.geometry.update( NULL, portal_triad.phenotype , 0 ); // self-consistency

	portal_triad.phenotype.volume.target_fluid_fraction = 0.75;
	portal_triad.phenotype.volume.target_solid_nuclear = portal_triad.phenotype.volume.nuclear_solid;
	portal_triad.phenotype.volume.cytoplasmic_to_nuclear_ratio = 0;
	portal_triad.phenotype.volume.target_solid_cytoplasmic = portal_triad.phenotype.volume.cytoplasmic_to_nuclear_ratio * portal_triad.phenotype.volume.target_solid_nuclear;

	// no birth or apoptosis

	portal_triad.phenotype.cycle.model().transition_rate( cycle_start_index , cycle_end_index ) = 0.0;
	portal_triad.phenotype.death.rates[apoptosis_index] = 0.0;

	// mechanics: turn down adhesion, so that these guys just are rendered but barely interact

	// portal_triad.phenotype.mechanics.cell_cell_adhesion_strength *= 0.1;
	// portal_triad.phenotype.mechanics.cell_cell_repulsion_strength *= 0.1;

	// now set upt the custom data (as needed)

	portal_triad.custom_data = liver_custom_data;

	
	// make portal_triad don't move 
	// portal_triad.phenotype.motility.is_motile = false; 
	// portal_triad.custom_data[spring_constant_index] *= 10.0;                                
	// portal_triad.custom_data[mechanical_relaxation_rate_index] *= 0.0;

	
	std::cout << portal_triad.custom_data << std::endl;

	// now set up functions
    
	// the phenotype function of portal triad is same as parenchyma
	portal_triad.functions.update_phenotype = parenchyme_phenotype_update_function;

*/

	HCT116 = cell_defaults;
	HCT116.name = "HCT 116";
	HCT116.type = 0;
	HCT116.phenotype.secretion.uptake_rates[oxygen_index] = 10;

	// set up the default functions
	HCT116.functions.cycle_model = Ki67_advanced;
	// HCT116.phenotype.cycle.data.transition_rate(0,0) = 1.0/(120.0);
	HCT116.phenotype.cycle.sync_to_cycle_model( HCT116.functions.cycle_model ); //

	HCT116.functions.custom_cell_rule = NULL;
	HCT116.functions.update_phenotype = HCT116_phenotype_update_function;
	// set up phenotype

	HCT116.parameters.pReference_live_phenotype = &( HCT116.phenotype );
	
	// set the tumor's birth and death rate as zero for the static studying
	HCT116.phenotype.cycle.data.transition_rate(0,1) = parameters.doubles( "tumor_transition_rate" ); //1.0 / ( 7.26 * 60.0); // so that the total cycle time is ~ 1 / .043
	HCT116.parameters.max_necrosis_rate = parameters.doubles( "tumor_max_necrosis_rate" );  //1.0 / ( 6.0 * 60.0 ); // 6 hour survival time
/*
	// figure out which death model is apoptosis.
    // int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );          

    // set apoptosis to trigger quickly
    HCT116.phenotype.death.rates[apoptosis_index] = 1.0 / ( 3.0 * 60.0 );   

    // prevent necrosis
    // HCT116.parameters.max_necrosis_rate = 0.0;

    // make apoptosis proceed very fast, once cells are apoptotic
    HCT116.phenotype.death.models[apoptosis_index]->transition_rate(0,1) = 1.0 / ( 5.0 * 60.0 );
*/  


	// parenchyme
/*


	int i = 0;
	// set normoxic (21% O2) phenotype
	std::unordered_map<int,int> phase_hashmap;//=new std::unordered_map<int,int>;

	// sizing information
	parenchyme.phenotype_names[i] = "viable";
	parenchyme.microenvironment_samples[i].variables.resize( number_of_substrates ) ;
	parenchyme.microenvironment_samples[i].variables[0] = "oxygen";
	parenchyme.microenvironment_samples[i].densities.resize( number_of_substrates ) ;
	parenchyme.microenvironment_samples[i].densities[0] = 38;

	parenchyme.phenotypes[i].secretion_rates.rates.resize( number_of_substrates, 0.0 );
	parenchyme.phenotypes[i].uptake_rates.rates.resize( number_of_substrates , 0.0 ) ;
	parenchyme.phenotypes[i].saturation_densities.densities.resize( number_of_substrates , 0.0 );


	// note: code does not allow cell that is zero nucleus
	// future: put a check into volume code, if target_CN > 9e9, assume all cytoplasmic ?

	// cell cycle information
	parenchyme.phenotypes[i].cycle.cycle_model = PhysiCell_constants::live_apoptotic_cycle_model;
	parenchyme.phenotypes[i].cycle.cycle_model_name = "live_apoptotic";
	parenchyme.phenotypes[i].cycle.phases.resize( 2 );
	parenchyme.phenotypes[i].cycle.phases[0].code = PhysiCell_constants::live;
	parenchyme.phenotypes[i].cycle.phases[0].name = "live";
	parenchyme.phenotypes[i].cycle.phases[0].elapsed_time = 0.0;
	parenchyme.phenotypes[i].cycle.phases[0].duration = 9e9; //
	parenchyme.phenotypes[i].cycle.phases[0].birth_rate = 0;
	parenchyme.phenotypes[i].cycle.phases[0].death_rate = 0  / 60; // arbitrarily set
	// parenchyme.phenotypes[i].cycle.phases[0].death_type = PhysiCell_constants::apoptotic;
	parenchyme.phenotypes[i].cycle.phases[0].arrest_density = pow( 0.006012 , 1.5 ); // cells per cubic micron
	parenchyme.phenotypes[i].cycle.phases[0].volume_change_timescale_N = 9e9; // 9.1*60.0;
	parenchyme.phenotypes[i].cycle.phases[0].volume_change_timescale_C = 9e9; // 11.1*60.0;
	parenchyme.phenotypes[i].cycle.phases[0].volume_change_timescale_F = 9e9; // 1*60.0;
	phase_hashmap[(int)parenchyme.phenotypes[i].cycle.phases[0].code]=0;

	parenchyme.phenotypes[i].cycle.phases[1].code = PhysiCell_constants::apoptotic;
	parenchyme.phenotypes[i].cycle.phases[1].name = "apoptotic";
	parenchyme.phenotypes[i].cycle.phases[1].elapsed_time = 0.0;
	parenchyme.phenotypes[i].cycle.phases[1].duration = 8.6 * 60.0; // minutes (MDA-MB-231)
	parenchyme.phenotypes[i].cycle.phases[1].birth_rate = 0;
	parenchyme.phenotypes[i].cycle.phases[1].death_rate = 0; // minutes (MDA-MB-231)
	// parenchyme.phenotypes[i].cycle.phases[1].death_type = PhysiCell_constants::apoptotic;
	parenchyme.phenotypes[i].cycle.phases[1].arrest_density = 9e99; // cells per cubic micron
	parenchyme.phenotypes[i].cycle.phases[1].volume_change_timescale_N = 8.6*60.0;
	parenchyme.phenotypes[i].cycle.phases[1].volume_change_timescale_C = 3*60.0;
	parenchyme.phenotypes[i].cycle.phases[1].volume_change_timescale_F = 1*60.0;
	phase_hashmap[(int)parenchyme.phenotypes[i].cycle.phases[1].code]=1;

	parenchyme.display_information( std::cout );
	std::cout << std::endl ;

	parenchyme.phenotypes[i].set_phase_maps(new std::unordered_map<int,int>(phase_hashmap));

	// parenchyme.phenotypes[i].set_phase_maps(&phase_hashmap);

	// HCT116;
	// portal_triad;
*/

	return;
}


void setup_liver( void )
{
	set_program_metadata();
	setup_custom_data();
	setup_liver_microenvironment( microenvironment, "../data/X.mat" , "../data/Y.mat", "../data/oxygen.mat");
	initialize_liver_cell_definitions();
	read_liver_cells( "../data/cells.mat" );

//	read_liver( microenvironment, "./config/X.mat" , "./config/Y.mat", "./config/oxygen.mat", "./config/cells.mat" );

	return;
}

void setup_liver( std::string& config_file )
{
	// the file should tell us where to find data.


	return;
}

// set up custom data
void setup_custom_data( void )
{
/*
	std::vector<double> ECM_attach;

	double spring_constant;
	double mechanical_relaxation_rate;

	double current_mechanical_strain;

	double max_mechanical_strain;
	double max_integrated_mechanical_strain;

*/
	empty_data = cell_defaults.custom_data;

	// mechanics

//	Vector_Variable myvec;
//	std::cout << myvec << std::endl;
	std::vector<double> myvec;
	myvec.resize(3,0.0);

	// assume elastic movement on the order of 10 min at maximum 10 micron elongation
	
	liver_custom_data.add_variable( "spring constant", "1/min" , parameters.doubles( "elastic_rate" ) );  // 0.05    // (1.0/10.0) * (1.0/10.0)
	// assume plastic movement on the order of 1 day at maximum 10 micron elongation
	liver_custom_data.add_variable( "mechanical relaxation rate", "1/min" , parameters.doubles( "plastic_rate" ) );  // 0.0005  // (1.0/10.0) * (1.0/(24.0*60.0)

	liver_custom_data.add_variable( "mechanical strain", "micron" , 0.0 );
	liver_custom_data.add_variable( "integrated mechanical strain", "micron*min" , 0.0 );
	
	liver_custom_data.add_variable( "pressure", "pascal" , 0.5 );  // add the pressure variable !!!!!!
	liver_custom_data.add_variable( "birth_rate", "1/min" , 0.05 );  // add the birth_rate variable !!!!!

	liver_custom_data.add_variable( "max mechanical strain", "micron" , parameters.doubles( "max_ECM_displacement" )); // 0.75       // 10.0
	liver_custom_data.add_variable( "max integrated mechanical strain", "micron*min" , parameters.doubles( "max_ECM_displacement" ) * 60.0 ); //0.75   // 1 hour of deformation by 10 microns  

	liver_custom_data.add_vector_variable( "ECM attachment point", "micron", myvec );
	
	// cell_defaults.custom_data = liver_custom_data; //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	liver_custom_data.add_vector_variable( "mechanical strain displacement" , "micron", myvec );
	
	cell_defaults.custom_data = liver_custom_data; //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	int i = liver_custom_data.find_variable_index( "elastic parameter" );
	std::cout << liver_custom_data[i] << std::endl;

	return;
}

void setup_liver_microenvironment( Microenvironment& M , std::string Xfile , std::string Yfile, std::string OxygenFile )
{
	// read X
	std::vector< std::vector<double> > X = read_matlab( Xfile );
	// determine dx
	double dx = X[0][1] - X[0][0];
	int Xnodes = X[0].size();

	// read Y
	std::vector< std::vector<double> > Y = read_matlab( Yfile );
	// determine dy
	double dy = Y[0][1] - Y[0][0];
	int Ynodes = Y[0].size();

	// determine dz
	double dz = 0.5* ( dx + dy );

	// set the microenvironment options
	default_microenvironment_options.name = "liver";
	default_microenvironment_options.simulate_2D = true;
	default_microenvironment_options.outer_Dirichlet_conditions = false;
	default_microenvironment_options.calculate_gradients = false;

	default_microenvironment_options.dx = dx;
	default_microenvironment_options.dy = dy;
	default_microenvironment_options.dz = dz;

	default_microenvironment_options.X_range[0] = X[0][0] - dx/2.0;
	default_microenvironment_options.X_range[1] = X[0][Xnodes-1] + dx/2.0;

	default_microenvironment_options.Y_range[0] = Y[0][0] - dy/2.0;
	default_microenvironment_options.Y_range[1] = Y[0][Ynodes-1] + dy/2.0;

	default_microenvironment_options.Z_range[0] = -dz/2.0;
	default_microenvironment_options.Z_range[1] = dz/2.0;


	// apply the microenvironment options
	initialize_microenvironment();

	// set up the microenvironment variable (oxygen)
	microenvironment.set_density( 0 , "oxygen" , "mmHg" );
	microenvironment.diffusion_coefficients[0] = 1e5;
	microenvironment.decay_rates[0] = 0.1; //  1 mm length scale

	// read oxygen
	std::vector< std::vector<double> > oxygen = read_matlab( OxygenFile );

	#pragma omp parallel for
	for( int n=0 ; n < microenvironment.number_of_voxels() ; n++ )
	{
		std::vector<int> indices = microenvironment.cartesian_indices( n );
		int i = indices[0];
		int j = indices[1];
		int k = indices[2];

		// use this syntax to access the coordinates (as a vector) of
		// the ith voxel;
		// microenvironment.mesh.voxels[i].center

		// use this access the jth substrate at the ith voxel
		microenvironment.density_vector(n)[0] = oxygen[i][j];
		// microenvironment.density_vector(n)[1] = oxygen[i][j];
	}

	// create Dirichlet nodes

	for( int i=0; i < microenvironment.number_of_voxels() ; i++ )
	{
		microenvironment.add_dirichlet_node( i, microenvironment.density_vector(i) );
	}

	// set up the diffusion solver, sources and sinks

	microenvironment.diffusion_decay_solver = diffusion_decay_solver__constant_coefficients_LOD_2D;
/*
	microenvironment.bulk_supply_rate_function = liver_supply_function;
	microenvironment.bulk_supply_target_densities_function = liver_supply_target_function;
	microenvironment.bulk_uptake_rate_function = liver_uptake_function;
*/

	return;
}

// read liver from a series of matlab files
void read_liver( Microenvironment& M , std::string Xfile , std::string Yfile, std::string OxygenFile, std::string CellsFile )
{
/*
	// read X
	std::vector< std::vector<double> > X = read_matlab( Xfile );
	// determine dx
	double dx = X[0][1] - X[0][0];
	int Xnodes = X[0].size();

	// read Y
	std::vector< std::vector<double> > Y = read_matlab( Yfile );
	// determine dy
	double dy = Y[0][1] - Y[0][0];
	int Ynodes = Y[0].size();

	// determine dz
	double dz = 0.5* ( dx + dy );

	// set the microenvironment options
	default_microenvironment_options.name = "liver";
	default_microenvironment_options.simulate_2D = true;
	default_microenvironment_options.outer_Dirichlet_conditions = false;
	default_microenvironment_options.calculate_gradients = true;

	default_microenvironment_options.dx = dx;
	default_microenvironment_options.dy = dy;
	default_microenvironment_options.dz = dz;

	default_microenvironment_options.X_range[0] = X[0][0] - dx/2.0;
	default_microenvironment_options.X_range[1] = X[0][Xnodes-1] + dx/2.0;

	default_microenvironment_options.Y_range[0] = Y[0][0] - dy/2.0;
	default_microenvironment_options.Y_range[1] = Y[0][Ynodes-1] + dy/2.0;

	default_microenvironment_options.Z_range[0] = -dz/2.0;
	default_microenvironment_options.Z_range[1] = dz/2.0;


	// apply the microenvironment options
	initialize_microenvironment();

	// set up the microenvironment variable (oxygen)
	microenvironment.set_density( 0 , "oxygen" , "mmHg" );
	microenvironment.diffusion_coefficients[0] = 1e5;
	microenvironment.decay_rates[0] = 0.1; //  1 mm length scale

	// read oxygen
	std::vector< std::vector<double> > oxygen = read_matlab( OxygenFile );

	#pragma omp parallel for
	for( int n=0 ; n < microenvironment.number_of_voxels() ; n++ )
	{
		std::vector<int> indices = microenvironment.cartesian_indices( n );
		int i = indices[0];
		int j = indices[1];
		int k = indices[2];

		// use this syntax to access the coordinates (as a vector) of
		// the ith voxel;
		// microenvironment.mesh.voxels[i].center

		// use this access the jth substrate at the ith voxel
		microenvironment.density_vector(n)[0] = oxygen[i][j];
		// microenvironment.density_vector(n)[1] = oxygen[i][j];
	}

	// create Dirichlet nodes

	for( int i=0; i < microenvironment.number_of_voxels() ; i++ )
	{
		microenvironment.add_dirichlet_node( i, microenvironment.density_vector(i) );
	}

	// set up the diffusion solver, sources and sinks

	microenvironment.diffusion_decay_solver = diffusion_decay_solver__constant_coefficients_LOD_2D;

	microenvironment.bulk_supply_rate_function = liver_supply_function;
	microenvironment.bulk_supply_target_densities_function = liver_supply_target_function;
	microenvironment.bulk_uptake_rate_function = liver_uptake_function;
*/
	// PhysiCell setup
	cell_container = new Cell_Container;
	int mechanics_voxel_size = 30.0;
	cell_container->initialize( microenvironment.mesh.bounding_box[0] , microenvironment.mesh.bounding_box[3] ,
		microenvironment.mesh.bounding_box[1] , microenvironment.mesh.bounding_box[4] ,
		microenvironment.mesh.bounding_box[2] , microenvironment.mesh.bounding_box[5] , mechanics_voxel_size );
	microenvironment.agent_container = (Agent_Container*) cell_container;

	// read cells

	std::cout << __LINE__ << std::endl;
	std::vector<double> template_vector3(3,0.0);
	central_vein_positions.resize( 1 ,  template_vector3 );

	std::cout << __LINE__ << std::endl;
	std::vector< std::vector<double> > cells = read_matlab( CellsFile );
	// columns: x , y ,z , radius , type
	// (0,1,2) : (x,y,z)
	// 3 : radius
	// 4: type
	// type 0 : tumor
	// type 1 : central vein (remove!)
	// type 2 : portal triad
	// type 3 : parenchyme

	// create basic agents (later swap out PhysiCell

	std::cout << cells.size() << " x " << cells[0].size() << std::endl;

	std::cout << __LINE__ << std::endl;
	for( int n=0; n < cells.size() ; n++ )
	{
		if( fabs( cells[n][4] - 1.0 ) < 1e-10 ) // type 1 : central vein
		{
			template_vector3[0] = cells[n][0];
			template_vector3[1] = cells[n][1];
			central_vein_positions.push_back( template_vector3 );
		}

		if( fabs( cells[n][4] - 2.0 ) < 1e-10 ) // type 2 : portal triad
		{
			template_vector3[0] = cells[n][0];
			template_vector3[1] = cells[n][1];
		}

		if( fabs( cells[n][4] - 3.0 ) < 1e-10 ) // type 3 : parenchyme
		{
			// template_vector3[0] = cells[n][0];
			// template_vector3[1] = cells[n][1];

			Cell* pCell = create_cell( parenchyme );
			pCell->assign_position( cells[n][0] , cells[n][1] , 0.0 );

			/*
			pCell->set_phenotype( parenchyme.phenotypes[0] );
			pCell->advance_cell_current_phase = do_nothing;

			pCell->phenotype.set_current_phase( PhysiCell_constants::live );
			pCell->set_internal_uptake_constants(dt);

			pCell->custom_data = default_custom_data;

			pCell->custom_data.ECM_attach = pCell->position;
			pCell->set_internal_uptake_constants(dt);

			pCell->custom_cell_rule = parenchyme_cell_rule;
			pCell->custom_data.cell_type = (int) cells[n][4] ;
			*/

		}

		if( fabs( cells[n][4] - 0 ) < 1e-10 ) // type  :  tumor
		{
			template_vector3[0] = cells[n][0];
			template_vector3[1] = cells[n][1];

			Cell* pCell = create_cell( HCT116 );
			pCell->assign_position( cells[n][0] , cells[n][1] , 0.0 );

/*
			Basic_Agent* pBA = create_basic_agent();
			pBA->assign_position( cells[n][0] , cells[n][1] , 0.0 );
			pBA->set_total_volume( 4.188790204786391 * cells[n][3] * cells[n][3] * cells[n][3] );
			pBA->set_internal_uptake_constants( 0.01 );

			pBA->type = (int) cells[n][4] ;
*/
		}
/*
	double get_total_volume();
	void set_total_volume(double);
	void update_voxel_index();

	void set_internal_uptake_constants( double dt ); // any time you update the cell volume or rates, should call this function.

	void register_microenvironment( Microenvironment* );
	Microenvironment* get_microenvironment( void );

	int ID;
	int index;
	int type;
*/

	}
	std::cout << __LINE__ << std::endl;
	std::cout << "basic agents: " << all_basic_agents.size() << std::endl;


//	microenvironment.

//resize( double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , double dx, double dy, double dz );
	/*
					// determine the number of voxels
				int rows;
				int columns;
				FILE* fp = read_matlab_header( &rows, &columns, node.text().get() );
				int voxel_count = columns;

				// resize the appropriate data structure
				M_destination.resize_voxels( voxel_count );


				// read the data directly into the voxels
				for( int j=0; j < columns ; j++ )
				{
					double temp;
					// read x, y, z, dV
					fread( (char*) & (M_destination.mesh.voxels[j].center[0])   , sizeof(double) , 1 , fp );
					fread( (char*) & (M_destination.mesh.voxels[j].center[1])   , sizeof(double) , 1 , fp );
					fread( (char*) & (M_destination.mesh.voxels[j].center[2])   , sizeof(double) , 1 , fp );
					fread( (char*) & (M_destination.mesh.voxels[j].volume)   , sizeof(double) , 1 , fp );
				}
				fclose( fp );
				*/






	return;
}

// read liver from a series of matlab files
void read_liver_cells( std::string CellsFile )
{
	// PhysiCell setup
	cell_container = new Cell_Container;
	cell_container->initialize( microenvironment.mesh.bounding_box[0] , microenvironment.mesh.bounding_box[3] ,
		microenvironment.mesh.bounding_box[1] , microenvironment.mesh.bounding_box[4] ,
		microenvironment.mesh.bounding_box[2] , microenvironment.mesh.bounding_box[5] , 30.0 );
	microenvironment.agent_container = (Agent_Container*) cell_container;

	// read cells

	std::cout << __LINE__ << std::endl;
	std::vector<double> template_vector3(3,0.0);
	central_vein_positions.resize( 1 ,  template_vector3 );

	std::cout << __LINE__ << std::endl;
	std::vector< std::vector<double> > cells = read_matlab( CellsFile );
	// columns: x , y ,z , radius , type
	// (0,1,2) : (x,y,z)
	// 3 : radius
	// 4: type
	// type 0 : tumor
	// type 1 : central vein (remove!)
	// type 2 : portal triad
	// type 3 : parenchyme

	// create basic agents (later swap out PhysiCell

	std::cout << cells.size() << " x " << cells[0].size() << std::endl;

	std::cout << __LINE__ << std::endl;
	for( int n=0; n < cells.size() ; n++ )
	{
		if( fabs( cells[n][4] - 1.0 ) < 1e-10 ) // type 1 : central vein
		{

			std::cout << "central vein radius: " << cells[n][3] << std::endl;

			template_vector3[0] = cells[n][0];
			template_vector3[1] = cells[n][1];
			central_vein_positions.push_back( template_vector3 );
		}


		if( fabs( cells[n][4] - 2.0 ) < 1e-10 ) // type 2 : portal triad
		{
			std::cout << "portal triad radius: " << cells[n][3] << std::endl;

			Cell* pCell = create_cell( portal_triad );
			pCell->assign_position( cells[n][0] , cells[n][1] , 0.0 );
			pCell->custom_data.vector_variables[0].value = pCell->position;
		}


		if( fabs( cells[n][4] - 3.0 ) < 1e-10 ) // type 3 : parenchyme
		{
			// template_vector3[0] = cells[n][0];
			// template_vector3[1] = cells[n][1];

			// create a parenchyma cell, and set its ECM attachment point to
			// its current location

			Cell* pCell = create_cell( parenchyme );
			pCell->assign_position( cells[n][0] , cells[n][1] , 0.0 );
			pCell->custom_data.vector_variables[0].value = pCell->position;

			/*
			pCell->set_phenotype( parenchyme.phenotypes[0] );
			pCell->advance_cell_current_phase = do_nothing;

			pCell->phenotype.set_current_phase( PhysiCell_constants::live );
			pCell->set_internal_uptake_constants(dt);

			pCell->custom_data = default_custom_data;

			pCell->custom_data.ECM_attach = pCell->position;
			pCell->set_internal_uptake_constants(dt);

			pCell->custom_cell_rule = parenchyme_cell_rule;
			pCell->custom_data.cell_type = (int) cells[n][4] ;
			*/

		}

		if( fabs( cells[n][4] - 0 ) < 1e-10 ) // type  :  tumor
		{
			template_vector3[0] = cells[n][0];
			template_vector3[1] = cells[n][1];

			Cell* pCell = create_cell( HCT116 );
			pCell->assign_position( cells[n][0] , cells[n][1] , 0.0 );

/*
			Basic_Agent* pBA = create_basic_agent();
			pBA->assign_position( cells[n][0] , cells[n][1] , 0.0 );
			pBA->set_total_volume( 4.188790204786391 * cells[n][3] * cells[n][3] * cells[n][3] );
			pBA->set_internal_uptake_constants( 0.01 );

			pBA->type = (int) cells[n][4] ;
*/
		}

	}
	std::cout << __LINE__ << std::endl;
	std::cout << "basic agents: " << all_basic_agents.size() << std::endl;

	return;
}



void liver_supply_function( Microenvironment* microenvironment, int voxel_index, std::vector<double>* write_here )
{
	// use this syntax to access the jth substrate write_here
	// (*write_here)[j]
	// use this syntax to access the jth substrate in voxel voxel_index of microenvironment:
	// microenvironment->density_vector(voxel_index)[j]

	return;
}

void liver_supply_target_function( Microenvironment* microenvironment, int voxel_index, std::vector<double>* write_here )
{
	// use this syntax to access the jth substrate write_here
	// (*write_here)[j]
	// use this syntax to access the jth substrate in voxel voxel_index of microenvironment:
	// microenvironment->density_vector(voxel_index)[j]

	return;
}

void liver_uptake_function( Microenvironment* microenvironment, int voxel_index, std::vector<double>* write_here )
{
	// use this syntax to access the jth substrate write_here
	// (*write_here)[j]
	// use this syntax to access the jth substrate in voxel voxel_index of microenvironment:
	// microenvironment->density_vector(voxel_index)[j]

	return;
}

void setup( void )
{
	// read input file

	// set metadata

	// create digital cell lines (later, load these)

	// get liver data file names

	// read in liver








}
