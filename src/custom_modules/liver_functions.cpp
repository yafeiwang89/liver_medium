/*
#############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the ver-  #
# sion number, such as below:                                               #
#                                                                           #
# We implemented and solved the model using PhysiCell (Version 0.0.5) [1].  #
#                                                                           #
# [1] A Ghaffarizadeh, SH Friedman, and P Macklin, PhysiCell: an open       #
#    source physics-based simulator for multicellular systemssimulator, 	#
#	 J. Comput. Biol., 2016 (submitted).                                    #
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

#include "./liver_functions.h"
#include "../modules/PhysiCell_pathology.h"

void parenchyme_phenotype_update_function( Cell* pCell , Phenotype& phenotype, double dt )
{
	static double t = 0.0;

	extern Custom_Cell_Data liver_custom_data;

	static int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );

	// // oxygen-based necrosis
	// update_cell_and_death_parameters_O2_based( pCell, phenotype, dt );

	// update the cell strain and integrated strain

	static int strain_index = liver_custom_data.find_variable_index( "mechanical strain" );
	static int integrated_strain_index = liver_custom_data.find_variable_index( "integrated mechanical strain" );

	// displacement
	pCell->custom_data.vector_variables[1].value = pCell->custom_data.vector_variables[0].value;
	pCell->custom_data.vector_variables[1].value -= pCell->position;

	// pCell->custom_data[strain_index] = norm( pCell->position - pCell->custom_data.vector_variables[0].value );
	pCell->custom_data[strain_index] = norm( pCell->custom_data.vector_variables[1].value );
	pCell->custom_data[integrated_strain_index] += ( dt * pCell->custom_data[strain_index] );

	// check for the custom apoptosis
	// Set to very high apoptosis rates
	// And get rid of the Dirichlet node


	if( strain_based_apoptosis( pCell ) )
	{
		pCell->phenotype.death.rates[apoptosis_index] = 9e9;
	}
/*
	if( strain_based_apoptosis( pCell ) )
	{
		pCell->phenotype.death.rates[apoptosis_index] = 9e9;

		int voxel_index = pCell->get_current_voxel_index();
//		microenvironment.remove_dirichlet_node( pCell->get_current_voxel_index() );

		if( pCell->type == 3 ) // parenchyma
		{
			microenvironment.remove_dirichlet_node( voxel_index );
		}
		else
		{
			std::vector<int> indices = microenvironment.mesh.cartesian_indices( voxel_index );
			for( int i= indices[0]-1; i <= indices[0]+1; i++ )
			{
				for( int j=indices[1]-1; j <= indices[1]+1 ; j++ )
				{
					int n = microenvironment.mesh.voxel_index(i,j,0);
					microenvironment.remove_dirichlet_node( n );
				}
			}
		}

		// microenvironment.remove_dirichlet_node( voxel_index );

//		microenvironment.dirichlet_value_vectors[ pCell->get_current_voxel_index() ][0] = 100.0;
	}
*/

	t += dt;

	return;
}


bool strain_based_apoptosis( Cell* pCell )
{
	if( pCell->phenotype.death.dead )
	{ return false; }

	extern Custom_Cell_Data liver_custom_data;

	static int strain_index = liver_custom_data.find_variable_index( "mechanical strain" );
	static int max_strain_index = liver_custom_data.find_variable_index( "max mechanical strain" );

	if( pCell->custom_data[strain_index] > pCell->custom_data[max_strain_index] )
	{
		return true;
	}
	return false;
}


bool integrated_strain_based_apoptosis( Cell* pCell )
{
	if( pCell->phenotype.death.dead )
	{ return false; }

	extern Custom_Cell_Data liver_custom_data;

	static int integrated_strain_index = liver_custom_data.find_variable_index( "integrated mechanical strain" );
	static int integrated_max_strain_index = liver_custom_data.find_variable_index( "max integrated mechanical strain" );

	if( pCell->custom_data[integrated_strain_index] > pCell->custom_data[integrated_max_strain_index] )
	{
		return true;
	}
	return false;
}

std::vector< std::string > liver_coloring_function( Cell* pCell )
{
	if( pCell->type == 0 ) // tumor cell
	{ return false_cell_coloring_live_dead( pCell ); }

	std::vector< std::string > output( 4, "rgb(217,128,64)" );
	// cyto_color, cyto_outline , nuclear_color, nuclear_outline

	if( pCell->type == 2 ) // portal triad
	{
		output.assign( 4, "rgb(128,128,128)" );
		return output;
	}

	if( pCell->type == 1 ) // central vein
	{
		output.assign( 4, "none" );
		return output;
	}

	bool is_apoptotic = false;
	if( pCell->phenotype.death.dead )
	{ is_apoptotic = true; }

	if( pCell->type == 3 )
	{
		if( is_apoptotic )
		{
			output.assign( 4, "black" );
		}
	}

	return output;
}

std::vector< std::string > liver_strain_coloring_function( Cell* pCell )
{
	static std::string str_liver_color = "rgb(217,128,64)";
	static std::vector<double> liver_color_vector = {217, 128, 64};
	// static std::vector<double> fuchsia = {255,0,128};
    static std::vector<double> fuchsia = {196,196,196};                      // in order to make the high pressure parenchyma with gray color

	static std::string str_yellow = "rgb(255,255,0)";
	std::vector< std::string > output( 4, str_liver_color );
/*
	if( pCell->type == 0 ) // tumor cell
	{ return false_cell_coloring_live_dead( pCell ); }
*/

	if( pCell->type == 0 ) // tumor cell
	{       
		//start with standard tumor cell coloring
		output = false_cell_coloring_Ki67( pCell );		
		// if tumor cell is dead, exits, as pressure check is unnecessary                     !!!!!!!!!!!!!!!!!!!!!!
		if ( pCell->phenotype.death.dead == true )
		{
			return output;                      
		}
		             
		//modify cell color based on pressure
		double pressure_test_value = 1.0;                                      // the value should be equal to p2
		if (pCell->state.simple_pressure > pressure_test_value)
		{
			output[0] = str_yellow;
			output[2] = str_yellow;                        
		}
                
	 	return output;
	}

	// cyto_color, cyto_outline , nuclear_color, nuclear_outline

	if( pCell->type == 2 ) // portal triad
	{
		output.assign( 4, "rgb(128,0,0)" );
		return output;
	}

	if( pCell->type == 1 ) // central vein
	{
		output.assign( 4, "none" );
		return output;
	}

	static int strain_index = pCell->custom_data.find_variable_index( "mechanical strain" );
	static int max_strain_index = pCell->custom_data.find_variable_index( "max mechanical strain" );

	double relative_strain = pCell->custom_data[strain_index] / ( pCell->custom_data[max_strain_index] + 1e-10 );
	double param = 1.0 - relative_strain;

	std::string color = "rgb(" + std::to_string( (int)round( relative_strain*fuchsia[0] + param*liver_color_vector[0] ) ) + "," +
		std::to_string( (int) round( relative_strain*fuchsia[1] + param*liver_color_vector[1] ) ) + "," +
		std::to_string( (int) round( relative_strain*fuchsia[2] + param*liver_color_vector[2] ) ) + ")";
/*
	std::string color = "rgb(";
		color << std::to_string( (int)round( param*liver_color_vector[0] ) ) << "," <<
		std::to_string( (int) round( param*liver_color_vector[1] ) ) << "," <<
		std::to_string( (int) round( param*liver_color_vector[2] ) ) << ")";
*/
	output.assign( 4 , color );
	
	// outstand the outline of the parenchyma with a relative_strain threshold                   !!!!!!!!!!!!!!!!!!!!
	// which is different from quiescent tumor cell 
	if( relative_strain > 0.01 )
	{
		output[3] = "black";  
	}

	bool is_apoptotic = false;
	if( pCell->phenotype.death.dead )
	{ is_apoptotic = true; }

	if( pCell->type == 3 )
	{
		if( is_apoptotic )
		{
			output.assign( 4, "black" );
		}
	}

	return output;
}

Cartesian_Mesh coarse_mesh;
std::vector<int> tumor_cell_counts;
std::vector<bool> apoptotic_parenchyma;
std::vector<bool> blocked_parenchyma;

void check_for_flow_disruption( void )
{
/*
	static Cartesian_Mesh coarse_mesh;
	static std::vector<int> tumor_cell_counts;
	static std::vector<bool> apoptotic_parenchyma;
	static std::vector<bool> blocked_parenchyma;
*/

	int coarsening_factor = 2; // 2 or 4 are best

	static bool setup_done = false;
	if( setup_done == false )
	{
		std::cout << "setup!" << std::endl;
		coarse_mesh.resize( microenvironment.mesh.bounding_box[0] , microenvironment.mesh.bounding_box[3],
			microenvironment.mesh.bounding_box[1] , microenvironment.mesh.bounding_box[4],
			coarsening_factor*microenvironment.mesh.bounding_box[2] , coarsening_factor*microenvironment.mesh.bounding_box[5],
			microenvironment.mesh.dx * coarsening_factor, microenvironment.mesh.dy * coarsening_factor , microenvironment.mesh.dz * coarsening_factor );

		tumor_cell_counts.assign( coarse_mesh.voxels.size() , 0 );
		apoptotic_parenchyma.assign( coarse_mesh.voxels.size() , false );
		blocked_parenchyma.assign( coarse_mesh.voxels.size() , false );

		setup_done = true;
	}

	// reset the tumor cell counts
	tumor_cell_counts.assign( coarse_mesh.voxels.size() , 0 );

	// count up the number of tumor cells in each coarsened voxel_index
	// also check for apoptotic parenchyma



	//Cell* pCell = NULL;
	#pragma omp parallel for
	for( int i=0; i < (*all_cells).size() ;i++ )
	{
		Cell* pCell = (*all_cells)[i];
		int n = coarse_mesh.nearest_voxel_index( pCell->position );
		if( pCell->type == 0 ) // tumor cell
		{
			blocked_parenchyma[n] = true;
		}
/*
		else // non-tumor
		{
			if( pCell->phenotype.death.dead == true )
			{
				apoptotic_parenchyma[n] = true;
			}
		}

*/
		if( pCell->type != 0 && pCell->phenotype.death.dead == true )           // not == !!!!!!!!!!!!!!!!!!!!!!!!!!
		{
			apoptotic_parenchyma[n] = true;
			blocked_parenchyma[n] = true;
		}
	}

	// non-isolated tumor cells disrupt flow.
	static int tumor_cell_disruption_threshold = 2;

	// now, go through the microenvironment, and remove Dirichlet nodes as needed

	#pragma omp parallel for
	for( int j=0; j < coarse_mesh.y_coordinates.size() ; j++ )
	{
		for( int i=0; i < coarse_mesh.x_coordinates.size() ; i++ )
		{
			int n = coarse_mesh.voxel_index(i,j,0);

			/* if( tumor_cell_counts[n] >= tumor_cell_disruption_threshold ||
				apoptotic_parenchyma[n] == true ) */
			if( blocked_parenchyma[n] == true )
			{

				// projecting on the finer microenvironment mesh
				for( int jj=coarsening_factor*j ; jj < coarsening_factor*(j+1) ; jj++ )
				{
					for( int ii=coarsening_factor*i ; ii < coarsening_factor*(i+1) ; ii++ )
					{
						int nn = microenvironment.mesh.voxel_index(ii,jj,0);
						microenvironment.mesh.voxels[nn].is_Dirichlet = false;
					}
				}

			}


		}
	}



	return;
}

void random_metastatic_seeding( double dt )
{
	static double seeding_rate_per_voxel = 0.1 *  1.0 / ( 125.0 * 24.0 * 60.0 );

	extern Cell_Definition HCT116;

	static int mets_count = 0;

	#pragma omp parallel for
	for( int n=0 ; n < microenvironment.number_of_voxels() ; n++ )
	{
		if( microenvironment.mesh.voxels[n].is_Dirichlet == true )
		{
			double probability = seeding_rate_per_voxel * dt;
			if( uniform_random() < probability )
			{
				#pragma omp critical
				{
					mets_count++;

					std::cout << "\n\n\t\tNEW METASTASIS! (number " << mets_count << ") \t\t at ";
					Cell* pCell = create_cell( HCT116 );
					pCell->assign_position( microenvironment.mesh.voxels[n].center );
					std::cout << pCell->position << "\n\n";
				}
			}
		}
	}

	return;
}

void advance_liver_model( double dt )
{
	static double elapsed_phenotype_time_from_last_update = 0.0;
	static double update_interval = phenotype_dt;

	if( elapsed_phenotype_time_from_last_update > update_interval )
	{
		check_for_flow_disruption();
		random_metastatic_seeding( dt );              // with random seeding    

		elapsed_phenotype_time_from_last_update = 0.0;
	}
	else
	{
		elapsed_phenotype_time_from_last_update += dt;
	}

	return;
}

void liver_plastoelastic_mechanics( Cell* pCell, Phenotype& phenotype , double dt )
{

	static int strain_index = pCell->custom_data.find_variable_index( "mechanical strain" );
	static int integrated_strain_index = pCell->custom_data.find_variable_index( "integrated mechanical strain" );


	static int spring_constant_index = pCell->custom_data.find_variable_index( "spring constant" );
	static int relaxation_constant_index = pCell->custom_data.find_variable_index( "mechanical relaxation rate" );

	// first, update the cell's velocity based upon the elastic model
	axpy( &( pCell->velocity ) , pCell->custom_data[spring_constant_index] , pCell->custom_data.vector_variables[1].value );

	// now, plastic mechanical relaxation

	static double plastic_temp_constant = -dt * pCell->custom_data[relaxation_constant_index];

	axpy( &(pCell->custom_data.vector_variables[0].value) , plastic_temp_constant , pCell->custom_data.vector_variables[1].value );

	return;
}

void HCT116_phenotype_update_function( Cell* pCell , Phenotype& phenotype, double dt )
{
	static int pressure_index = pCell->custom_data.find_variable_index("pressure");
	
	//std::cout<< "np=" << np << std::endl;
	//std::cout<< "pCell->state.simple_pressure=" << pCell->state.simple_pressure << std::endl;
	//std::cout << "var size" << pCell->custom_data.variables.size() << std::endl; 
	//std::cout << "var name: " << pCell->custom_data.variables[np].name << std::endl; 
	
    pCell->custom_data[pressure_index] = pCell->state.simple_pressure; //  save the pressure data
	
	
    // for now, put in the oxygen-based one.
    update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
	
	// add a check -- exit if the cell is dead 
	if( phenotype.death.dead == true )
	{ return; }
	
	// pressure-dependent phenotype
	static int start_phase_index = 0;
	static int end_phase_index = 0;
	double p1 = 0.0;
	double p2 = parameters.doubles( "tumor_max_pressue" );  //1.0 ;   
	
	double p = pCell->custom_data[pressure_index];
	double b_star = pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index);
	double multiplier = (p2 - p)/(p2 - p1);
	
	if (multiplier < 0.0)
	{
		multiplier = 0.0;
	}
	
	if (multiplier > 1.0)
	{
		multiplier = 1.0;
	}
	
	// phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) = multiplier*b_star; // birth_rate = b_star * pressure_stuff
	phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) *= multiplier; // birth_rate = b_star * oxygen_stuff * pressure_stuff 
	
	//std::cout<<"multiplier*b_star="<<multiplier*b_star<<std::endl;
	
	static int birth_rate_index = pCell->custom_data.find_variable_index("birth_rate");
    pCell->custom_data[birth_rate_index] = phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index); //  save the birth rate data

	
	return;
	
	
}
