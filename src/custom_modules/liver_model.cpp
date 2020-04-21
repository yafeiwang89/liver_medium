/*
#############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the ver-  #
# sion number, such as below:                                               #
#                                                                           #
# We implemented and solved the model using PhysiCell (Version 0.0.5) [1].  #
#                                                                           #
# [1] A Ghaffarizadeh, SH Friedman, and P Macklin, PhysiCell: an open       #
#    source physics-based simulator for multicellular systemssimulator, 	#
#	 J. Comput. Biol., 2016 (submitted). 									# 
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

#include "./liver_model.h"

// std::vector<int> dirichlet_voxels_to_remove; 
std::vector<bool> dirichlet_voxels_to_remove; 

void update_liver( double dt )
{
	if( dirichlet_voxels_to_remove.size() != M.mesh.voxels.size() )
	{ dirichlet_voxels_to_remove.resize( M.mesh.voxels.size() , false ); }

	// remove Dirichlet conditions in disrupted regions 
	static int total_removed = 0; 
	int remove_count = 0; 
/*	
	for( int i=0; i < dirichlet_voxels_to_remove.size() ; i++ )
	{
		int j = dirichlet_voxels_to_remove[i]; 
		if( M.mesh.voxels[j].is_Dirichlet == true )
		{
			M.remove_dirichlet_node( j ); 			
			remove_count++; 
		}
	}
*/	
	for( int i=0; i < dirichlet_voxels_to_remove.size() ; i++ )
	{
		if( M.mesh.voxels[i].is_Dirichlet == true && dirichlet_voxels_to_remove[i] == true )
		{
			M.remove_dirichlet_node( i ); 			
			remove_count++; 
		}
	}
	total_removed += remove_count; 
	if( remove_count > 0 )
	{
	//	std::cout << "removed " << remove_count << " dirichlet nodes (" << total_removed << " lifetime total)" << std::endl; 
	}
	dirichlet_voxels_to_remove.resize( M.mesh.voxels.size() , false ); 

/*	
	#pragma omp parallel for 
	for( int i=0; i < (*all_cells).size() ; i++ )
	{
//		if( (*all_cells)[i]->custom_data.cell_type == 3 )
//		{ parenchyme_cell_rule( (*all_cells)[i] , dt ); }
		(*all_cells)[i]->custom_cell_rule( (*all_cells)[i] , dt ) ; 
	}
*/
	
	return; 
}


void parenchyme_cell_rule( Cell* pCell, double dt )
{
	// apoptosis rule 
	
	parenchyme_apoptosis_model1( pCell, dt ); 
	// parenchyme_apoptosis_model2( pCell, dt ); 
	// parenchyme_apoptosis_model3( pCell, dt ); 

	return; // if the code below is left on, PhysiCell crashes very early 
	
	// ellasticity 
	// (d/dt) x = r*( ECM_attach - x )

	// backwards Euler 
	std::vector<double> vec1 = pCell->position; 
	
	double const1 = pCell->custom_data.spring_constant; 
	const1 *= dt; // dt*spring_constant;
	double const2 = 1.0; 
	const2 += const1; // 1 + dt*spring_constant

	// y = y + a.*x
	// void axpy( std::vector<double>* y, std::vector<double>& a , std::vector<double>& x ); 	
	axpy( &vec1 , const1, 	pCell->custom_data.ECM_attach ); // x = x + const1*ECM_attach 
	vec1 /= const2; // backward Euler complete 
	pCell->assign_position( vec1 ); // safely assign this new position 
	
	// plasticity 
	// (d/dt) ECM_attach = r*( x - ECM_attach )	

	// backwards Euler
	const1 = pCell->custom_data.mechanical_relaxation_rate; 
	const1 *= dt; // dt*mechanical_relaxation_rate;
	const2 = 1.0; 
	const2 += const1; // 1 + dt*mechanical_relaxation_rate
	
	// y = y + a.*x
	// void axpy( std::vector<double>* y, std::vector<double>& a , std::vector<double>& x ); 	
	axpy( &(pCell->custom_data.ECM_attach) , const1, 	pCell->position ); // ECM_attach = ECM_attach + const1*x 
	pCell->custom_data.ECM_attach /= const2; // backward Euler complete 
	
	vec1 -= pCell->custom_data.ECM_attach; 
	pCell->custom_data.current_mechanical_strain = norm( vec1 ); 
	
	return; 
}

void parenchyme_apoptosis_model1( Cell* pCell, double)
{
	// if already apoptotic, don't bother 
	if( pCell->phenotype.cycle.phases[pCell->phenotype.current_phase_index].code == PhysiCell_constants::PhysiCell_constants::apoptotic )
	{ return; }
	
	double rn = uniform_random(); 
	if( rn < 5e-5 ) // 1e-5 ) // replace with correct condition 
	{
		// std::cout << "death!" << std::endl; 
		// set the current phase to apoptotic 
		pCell->phenotype.set_current_phase(PhysiCell_constants::apoptotic); 
		//pCell->phenotype.phase_model_initialized = false; 
		pCell->advance_cell_current_phase = death_apoptosis_model_stochastic; 
		pCell->phenotype.cycle.phases[pCell->phenotype.current_phase_index].elapsed_time = 0.0; 
		
		// if apoptotic, flag the node for removing the far-field Dirichlet condition 
//		int n = pCell->get_nearest_voxel_index();
		int n = pCell->get_current_voxel_index();
		// dirichlet_voxels_to_remove.push_back( n );
		if(n!=-1)
			dirichlet_voxels_to_remove[n] = true; 
	}	
	
	return; 
}

void parenchyme_apoptosis_model2( Cell* pCell, double)
{
	
	
	
}

void parenchyme_apoptosis_model3( Cell* pCell, double)
{
	
	
	
	
}
