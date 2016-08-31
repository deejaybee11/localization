/**The MIT License (MIT)
*
*Copyright (c) 2016 Dylan
*
*Permission is hereby granted, free of charge, to any person obtaining a copy
*of this software and associated documentation files (the "Software"), to deal
*in the Software without restriction, including without limitation the rights
*to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*copies of the Software, and to permit persons to whom the Software is
*furnished to do so, subject to the following conditions:
*
*The above copyright notice and this permission notice shall be included in all
*copies or substantial portions of the Software.
*
*THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
*SOFTWARE.
*/

#include "../include/solve.hpp"

#include <stdlib.h>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdexcept>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

#include "mkl.h"

#include "../include/simulationdata.hpp"
#include "../include/wavefunction.hpp"
#include "../include/potentialdata.hpp"
#include "../include/solve.hpp"
#include "../include/savedata.hpp"

void calculate_ground_state(SimulationData &sim_data, WaveFunction &psi, PotentialData &pot_data) {

	struct stat sb;
	if (stat("GroundState.fit", &sb) == 0){
		system("rm GroundState.fit");
		std::cout << "GroundState.fit deleted" << std::endl;
	}	
	//Create handle and descriptor for the FFT routine
	MKL_LONG status = 0;
	DFTI_DESCRIPTOR_HANDLE handle = 0;
	MKL_LONG N[2]; N[0] = sim_data.get_num_x(); N[1] = sim_data.get_num_y();
	status = DftiCreateDescriptor(&handle, DFTI_DOUBLE, DFTI_COMPLEX, 2, N);
	status = DftiSetValue(handle, DFTI_BACKWARD_SCALE, (1.0 / (N[0] * N[1])));
	status = DftiCommitDescriptor(handle);

	pot_data.assign_momentum_time_evolution(sim_data, psi, false);
	//Time to find the ground state!
	for (int i = 0; i < sim_data.num_imaginary_steps; ++i) {
		if (i%1000 == 0) {
			std::cout << "Imaginary step " << i << " of " << sim_data.num_imaginary_steps << std::endl;
		}
		
		//calculate abs(psi)**2
		psi.calc_abs_psi(sim_data.get_N());
		psi.calculate_norm(sim_data);
		//Calculate nonlinear energy term
		pot_data.calculate_non_linear(sim_data, psi);
		//Assign the position time evolution operator values
		pot_data.assign_position_time_evolution(sim_data, psi, true, false, false);
		//multiply by position operator with t/2
		vzMul(sim_data.get_N(), psi.psi, pot_data.pos_time_evolution, psi.psi);
 		//Perform forward transform and multiply by evolution operator	
		status = DftiComputeForward(handle, psi.psi, psi.psi);
		vzMul(sim_data.get_N(), psi.psi, pot_data.mom_time_evolution, psi.psi);
		//Transform back then multiply by other half of position operator
		status = DftiComputeBackward(handle, psi.psi, psi.psi);
	    	vzMul(sim_data.get_N(), psi.psi, pot_data.pos_time_evolution, psi.psi);
	}
	DftiFreeDescriptor(&handle);
	save_fits_image_wavefunction(sim_data, psi, "GroundState.fit");
}


void calculate_time_evolution(SimulationData &sim_data, WaveFunction &psi, PotentialData &pot_data) {

	struct stat sb;
	if (stat("FinalState.fit", &sb) == 0){
		system("rm FinalState.fit");
		std::cout << "FinalState.fit deleted" << std::endl;
	}	
	//Create handle and descriptor for the FFT routine
	MKL_LONG status = 0;
	DFTI_DESCRIPTOR_HANDLE handle = 0;
	MKL_LONG N[2]; N[0] = sim_data.get_num_x(); N[1] = sim_data.get_num_y();
	status = DftiCreateDescriptor(&handle, DFTI_DOUBLE, DFTI_COMPLEX, 2, N);
	status = DftiSetValue(handle, DFTI_BACKWARD_SCALE, (1.0 / (N[0] * N[1])));
	status = DftiCommitDescriptor(handle);

	pot_data.assign_momentum_time_evolution(sim_data, psi, true);
	//Time to find the ground state!
	for (int i = 0; i < sim_data.num_real_steps; ++i) {
		//calculate abs(psi)**2
		psi.calc_abs_psi(sim_data.get_N());
		if (i%10000 == 0) {
			sim_data.current_step = i;
			std::cout << "Real step " << i << " of " << sim_data.num_real_steps << std::endl;
			char buf1[200];
			char buf2[200];
			strcpy(buf1, sim_data.folder);
			sprintf(buf2, "/psi%d.fit\0", i/10000);
			strcat(buf1, buf2);
			int length = strlen(buf1);
			char *full_filename;
			full_filename = (char*)mkl_malloc(length*sizeof(char), 64);
			strcpy(full_filename, buf1);
			save_fits_image_wavefunction(sim_data, psi, full_filename);
			mkl_free(full_filename);
		}
		//Calculate nonlinear energy term
		pot_data.calculate_non_linear(sim_data, psi);
		//Assign the position time evolution operator values
		pot_data.assign_position_time_evolution(sim_data, psi, false, true, true);
		//multiply by position operator with t/2
		vzMul(sim_data.get_N(), psi.psi, pot_data.pos_time_evolution, psi.psi);
		//Perform forward transform and multiply by evolution operator	
		status = DftiComputeForward(handle, psi.psi, psi.psi);
		vzMul(sim_data.get_N(), psi.psi, pot_data.mom_time_evolution, psi.psi);
		//Transform back then multiply by other half of position operator
		status = DftiComputeBackward(handle, psi.psi, psi.psi);
		vzMul(sim_data.get_N(), psi.psi, pot_data.pos_time_evolution, psi.psi);
	}
	DftiFreeDescriptor(&handle);
	save_fits_image_wavefunction(sim_data, psi, "FinalState.fit");
}



