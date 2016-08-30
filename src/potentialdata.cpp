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

#include "../include/potentialdata.hpp"

#include <stdlib.h>
#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <sys/stat.h>

#include "mkl.h"

#include "../include/simulationdata.hpp"
#include "../include/wavefunction.hpp"
#include "../include/savedata.hpp"

PotentialData::PotentialData(SimulationData &sim_data) {


	//Allocate memory arrays
	this->harmonic_trap = (double*)mkl_malloc(sim_data.get_N() * sizeof(double), 64);
	this->nonlinear = (double*)mkl_malloc(sim_data.get_N() * sizeof(double), 64);
	this->kinetic_energy = (double*)mkl_malloc(sim_data.get_N() * sizeof(double), 64);
	this->pos_time_evolution = (MKL_Complex16*)mkl_malloc(sim_data.get_N() * sizeof(MKL_Complex16), 64);
	this->mom_time_evolution = (MKL_Complex16*)mkl_malloc(sim_data.get_N() * sizeof(MKL_Complex16), 64);

	//Populate arrays with numbers
	double h_pot_val = 0;
	double k_en_val = 0;
	int index;

	#pragma omp parallel for private(index, h_pot_val)
	for (int i = 0; i < sim_data.get_num_x(); ++i) {
		for (int j = 0; j < sim_data.get_num_y(); ++j) {
			index = i * sim_data.get_num_y() + j;
			h_pot_val = 0.5 * (pow(sim_data.gamma_x, 2.0) * pow(sim_data.x[i], 2.0) + pow(sim_data.gamma_y, 2.0) * pow(sim_data.y[j], 2.0));
			k_en_val = 0.5 * (pow(sim_data.px[i], 2.0) + pow(sim_data.py[j], 2.0));
			this->harmonic_trap[index] = h_pot_val;
			this->kinetic_energy[index] = k_en_val;
		}
	}
	std::cout << this->harmonic_trap[120] << std::endl;
};

void PotentialData::calculate_green(SimulationData &sim_data) {

	struct stat sb;
	if (stat("GreenPotential.fit", &sb) == 0){
		system("rm GreenPotential.fit");
		std::cout << "GreenPotential.fit deleted" << std::endl;
	}	


	double radius_squared = 0;
	double radius_squared2 = 0;
	double dumbell_radius_squared = pow(sim_data.dumbell_radius, 2.0);
	this->green_potential = 0;
	this->green_potential = (double*)mkl_malloc(sim_data.get_N() * sizeof(double), 64);


	for (int i = 0; i < sim_data.get_N(); ++i) {
		this->green_potential[i] = 5000;
	}

	int index;
	for (int i = 0; i < sim_data.get_num_x(); ++i) {
		for (int j = 0; j < sim_data.get_num_y(); ++j) {
			index = i*sim_data.get_num_y() + j;
			radius_squared = pow(sim_data.x[i] - sim_data.x_offset, 2.0) + pow(sim_data.y[j] - sim_data.y_offset, 2.0);
			if (radius_squared <= dumbell_radius_squared) {
				this->green_potential[index] = 0;
			}
			radius_squared2 = pow(sim_data.x[i] - sim_data.x_offset, 2.0) + pow(sim_data.y[j] - (sim_data.y_offset + sim_data.channel_length + 2 * sim_data.dumbell_radius), 2.0);
			if (radius_squared2 <= dumbell_radius_squared) {
				this->green_potential[index] = 0;
			}
			if ((fabs(sim_data.x[i] - sim_data.x_offset) <= sim_data.channel_width/2.0) && ((sim_data.y[j] - sim_data.y_offset) >= 0) && ((sim_data.y[j] - sim_data.y_offset) <= (2*sim_data.dumbell_radius + sim_data.channel_length))) {
				this->green_potential[index] = 0;
			}
		}
	}
	save_fits_image_potential(sim_data, this->green_potential, "GreenPotential.fit");	
}

void PotentialData::assign_position_time_evolution(SimulationData &sim_data, WaveFunction &psi, bool harmonic_on, bool green_on, bool is_real){
	double theta = 0;
	double harm = 0;
	double green = 0;

	if (harmonic_on) harm = 1;
	if (green_on) green = 1;

	if (is_real) {
		#pragma omp parallel for private(theta)
		for (int i = 0; i < sim_data.get_N(); ++i) {
			theta = (this->nonlinear[i] + harm * this->harmonic_trap[i] + green * this->green_potential[i]) * 0.5 * sim_data.get_dt();
			this->pos_time_evolution[i].real = cos(theta);
			this->pos_time_evolution[i].imag = -1.0 * sin(theta);
		}
	}
	else {
		#pragma omp parallel for private(theta)
		for (int i = 0; i < sim_data.get_N(); ++i) {
			theta = (this->nonlinear[i] + harm * this->harmonic_trap[i]) * 0.5 * sim_data.get_dt();
			this->pos_time_evolution[i].real = exp(-1.0 * theta);
			this->pos_time_evolution[i].imag = 0;
		}
	}
}

void PotentialData::assign_momentum_time_evolution(SimulationData &sim_data, WaveFunction &psi, bool is_real){
	double theta = 0;
	if (is_real) {
		#pragma omp parallel for private(theta)
		for (int i = 0; i < sim_data.get_N(); ++i) {
			theta = this->kinetic_energy[i] * sim_data.get_dt();
			this->mom_time_evolution[i].real = cos(theta);
			this->mom_time_evolution[i].imag = -1.0 * sin(theta);
		}
	}
	else {
		#pragma omp parallel for private(theta)
		for (int i = 0; i < sim_data.get_N(); ++i) {
			theta = this->kinetic_energy[i] * sim_data.get_dt();
			this->mom_time_evolution[i].real = exp(-1.0 * theta);
			this->mom_time_evolution[i].imag = 0;
		}
	}
}

void PotentialData::calculate_non_linear(SimulationData &sim_data, WaveFunction &psi){
	double nonlinearval = 0;
	#pragma omp parallel for private(nonlinearval)
	for (int i = 0; i < sim_data.get_N(); ++i) {
		nonlinearval = sim_data.beta * psi.abs_psi[i];
		this->nonlinear[i] = nonlinearval;
	}

}


PotentialData::~PotentialData() {
	mkl_free(harmonic_trap);
	mkl_free(nonlinear);
	mkl_free(green_potential);
	mkl_free(mom_time_evolution);
	mkl_free(pos_time_evolution);
}
	
