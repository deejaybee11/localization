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
#include <iomanip>
#include <random>

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
			if ((fabs(sim_data.x[i] - sim_data.x_offset) <= sim_data.channel_width/2.0) && ((sim_data.y[j] - sim_data.y_offset - sim_data.dumbell_radius) >= 0) && ((sim_data.y[j] - sim_data.y_offset) <= (sim_data.dumbell_radius + sim_data.channel_length))) {
				sim_data.num_pixels_in_channel_total += 1;
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

void PotentialData::smooth_edges_green(SimulationData &sim_data, int num_iterations) {
	double sigma = 7.0;
	double r, s = 2.0 * sigma * sigma;
	double sum_kernel = 0.0;

	struct stat sb;
	if (stat("SmoothedPotential.fit", &sb) == 0){
		system("rm SmoothedPotential.fit");
		std::cout << "SmoothedPotential.fit deleted" << std::endl;
	}	
	
	int gridsize = 11;
	int lower_bound = floor(gridsize/2);
	int upper_bound = lower_bound + 1;

	double gauss_kernel[gridsize*gridsize];
	double cosine_kernel[gridsize*gridsize];
	double omega = 0;
	double t = 0;
	int index;
	int index2;
	for (int i = -lower_bound; i < upper_bound; ++i) {
		for (int j = -lower_bound; j < upper_bound; ++j) {
			index = (i+lower_bound)*gridsize + (j+lower_bound);
			r = sqrt(i*i + j*j);
			gauss_kernel[index] = exp(-(r*r)/s)/(M_PI * s);
			cosine_kernel[index] = cos(omega * t);
			sum_kernel += gauss_kernel[index];
		}
	}

	for (int i = 0; i < gridsize*gridsize; ++i) {
		gauss_kernel[i] /= sum_kernel;
	}

	double *green_copy = 0;
	green_copy = (double*)mkl_malloc(sim_data.get_N() * sizeof(double), 64);
	for (int i = 0; i < sim_data.get_N(); ++i) {
		green_copy[i] = 0;
	}

	double mysum = 0;
	int x1, y1;
	for (int iteration = 0; iteration < num_iterations; ++iteration) {

		for (int i = 0; i < sim_data.get_num_x(); ++i) {
			for (int j = 0; j < sim_data.get_num_y(); ++j) {
				index2 = i*sim_data.get_num_x() + j;
				for (int k = -lower_bound; k < upper_bound; ++k) {
					for (int l = -lower_bound; l < upper_bound; ++l) {
						index = (k+lower_bound)*gridsize + (l+lower_bound);
						if ((i - k) < 0) {
							x1 = (i - k) + sim_data.get_num_x();
						}
						else if ((i-k) >= sim_data.get_num_x()){
							x1 = (i-k) - sim_data.get_num_x();
						}
						else {
							x1 = (i-k);
						}
						if ((j - l) < 0) {
							y1 = (j - l) + sim_data.get_num_y();
						}
						else if ((j-l) >= sim_data.get_num_y()){
							y1 = (j-l) - sim_data.get_num_y();
						}
						else {
							y1 = (j-l);
						}

						green_copy[index2] += gauss_kernel[index] * this->green_potential[x1*sim_data.get_num_y() + y1];
//						green_copy[index2] += gauss_kernel[index] * this->green_potential[(i - k)*sim_data.get_num_y() + (j-l)];
					}
				}
			}
		}
	}


	for (int i = 0; i < sim_data.get_N(); ++i) {
		this->green_potential[i] = green_copy[i];
	}

	
	std::random_device rnd;
	std::mt19937 mt(rnd());
	std::uniform_real_distribution<double> dist(0.0,1.0);

	int num_pixels_in_channel = floor(sim_data.num_pixels_in_channel_total * sim_data.fill_factor);
	int num_pixels_left = floor(sim_data.num_pixels_in_channel_total * sim_data.fill_factor);

	int xvalue = 0;
	int yvalue = 0;
	index = 0;

	while (num_pixels_left >= 1) {
		
		for (int i = 0; i < sim_data.get_num_x(); ++i) {
			for (int j = 0; j < sim_data.get_num_y(); ++j) {
				index = i * sim_data.get_num_y() + j;
				if ((fabs(sim_data.x[i] - sim_data.x_offset) <= sim_data.channel_width/2.0) && ((sim_data.y[j] - sim_data.y_offset - sim_data.dumbell_radius) >= 0) && ((sim_data.y[j] - sim_data.y_offset) <= (sim_data.dumbell_radius + sim_data.channel_length))) {// && ((sim_data.y[i] - sim_data.y_offset) >= sim_data.dumbell_radius))  {
					if (dist(mt) < sim_data.fill_factor) {
						this->green_potential[index] = 5000;
						this->green_potential[index+1] = 2500;
						this->green_potential[index-1] = 2500;
						this->green_potential[index+sim_data.get_num_y()] = 2500;
						this->green_potential[index-sim_data.get_num_y()] = 2500;
						num_pixels_left -= 1;
						if (num_pixels_left == 0) {
							goto BREAKLOOP;
						}
					}
				}
			}
		}
	}

	BREAKLOOP:
		printf("While broken\n");
	
	save_fits_image_potential(sim_data, this->green_potential, "SmoothedPotential.fit");	

}


PotentialData::~PotentialData() {
	mkl_free(harmonic_trap);
	mkl_free(nonlinear);
	mkl_free(green_potential);
	mkl_free(mom_time_evolution);
	mkl_free(pos_time_evolution);
}
	
