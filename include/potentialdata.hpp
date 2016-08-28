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

#ifndef _POTENTIAL_DATA_H
#define _POTENTIAL_DATA_H

#include <stdlib.h>

#include "mkl.h"

#include "simulationdata.hpp"
#include "wavefunction.hpp"

class PotentialData {
public:
	PotentialData(SimulationData &sim_data);
	~PotentialData();

	//Arrays for holding potential information
	double *harmonic_trap;
	double *green_potential;
	double *kinetic_energy;
	double *nonlinear;
	MKL_Complex16 *pos_time_evolution;
	MKL_Complex16 *mom_time_evolution;

	//Allocation and assignment functions
	void calculate_non_linear(SimulationData &sim_data, WaveFunction &psi);
	void calculate_green(SimulationData &sim_data);
	void assign_position_time_evolution(SimulationData &sim_data, WaveFunction &psi, bool harmonic_on, bool green_on, bool is_real);
	void assign_momentum_time_evolution(SimulationData &sim_data, WaveFunction &psi, bool is_real);
};

#endif    //    _POTENTIAL_DATA_H
