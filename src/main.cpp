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

#include <stdlib.h>
#include <iostream>

#include "mkl.h"

#include "../include/simulationdata.hpp"
#include "../include/potentialdata.hpp"
#include "../include/wavefunction.hpp"
#include "../include/solve.hpp"

#if !defined(MKL_ILP64)
	#define LI "%li"
	#else
	#define LI "%lli"
#endif

int main() {
	
	//Sets up the multi-threading aspects
	putenv("KMP_BLOCKTIME=0");
	putenv("KMP_AFFINITY=verbose,granularity=fine,compact,norespect");
	system("echo KMP_BLOCKTIME = $KMP_BLOCKTIME");
	system("echo KMP_AFFINITY = $KMP_AFFINITY");
	mkl_set_num_threads(mkl_get_max_threads());
	mkl_disable_fast_mm();

	//Initiate classes
	SimulationData sim_data(1024, 1024);
	PotentialData pot_data(sim_data);
	WaveFunction psi(sim_data, pot_data.harmonic_trap);
	pot_data.calculate_green(sim_data);

	//Calculate the ground state
	calculate_ground_state(sim_data, psi, pot_data);
	calculate_time_evolution(sim_data, psi, pot_data);

	return 0;
}


