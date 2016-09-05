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

#include "mkl.h"

#include "../include/savedata.hpp"

void save_fits_image_wavefunction(SimulationData &sim_data, WaveFunction &psi, const char *fits_file_name) {
	//Stuff needed for fits files
	double *save_data;
	save_data = (double*)mkl_malloc(sim_data.get_N() * sizeof(double), 64);
	fitsfile *fptr;
	int status = 0;
	long fpixel = 1, naxis = 2, nelements;
	long naxes[2] = {sim_data.get_num_y(), sim_data.get_num_x()};
	//Absolute value of psi calculated
	psi.calc_abs_psi(sim_data.get_N());
	//Add psi values to savedata
	#pragma omp parallel for
	for (int i = 0; i < sim_data.get_N(); ++i) {
		save_data[i] = psi.abs_psi[i];
	}
	int numx = sim_data.get_num_x();
	int numy = sim_data.get_num_y();
	int step = sim_data.current_step;
	char *date = sim_data.date;
	double dumbell_radius = sim_data.dumbell_radius;
	double channel_length = sim_data.channel_length;
	double channel_width = sim_data.channel_width;
	double fill_factor = sim_data.fill_factor;
	//Create fits file
	fits_create_file(&fptr, fits_file_name, &status);
	fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);
	nelements = naxes[0] * naxes[1];
	fits_update_key(fptr, TINT, "NUMX", &numx, "Number of points in x", &status);
	fits_update_key(fptr, TINT, "NUMY", &numy, "Number of points in y", &status);
	fits_update_key(fptr, TINT, "TSTEP", &step, "Current time step.", &status);
	fits_update_key(fptr, TSTRING, "DATE", &date, "Date", &status);
	fits_update_key(fptr, TDOUBLE, "DRADUIS", &dumbell_radius, "Radius of dumbell potential", &status);
	fits_update_key(fptr, TDOUBLE, "CWIDTH", &channel_width, "Width of channel", &status);
	fits_update_key(fptr, TDOUBLE, "CLENGTH", &channel_length, "Length of channel", &status);
	fits_update_key(fptr, TDOUBLE, "FILLFAC", &fill_factor, "Fill factor", &status);
	fits_write_img(fptr, TDOUBLE, fpixel, nelements, save_data, &status);
	fits_close_file(fptr, &status);
	fits_report_error(stderr, status);

	mkl_free(save_data);

}


void save_fits_image_potential(SimulationData &sim_data, double *potential, const char *fits_file_name) {
	//Stuff needed for fits files
	double *save_data;
	save_data = (double*)mkl_malloc(sim_data.get_N() * sizeof(double), 64);
	fitsfile *fptr;
	int status = 0;
	long fpixel = 1, naxis = 2, nelements;
	long naxes[2] = {sim_data.get_num_y(), sim_data.get_num_x()};
	//add values to save_data
	#pragma omp parallel for
	for (int i = 0; i < sim_data.get_N(); ++i) {
		save_data[i] = potential[i];
	}
	//Create fits file
	fits_create_file(&fptr, fits_file_name, &status);
	fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);
	nelements = naxes[0] * naxes[1];
	fits_write_img(fptr, TDOUBLE, fpixel, nelements, save_data, &status);
	fits_close_file(fptr, &status);
	fits_report_error(stderr, status);

	mkl_free(save_data);

}
