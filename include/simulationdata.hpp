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


#ifndef _LOCALIZATION_SIM_DATA_H
#define _LOCALIZATION_SIM_DATA_H

#include <stdlib.h>
#include <iostream>

#include "mkl.h"

/**
 * SimulationData class
 *
 * Stores data and arrays used in the simulation
 */

class SimulationData {
public:

	SimulationData(int num_x, int num_y);
	~SimulationData();

	//Time and folder to save data
	time_t curr_time;
	struct tm *mytime;
	char folder[80];
	char command[80];
	//Array pointers 
	double *x;
	double *y;
	double *px;
	double *py;
	//Length scales
	double length_x;
	double length_y;
	double dx;
	double dy;
	double dpx;
	double dpy;
	//BEC parameters
	double sigma_x;
	double sigma_y;
	double beta;
	//Harmonic potential
	double gamma_x;
	double gamma_y;
	//Green potential
	double fill_factor;
	double scatter_height;
	double dumbell_radius;
	double channel_width;
	double channel_length;
	double x_offset;
	double y_offset;
	//Iteration stuff
	int num_real_steps;
	int num_imaginary_steps;
	const char *R = "REAL";
	const char *I = "IMAG";
	//Private setters and getters
	int get_num_x() { return this->num_x; };
	int get_num_y() { return this->num_y; };
	int get_N() { return this->N; };
	double get_dt() { return this->dt; };

private:
	//Number of points
	int num_x;
	int num_y;
	int N;
	//Time step
	double dt;

};

#endif    //    _LOCALIZATION_SIM_DATA_H_
