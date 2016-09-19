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

#include "../include/simulationdata.hpp"

#include <stdlib.h>
#include <iostream>
#include <sys/stat.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "../inih/ini.h"
#include "../inih/INIReader.h"


//Class Constructor

SimulationData::SimulationData(int num_x, int num_y) {

	//Builds folder with name as current date
	struct stat sb;

	if (!(stat("fits", &sb) == 0 && S_ISDIR(sb.st_mode))){
		system("mkdir fits");
	}	
	time(&this->curr_time);
	this->mytime = localtime(&this->curr_time);
	sprintf(this->command, "mkdir fits/%.4d%.2d%.2d%.2d%.2d", 1900 + this->mytime->tm_year, 1 + this->mytime->tm_mon, this->mytime->tm_mday, this->mytime->tm_hour, this->mytime->tm_min);
	system(this->command);
	printf("Directory Created\n");
	sprintf(this->folder, "fits/%.4d%.2d%.2d%.2d%.2d", 1900 + this->mytime->tm_year, 1 + this->mytime->tm_mon, this->mytime->tm_mday, this->mytime->tm_hour, this->mytime->tm_min);
	sprintf(this->date, "%.4d%.2d%.2d%.2d%.2d", 1900 + this->mytime->tm_year, 1 + this->mytime->tm_mon, this->mytime->tm_mday, this->mytime->tm_hour, this->mytime->tm_min);

	//Loads INI file
	INIReader reader("conf.ini");
	if (reader.ParseError() <0) {
		std::cout << "Error loading config file 'conf.ini', using default parameters." << std::endl;
	}

	//Simulation Data
	this->num_x = num_x;
	this->num_y = num_y;
	this->N = num_x*num_y;
	//Length scales
	this->length_x = reader.GetInteger("config", "length_x", 150);
	std::cout << "The ini file has length_x = " << reader.GetInteger("config", "length_x", 150) << " but the value in memory is = " << this->length_x  << std::endl;
	this->length_y = reader.GetInteger("config", "length_y", 150);
	//BEC parameters
	this->sigma_x = reader.GetReal("config", "sigma_x", 1);
	this->sigma_y = reader.GetReal("config", "sigma_y", 1.2);
	this->beta = reader.GetReal("config", "beta", 1);
	//Harmonic Trap
	this->gamma_x = reader.GetReal("config", "gamma_x", 1);
	this->gamma_y = reader.GetReal("config", "gamma_y", 1.2);
	//Green parameters
	this->fill_factor = reader.GetReal("config", "fill_factor", 0.001);
	this->scatter_height = reader.GetReal("config", "scatter_height", 5000);
	this->dumbell_radius = reader.GetReal("config", "dumbell_radius", 15);
	this->channel_width = reader.GetReal("config", "channel_width", 10);
	this->channel_length = reader.GetReal("config", "channel_length", 25);
	this->num_pixels_in_channel_total = 0;
	this->x_offset = reader.GetReal("config", "x_offset", 0);
	this->y_offset = reader.GetReal("config", "y_offset", 0);
	//Array memory
	this->x = (double*)mkl_malloc(this->num_x * sizeof(double), 64);
	this->y = (double*)mkl_malloc(this->num_y * sizeof(double), 64);
	this->px = (double*)mkl_malloc(this->num_x * sizeof(double), 64);
	this->py = (double*)mkl_malloc(this->num_y * sizeof(double), 64);
	//Populate position arrays
	for (int i = 0; i < this->num_x; ++i) {
		this->x[i] = -0.5 * this->length_x + i * this->length_x / ((double)this->num_x);
	}
	for (int i = 0; i < this->num_y; ++i) {
		this->y[i] = -0.5 * this->length_y + i * this->length_y / ((double)this->num_y);
	}
	this->dx = this->x[1] - this->x[0];
	this->dy = this->y[1] - this->y[0];
	//Tine steps
	this->dt = this->dx * 0.001;
	this->num_imaginary_steps = reader.GetInteger("config", "num_imaginary_steps", 100000);
	this->num_real_steps = reader.GetInteger("config", "num_real_steps", 100000000);
	//Populate momentum arrays
	double ax = -0.5 * this->num_x;
	double bx = 0.5 * this->num_x - 1.0;
	double ay = -0.5 * this->num_y;
	double by = 0.5 * this->num_y - 1.0;
	double step_x = (2 * M_PI / this->length_x) * ((bx - ax) / (this->num_x - 1.0));
	double step_y = (2 * M_PI / this->length_y) * ((by - ay) / (this->num_y - 1.0));
	for (int i = 0; i < this->num_x; ++i) {
		this->px[i] = (2 * M_PI / this->length_x) * ax + i * step_x;
	}
	for (int i = 0; i < this->num_y; ++i) {
		this->py[i] = (2 * M_PI / this->length_y) * ax + i * step_y;
	}
	//Perform fftshift of momentum arrays
	double temp;
	int n2x = this->get_num_x() / 2.0;
	int n2y = this->get_num_y() / 2.0;
	for (int i = 0; i < n2x; ++i) {
		temp = this->px[i];
		this->px[i] = this->px[i + n2x];
		this->px[i + n2x] = temp;
	}
	for (int i = 0; i < n2y; ++i) {
		temp = this->py[i];
		this->py[i] = this->py[i + n2y];
		this->py[i + n2y] = temp;
	}

};

//Class Destructor
SimulationData::~SimulationData() {
	mkl_free(x);
	mkl_free(y);
	mkl_free(px);
	mkl_free(py);
};
