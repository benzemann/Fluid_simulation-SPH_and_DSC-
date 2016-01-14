#pragma once
#include "util.h"

class SPH
{
public:
	SPH();
	/*
	Initializes SPH
	*/
	void init();
	/*
	Update all particles in SPH, the input delta time is the time since last update in ms
	*/
	void update(double delta_time);
	/*
	Calculate the density for the input particle
	*/
	void calculate_density();
private:
	/*
	Here all constants used in the SPH calculations
	*/

};

