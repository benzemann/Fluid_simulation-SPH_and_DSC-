#pragma once

#include "velocity_function.h"

class track_particles_function : public DSC::VelocityFunc<>
{
public:
	track_particles_function(real velocity, real accuracy, int max_time_steps = 500) :
		VelocityFunc<>(velocity, accuracy, max_time_steps)
	{

	}

	virtual std::string get_name() const
	{
		return std::string("TRACK PARTICLES MOTION");
	}

	/**
	Computes the motion of each interface vertex and stores the new position in new_pos in the simplicial complex class.
	*/
	virtual void deform(DSC::DeformableSimplicialComplex<>& dsc);
};

