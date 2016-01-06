#pragma once
#include "velocity_function.h"

class track_particle_function : public DSC2D::VelocityFunc<>
{
public:
	track_particle_function(double velocity, double accuracy) : VelocityFunc(velocity / 100., accuracy / 100.)
	{
		
	}

	/**
	Returns the name of the velocity function.
	*/
	virtual std::string get_name() const
	{
		return std::string("PARTICLE TRACK FUNCTION");
	}
	/**
	Computes the motion of each interface vertex and stores the destination in the simplicial complex class.
	*/
	virtual void deform(DSC2D::DeformableSimplicialComplex& dsc);

};

