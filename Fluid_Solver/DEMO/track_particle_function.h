#pragma once
#include "velocity_function.h"
#include "SPH.h"
//#include "DSC.h"

class track_particle_function : public DSC2D::VelocityFunc<>
{
public:
	track_particle_function(double velocity, double accuracy, SPH& sph) : VelocityFunc(velocity / 100., accuracy / 100.), sph(sph)
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
	virtual void init(DSC2D::DeformableSimplicialComplex& dsc);
	virtual void deform_old(DSC2D::DeformableSimplicialComplex& dsc);
	virtual void deform_old_old(DSC2D::DeformableSimplicialComplex& dsc);
	virtual void deform(DSC2D::DeformableSimplicialComplex& dsc);
	virtual void done(DSC2D::DeformableSimplicialComplex& dsc);
	void solve_for_incompressibility(std::vector<HMesh::FaceID> face_ids,
		std::vector<HMesh::VertexID> verts_ids,
		HMesh::VertexAttributeVector<DSC2D::vec2> &velocities,
		DSC2D::DeformableSimplicialComplex& dsc);
	double get_curent_volume(DSC2D::DeformableSimplicialComplex& dsc);
	void solve_for_incompressibility_test(HMesh::VertexAttributeVector<DSC2D::vec2> &velocities, DSC2D::DeformableSimplicialComplex& dsc);
private:
	SPH& sph;
	vector<int> p_verts;
	// for debugging
	vector<double> vol_log;
	int tmp = 0;
	double m_V0 = 90000.0;
};

