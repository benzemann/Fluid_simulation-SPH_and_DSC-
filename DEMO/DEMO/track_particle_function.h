#pragma once
#include "velocity_function.h"
#include "ParticleSystem.h"

class track_particle_function : public DSC2D::VelocityFunc<>
{
public:
	track_particle_function(double velocity, double accuracy,DSC2D::ParticleSystem& p_sys ) : VelocityFunc(velocity / 100., accuracy / 100.), particle_system(p_sys)
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
	void deform_old(DSC2D::DeformableSimplicialComplex& dsc);
	virtual void deform(DSC2D::DeformableSimplicialComplex& dsc);

	virtual void calculate_forces(DSC2D::vec2 pos1, DSC2D::vec2 pos2, DSC2D::vec2 *force1, DSC2D::vec2 *force2, DSC2D::DeformableSimplicialComplex& dsc, HMesh::VertexID vert1, HMesh::VertexID vert2);

	virtual void track_particle_function::solve_for_incompressibility(std::vector<HMesh::FaceID> face_ids, std::vector<HMesh::VertexID> verts_ids, HMesh::VertexAttributeVector<DSC2D::vec2> &velocities,DSC2D::DeformableSimplicialComplex& dsc);

	void init();
public:
	// Tuan
	DSC2D::DeformableSimplicialComplex * dsc_ptr;

	// Constant volume
	double m_V0;

	// Boundary
	DSC2D::vec2 center_bound;
	double r_bound;

	double get_curent_volume();
	bool is_outside(DSC2D::vec2 pt, DSC2D::vec2 &projectionPoint);
	bool is_on_boundary(DSC2D::vec2 pt, DSC2D::vec2 &projectionPoint);
private:
	int solver_delay = 100;
	float last_area = 0.0f;
	float max_diff = 0.0f;
	float start_area = 0.0f;
	int DSC_iterations = 0;
	float r;
	float alpha;
	float mean_density;
	int no_segments;
	DSC2D::ParticleSystem& particle_system;
};

