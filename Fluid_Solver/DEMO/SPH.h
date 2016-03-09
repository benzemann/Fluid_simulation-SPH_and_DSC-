#pragma once
#include "Particle_System.h"
#include "util.h"
#include "DSC.h"
#ifdef WIN32
#include <GL/glew.h>
#include <GL/glut.h>

#else
#include <GEL/GL/glew.h>
#include <GLUT/glut.h>
#include <GEL/GLGraphics/SOIL.h>
#endif

using namespace particlesystem;
//using namespace CGLA;
//using namespace DSC2D;

class SPH
{
	struct Collision_Box {
		Collision_Box(DSC2D::vec2 p, float h, float w) :
			pos(p),
			height(h),
			width(w),
			color(DSC2D::vec3(0.0,0.0,1.0))
		{
		}
		DSC2D::vec2 pos;
		float height;
		float width;
		DSC2D::vec3 color;
		vector<int> ghost_particles;
		void set_color(DSC2D::vec3 c) {
			color = c;
		}
		// Checks if position p is inside the collision box
		bool is_inside(DSC2D::vec2 p) {
			if (p[0] > (pos[0] - (width*0.5)) &&
				p[0] < (pos[0] + (width*0.5)) &&
				p[1] > (pos[1] - (height*0.5)) &&
				p[1] < (pos[1] + (height*0.5))
				) {
				return true;
			}
			return false;
		}
	};
public:
	SPH(int tmp);
	/*
	Initialization of SPH
	*/
	void init();
	/*
	Reset particles and collision boxes
	*/
	void reset();
	/*
	Update all particles in SPH, the input delta time is the time since last update in ms
	*/
	void update(double delta_time);
	/*
	Calculations for the different SPH components
	*/
	double calculate_density(Particle p, vector<Particle> close_particles);
	double calculate_local_pressure(Particle p);
	DSC2D::vec2 calculate_pressure(Particle p, vector<Particle> close_particles);
	DSC2D::vec2 calculate_viscocity(Particle p, vector<Particle> close_particles);
	DSC2D::vec2 calculate_external_forces(Particle p);
	DSC2D::vec2 calculate_surface_tension(Particle p, vector<Particle> close_particles);
	/*
	Find close particles to a given pos or particle
	*/
	vector<Particle> get_close_particles(Particle* p_ptr, double radius=-1);
	vector<Particle> get_close_particles_to_pos(DSC2D::vec2 pos, double radius = -1);
	Particle get_closest_particle(DSC2D::vec2 pos);
	/*
	Returns the number of particles in the system
	*/
	int get_no_of_particle();
	/*
	Returns position of particle pos
	*/
	DSC2D::vec2 get_particle_pos(int id);
	/*
	Drawing functions 
	*/
	void draw_particles();
	void draw_collision_boxes();
	void draw_kernel_radius();
	void draw_velocities();
	void draw_pressure();
	void draw_viscocity();
	void draw_surface_tension();
	void draw_dsc_velocities(DSC2D::DeformableSimplicialComplex& dsc);
	bool is_it_drawing_kernel() {
		return is_drawing_kernel;
	}
	bool is_it_drawing_velocities() {
		return is_drawing_velocities;
	}
	bool is_it_drawing_dsc_vel() {
		return is_drawing_dsc_vel;
	}
	/*
	Calculates new velocities and apply it to position
	*/
	void update_position(double delta_time);
	void update_velocity(double delta_time);
	/*
	Functions for creating or moving collision boxes and new particles
	*/
	void create_collision_box(DSC2D::vec2 pos, float height, float width, DSC2D::vec3 c=DSC2D::vec3(0.0,0.0,1.0));
	void create_collision_box_at_mouse_pos();
	void create_particle_at_mouse_pos();
	void move_collision_box(int id);
	/*
	Check if particle is colliding with a collision box.
	The particle are then moved to the point of collision and
	an impulse force is returned
	*/
	DSC2D::vec2 check_collision(Particle p, double delta_time);
	DSC2D::vec2 calculate_collision_impulse(Particle p, DSC2D::vec2 wall_normal);
	/*
	Signed distance from point p to line with a point on a line and a direction
	*/
	double signed_distance(DSC2D::vec2 p, DSC2D::vec2 plane_point, DSC2D::vec2 plane_dir);

	vector<double> signed_distance_to_walls(DSC2D::vec2 p, Collision_Box c_b);
	/*
	The three different smoothing kernels
	*/
	float poly6_kernel(float r, float d);

	DSC2D::vec2 gradient_kernel(float r, DSC2D::vec2 v, float d);

	float laplacian_kernel(float r, float d);
	/*
	Projecting of velocities from and to particles
	*/
	void project_velocities(DSC2D::DeformableSimplicialComplex& dsc);
	void project_velocities_to_particles(DSC2D::DeformableSimplicialComplex& dsc);
	HMesh::VertexAttributeVector<DSC2D::vec2> get_dsc_vert_velocities() {
		return dsc_vert_velocities;
	}
	void set_dsc_vert_velocities(HMesh::VertexAttributeVector<DSC2D::vec2> vels) {
		dsc_vert_velocities = vels;
	}
	DSC2D::vec2 get_dsc_vert_vel(HMesh::VertexID vi) {
		return dsc_vert_velocities[vi];
	}
	/*
	True if p is inside the triangle spanned by p0, p1, and p2
	*/
	bool in_triangle(DSC2D::vec2 p, DSC2D::vec2 p0, DSC2D::vec2 p1, DSC2D::vec2 p2);
	/* 
	Debugging helper functions
	*/
	vector<double*> get_user_variables_ptr();
	vector<bool*> get_user_flags_ptr() {
		vector<bool*> user_flags;
		user_flags.push_back(&is_drawing_kernel);
		user_flags.push_back(&is_drawing_velocities);
		user_flags.push_back(&is_drawing_viscocity);
		user_flags.push_back(&is_drawing_pressure);
		user_flags.push_back(&is_drawing_surface_tension);
		user_flags.push_back(&is_drawing_dsc_vel);
		user_flags.push_back(&enabled_grid);
		user_flags.push_back(&project_velocities_to_dsc);
		user_flags.push_back(&enabled_dsc_tracking);
		user_flags.push_back(&enabled_fluid_detection);
		return user_flags;
	}
	bool is_it_projecting() {
		return project_velocities_to_dsc;
	}
	bool is_dsc_tracking() {
		return enabled_dsc_tracking;
	}
	bool is_fluid_detection() {
		return enabled_fluid_detection;
	}
	double get_right_wall() {
		return collision_boxes[0].pos[0] - (collision_boxes[0].width/2.0);
	}


	DSC2D::vec2 calculate_inward_normal(DSC2D::vec2 pos);
private:
	vector<Collision_Box> collision_boxes;
	HMesh::VertexAttributeVector<DSC2D::vec2> dsc_vert_velocities;
	bool is_drawing_kernel = false;
	bool is_drawing_velocities = false;
	bool is_drawing_pressure = false;
	bool is_drawing_viscocity = false;
	bool is_drawing_surface_tension = false;
	bool is_drawing_dsc_vel = false;
	bool enabled_grid = true;
	bool project_velocities_to_dsc = false;
	bool enabled_dsc_tracking = false;
	bool enabled_fluid_detection = false;
	/*
	These are user controlled constant to get the desired fluid
	*/
	float MAX_VELOCITY = 1.0;
	double COEFFICIENT_OF_RESTITUTION = 0.0;//0.4
	double COEFFICIENT_OF_FRICTION = 0.03;//0.03
	double KERNEL_RADIUS = 30.0;
	double GAS_CONSTANT = 18900000;
	double REST_DENSITY = 0.000103;
	double VISCOCITY_TERM = 10.0; 
	
};





