#pragma once
#include <vector>
#include <CGLA/Vec3d.h>
#include "Particle_System.h"
typedef CGLA::Vec3d vec3;
using namespace std;
using namespace particlesystem;
class SPH
{
	struct Collision_Box {
		Collision_Box(vec3 p, float h, float w, float l) :
			pos(p),
			height(h),
			width(w),
			length(l)
		{
		}
		vec3 pos;
		float height;
		float width;
		float length;
		// Checks if position p is inside the collision box
		bool is_inside(vec3 p) {

			if (p[0] > (pos[0] - (width*0.5)) &&
				p[0] < (pos[0] + (width*0.5)) &&
				p[1] > (pos[1] - (length*0.5)) &&
				p[1] < (pos[1] + (length*0.5)) &&
				p[2] > (pos[2] - (height*0.5)) &&
				p[2] < (pos[2] + (height*0.5))
				) {
				return true;
			}
			return false;
		}
	};
public:
	SPH();
	/*
	Initialization of SPH
	*/
	void init();
	/*
	Reset particles and collision boxes
	*/
	void reset();
	/*
	Find close particles to a given pos or particle
	*/
	vector<Particle> get_close_particles(Particle* p_ptr, double radius = -1);
	vector<Particle> get_close_particles_to_pos(vec3 pos, double radius = -1);
	Particle get_closest_particle(vec3 pos);

	void update(double delta_time);

	double calculate_density(Particle p, vector<Particle> close_particles);
	double calculate_local_pressure(Particle p);
	vec3 calculate_pressure(Particle p, vector<Particle> close_particles);
	vec3 calculate_viscocity(Particle p, vector<Particle> close_particles);
	vec3 calculate_external_forces(Particle p);

	double calculate_a(Particle p, vector<Particle> close_particles);
	void correct_density_error(double delta_time);

	double poly6_kernel(float r, float d);
	vec3 gradient_kernel(float r, vec3 v, float d);
	float laplacian_kernel(float r, float d);

	void update_position(double delta_time);
	void update_velocity(double delta_time);

	void draw_particles(vec3 light_pos);
	void draw_collision_boxes(vec3 light_pos);

	void create_collision_box(vec3 pos, float height, float width, float length);
	vec3 check_collision(Particle p, double delta_time);
	vector<double> signed_distance_to_walls(vec3 p, Collision_Box c_b);
	double signed_distance(vec3 p, vec3 plane_point, vec3 plane_dir);
	vec3 calculate_collision_impulse(Particle p,vec3 wall_normal);

	void move_wall() {
		collision_boxes[2].pos[1] = 800;
	}

	double CFL_correction(double delta_time);

private:
	vector<Collision_Box> collision_boxes;
	double down_scale = 200.0;
	double avg_density = 0.0;
	/*
	These are user controlled constant to get the desired fluid
	*/
	float MAX_VELOCITY = 1.0;

	double COEFFICIENT_OF_RESTITUTION = 0.5;//0.4
	double COEFFICIENT_OF_FRICTION = 0.00;//0.03
	double KERNEL_RADIUS = 30.0;
	double GAS_CONSTANT = 18900000;
	double REST_DENSITY = 0.000103;
	double delta = 0.008;
	double VISCOCITY_TERM = 10.0;
};

