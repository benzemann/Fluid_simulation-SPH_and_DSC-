#pragma once
#include "Particle_System.h"
#include "util.h"
#ifdef WIN32
#include <GL/glew.h>
#include <GL/glut.h>

#else
#include <GEL/GL/glew.h>
#include <GLUT/glut.h>
#include <GEL/GLGraphics/SOIL.h>
#endif

using namespace particlesystem;
using namespace CGLA;

class SPH
{
	struct Collision_Box {
		Collision_Box(vec2 p, float h, float w) :
			pos(p),
			height(h),
			width(w),
			color(vec3(0.0,0.0,1.0))
		{}
		vec2 pos;
		float height;
		float width;
		vec3 color;
		void set_color(vec3 c) {
			color = c;
		}
		// Checks if position p is inside the collision box
		bool is_inside(vec2 p) {
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
	void reset();
	/*
	Update all particles in SPH, the input delta time is the time since last update in ms
	*/
	void update(double delta_time);
	/*
	Calculate the density for the input particle
	*/
	double calculate_density(Particle p, vector<Particle> close_particles);
	double calculate_local_pressure(Particle p);
	vec2 calculate_pressure(Particle p, vector<Particle> close_particles);
	vec2 calculate_viscocity(Particle p, vector<Particle> close_particles);
	vec2 calculate_external_forces(Particle p);
	/*
	Draws all particles as red dots
	*/
	void draw_particles();
	/*
	Calculates new velocities and apply it to position
	*/
	void update_position(double delta_time);
	void update_velocity(double delta_time);
	/*
	Creates a collision box at the given postion with height and width
	*/
	void create_collision_box(vec2 pos, float height, float width, vec3 c=vec3(0.0,0.0,1.0));
	void create_collision_box_at_mouse_pos();
	/*
	Draw all collision boxes as squares
	*/
	void draw_collision_boxes();
	void draw_kernel_radius();
	void draw_velocities();
	/*
	Check if particle is colliding with a collision box.
	The particle are then moved to the point of collision and
	an impulse force is returned
	*/
	vec2 check_collision(Particle p, double delta_time);
	vec2 calculate_collision_impulse(Particle p, vec2 wall_normal);
	/*
	Signed distance from point p to line with a point on a line and a direction
	*/
	double signed_distance(vec2 p, vec2 plane_point, vec2 plane_dir);

	vector<double> signed_distance_to_walls(vec2 p, Collision_Box c_b);

	void set_kernel_radius(double radius) {
		KERNEL_RADIUS = radius;
	}

	void change_draw_kernel_radius() {
		is_drawing_kernel = !is_drawing_kernel;
	}

	void change_draw_velocities() {
		is_drawing_velocities = !is_drawing_velocities;
	}

	void move_collision_box(int id);
	/*
	The three different smoothing kernels
	*/
	float poly6_kernel(float r, float d);

	vec2 gradient_kernel(float r, DSC2D::vec2 v, float d);

	float laplacian_kernel(float r, float d);

	void add_costum() {
		GAS_CONSTANT += 0.1;
		cout << GAS_CONSTANT << endl;
	}
	void subtract_costum() {
		GAS_CONSTANT -= 0.1;
		cout << GAS_CONSTANT << endl;
	}
private:
	vector<Collision_Box> collision_boxes;
	bool is_drawing_kernel = false;
	bool is_drawing_velocities = false;
	/*
	These are user controlled constant to get the desired fluid
	*/
	float MAX_VELOCITY = 1.0;
	double COEFFICIENT_OF_RESTITUTION = 0.4;
	double COEFFICIENT_OF_FRICTION = 0.03;
	double KERNEL_RADIUS = 35.0;
	double GAS_CONSTANT = 0.85;
	double REST_DENSITY = 0.02;
	double VISCOCITY_TERM = 25;

};





