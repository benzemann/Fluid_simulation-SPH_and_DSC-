#pragma once
#include "util.h"
#include <vector>
#ifdef WIN32
#include <GL/glew.h>
#include <GL/glut.h>

#else
#include <GEL/GL/glew.h>
#include <GLUT/glut.h>
#include <GEL/GLGraphics/SOIL.h>
#endif

//using namespace DSC2D;
using namespace std;

namespace particlesystem {
	struct Particle
	{
		Particle(DSC2D::vec2 p, int id) :
			id(id),
			pos(p),
			pos_0(DSC2D::vec2(0.0)),
			vel(DSC2D::vec2(0.0)),
			vel_0(DSC2D::vec2(0.0)),
			a_0(DSC2D::vec2(0.0)),
			mass(1.0),
			density(0.0),
			local_pressure(0.0),
			pressure(DSC2D::vec2(0)),
			viscocity(DSC2D::vec2(0)),
			external_forces(DSC2D::vec2(0)),
			surface_tension(DSC2D::vec2(0)),
			is_fixed(false)
		{}
		// id to refer to a single particles in the particle system
		int id;
		// Position and velocity of one particle
		DSC2D::vec2 pos;
		DSC2D::vec2 vel;
		DSC2D::vec2 pos_0;
		DSC2D::vec2 vel_0;
		DSC2D::vec2 a_0;
		// Mass, density, and local pressure
		double mass;
		double density;
		double local_pressure;
		// The different forces applied to a particle
		DSC2D::vec2 pressure;
		DSC2D::vec2 viscocity;
		DSC2D::vec2 external_forces;
		DSC2D::vec2 surface_tension;
		bool is_fixed;
	};

	class Particle_System
	{
	public:
		Particle_System(int max);
		/*
		Creates a new particle with a position
		*/
		int create_particle(DSC2D::vec2 p);
		/*
		Returns the number of particles in the system.
		*/
		int get_number_of_particles() {
			return particles.size();
		};
		/*
		Returns the particle with the given id.
		*/
		Particle get_particle(int id) {
			return particles[id];
		};
		/*
		Returns a pointer to the particle with the given id.
		*/
		Particle* get_particle_ptr(int id) {
			return &particles[id];
		};
		/*
		Returns false if no close particles else true, the input closest will be assigned
		the closest particles. max_distance defines the maximum searching distance.
		*/
		bool get_closest_particle(Particle p, Particle &closest, double max_dist = 10.0);
		/*
		Returns a list of all particles in the circle with radius around the
		given particle.
		*/
		vector<Particle> get_particles_in_circle(Particle p, double radius);
		/*
		Renders all particles as small red dots
		*/
		void draw();
		/*
		Render a circle around each particle with radius r
		*/
		void draw_circles(double r);
		/*
		Deletes all particles in the particle system
		*/
		void clear() {
			particles.clear();
		}
		/*
		Update grid
		*/
		void update_grid();

		void create_grid(int width, int height, double cell_size);
	private:

		vector<Particle> particles;
		const int MAX_NUMBER_OF_PARTICLES;

	};
}
