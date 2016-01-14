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

using namespace DSC2D;
using namespace std;

namespace particlesystem {
	struct Particle
	{
		Particle(vec2 p, int id) :
			id(id),
			pos(p),
			vel(vec2(0.0)),
			mass(1.0),
			density(0.0),
			local_pressure(0.0),
			pressure(vec2(0)),
			viscocity(vec2(0)),
			external_forces(vec2(0)),
			surface_tension(vec2(0))
		{}
		// id to refer to a single particles in the particle system
		int id;
		// Position and velocity of one particle
		vec2 pos;
		vec2 vel;
		// Mass, density, and local pressure
		double mass;
		double density;
		double local_pressure;
		// The different forces applied to a particle
		vec2 pressure;
		vec2 viscocity;
		vec2 external_forces;
		vec2 surface_tension;
	};

	class Particle_System
	{
	public:
		Particle_System(int max);
		/*
		Creates a new particle with a position
		*/
		void create_particle(vec2 p);
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
		given particle. Optionally one can provide a list to get all the distances.
		*/
		vector<Particle> get_particles_in_circle(Particle p, double radius, vector<double> &distances);
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
	private:

		vector<Particle> particles;
		const int MAX_NUMBER_OF_PARTICLES;

	};
}
