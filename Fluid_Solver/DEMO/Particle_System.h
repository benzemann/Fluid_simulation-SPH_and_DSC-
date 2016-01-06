#pragma once
#include "util.h"
#include <vector>

using namespace DSC2D;
using namespace std;
namespace particlesystem {
	struct Particle
	{
		Particle(vec2 p, int id) :
			id(id),
			pos(p),
			vel(vec2(0)),
			mass(0.0),
			density(0.0),
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
		// Mass and density
		double mass;
		double density;
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
		void create_particle(vec2 p) {
			Particle new_particle = Particle(p, particles.size());
			particles.push_back(new_particle);
		}	
		/*
			Returns the number of particles in the system.
		*/
		int get_number_of_particles() {
			return number_of_particles;
		};
		/* 
			Returns the particle with the given id.
		*/
		Particle get_particle(int id) {
			return particles[id];
		};
		/*
			Returns false if no close particles else true, the input closest will be assigned
			the closest particles. max_distance defines the maximum searching distance.
		*/
		bool get_closest_particle(Particle p, Particle &closest, double max_dist = 10.0) {
			if (grid.get_closest_particle(p, closest, max_dist)) {
				return true;
			}
			return false;
		}
		/*
			Returns a list of pointers to all particles in the circle with radius around the
			given particle. Optionally one can provide a list to get all the distances.
		*/
		vector<Particle> get_particles_in_circle(Particle p, double radius, vector<double> &distances) {
			return grid.get_particles_in_circle(p, radius, distances);
		};
	private:

		vector<Particle> particles;
		int number_of_particles;
		const int MAX_NUMBER_OF_PARTICLES;
		Grid grid;
	};
}


