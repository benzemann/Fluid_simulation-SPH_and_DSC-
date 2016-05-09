#pragma once
#include <vector>
#include <CGLA/Vec3d.h>
#include <GL/glew.h>
#include <GL/glut.h>

typedef CGLA::Vec3d vec3;
using namespace std;

namespace particlesystem {

	struct Particle
	{
		Particle() :
			id(-1),
			pos(vec3(0.0)),
			pos_0(vec3(0.0)),
			vel(vec3(0.0)),
			vel_0(vec3(0.0)),
			a_0(vec3(0.0)),
			mass(1.0),
			density(0.0),
			local_pressure(0.0),
			pressure(vec3(0)),
			viscocity(vec3(0)),
			external_forces(vec3(0)),
			surface_tension(vec3(0)),
			is_fixed(false),
			is_inside(true),
			a(0.0),
			old_density(0.0),
			density_divergence(vec3(0))
		{}
		Particle(vec3 p, int id) :
			id(id),
			pos(p),
			pos_0(vec3(0.0)),
			vel(vec3(0.0)),
			vel_0(vec3(0.0)),
			a_0(vec3(0.0)),
			mass(1.0),
			density(0.0),
			local_pressure(0.0),
			pressure(vec3(0)),
			viscocity(vec3(0)),
			external_forces(vec3(0)),
			surface_tension(vec3(0)),
			is_fixed(false),
			is_inside(true),
			a(0.0),
			old_density(0.0),
			density_divergence(vec3(0))
		{}
		// id to refer to a single particles in the particle system
		int id;
		// Position and velocity of one particle
		vec3 pos;
		vec3 vel;
		vec3 pos_0;
		vec3 vel_0;
		vec3 a_0;
		// Mass, density, and local pressure
		double mass;
		double density;
		double local_pressure;
		// The different forces applied to a particle
		vec3 pressure;
		vec3 viscocity;
		vec3 external_forces;
		vec3 surface_tension;
		bool is_fixed;
		bool is_inside;
		double a;
		double old_density;
		vec3 density_divergence;
	};

	class Particle_System
	{
	public:
		Particle_System();
		/*
		Creates a new particle with a position
		*/
		int create_particle(vec3 p);
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
		void draw(vec3 light_pos);
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
	};
}


