#pragma once
#include <vector>
#include "util.h"
#include "DSC.h"


namespace DSC2D{
	class ParticleSystem
	{
	public:
		ParticleSystem();
		// struct of a particle
		struct Particle
		{
			Particle(vec3 p, int i) :
			pos(p),
			velocity(vec3(0)),
			mass(1),
			density(0),
			pressure(vec3(0)),
			localPressure(0),
			viscosity(vec3(0)),
			externalForces(vec3(0)),
			index(i)
			{}
			Particle(vec3 p, int i, vec3 v) :
				pos(p),
				velocity(v),
				mass(0.02),
				density(0),
				pressure(vec3(0)),
				localPressure(0),
				viscosity(vec3(0)),
				externalForces(vec3(0)),
				surface_tension(vec3(0)),
				index(i)
			{}
			vec3 pos;
			vec3 velocity;
			float mass;
			float density;
			vec3 pressure;
			float localPressure;
			vec3 viscosity;
			vec3 externalForces;
			int index;
			vec3 surface_tension;
			bool is_inside;
		};
		void update(DSC2D::DeformableSimplicialComplex& dsc);
		float smoothKernel(float r, float d);
		DSC2D::vec3 gradientKernel(float r, DSC2D::vec3 v, float d);
		float lapKernel(float r, float d);
		void boundingBox(int i);
		void openDoor(){
			doorCounter += 1;
		};
		void startWater(){
			shotWater = !shotWater;
		}
		int get_no_particles(){
			return particles.size();
		}
		vec2 get_particle_pos(int id){
			return vec2(particles[id].pos[0], particles[id].pos[1]);
		}
		float get_particle_mass(int id){
			return particles[id].mass;
		}
		float get_border_floor(){
			return border_floor;
		}
		float get_border_left(){
			return border_left;
		}
		float get_border_right(){
			return border_right;
		}
		bool is_particle_inside(int id) {
			return particles[id].is_inside;
		}
		float fluid_density(vec2 x, float radius);
	private:
		std::vector <Particle> particles;
		std::vector <Particle> border_particles_floor;
		std::vector <Particle> border_particles_right;
		std::vector <Particle> border_particles_left;
		int doorCounter;
		float border_floor;
		float border_right;
		float border_left = 80;
		bool shotWater = false;
	};

}