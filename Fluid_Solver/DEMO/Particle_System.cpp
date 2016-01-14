#include "Particle_System.h"
#include "Grid.h"

using namespace particlesystem;

Grid grid(10.0,50.0,50.0);

Particle_System::Particle_System(int max):
	MAX_NUMBER_OF_PARTICLES(max)
{
	
}

void Particle_System::create_particle(vec2 p) {
	Particle new_particle = Particle(p, particles.size());
	particles.push_back(new_particle);
	grid.insert_particle(new_particle);
}

void Particle_System::update_grid() {
	grid.clear_grid();
	for each(Particle p in particles) {
		grid.insert_particle(p);
	}
}

bool Particle_System::get_closest_particle(Particle p, Particle &closest, double max_dist) {
	if (grid.get_closest_particle(p, closest, max_dist)) {
		return true;
	}
	return false;
}

vector<Particle> Particle_System::get_particles_in_circle(Particle p, double radius, vector<double> &distances) {
	return grid.get_particles_in_circle(p, radius, distances);
};

void Particle_System::draw() {
	
	glPointSize(5.0);
	glBegin(GL_POINTS);
	glColor3d(1.0, 0.0, 0.0);
	for each(Particle p in particles) {
		glVertex3d(p.pos[0], p.pos[1], 0.0);
	}
	glEnd();
}

void Particle_System::draw_circles(double r) {
	const float DEG2RAD = 3.14159 / 180;
	glBegin(GL_LINE_LOOP);

	for (int i = 0; i < 360; i++)
	{
		float degInRad = i*DEG2RAD;
		glVertex2f(cos(degInRad)*r, sin(degInRad)*r);
	}

	glEnd();
}