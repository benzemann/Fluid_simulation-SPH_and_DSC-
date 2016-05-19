#include "Particle_System.h"
#include "Grid.h"

using namespace particlesystem;

Particle_System::Particle_System()
{
}

int Particle_System::create_particle(vec3 p) {
	Particle new_particle = Particle(p, particles.size());
	particles.push_back(new_particle);
	//grid.insert_particle(new_particle); // NO GRID IMPLEMENTED YET //
	return new_particle.id;
}


void Particle_System::update_grid() {
	grid->clear_grid();
	for each(Particle p in particles) {
		grid->insert_particle(p);
	}
}

void Particle_System::create_grid(int width, int length, int height, double cell_size) {
	grid = new Grid();
	grid->create_grid(width, length, height, cell_size);
}

bool Particle_System::get_closest_particle(Particle p, Particle &closest, double max_dist) {
	//if (grid.get_closest_particle(p, closest, max_dist)) {
		//return true;
	//}

	return false;
}

vector<Particle> Particle_System::get_particles_in_sphere(Particle p, double radius) {
	return grid->get_particles_in_sphere(p, radius);
};

void Particle_System::draw( vec3 light_pos) {
	
	GLfloat light_position[] = { -light_pos[0] , -light_pos[1], -light_pos[2] , 1.0 };
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	
	for each (Particle p in particles)
	{

		vec3 pos = p.pos_scaled;
		glPushMatrix();
		glTranslated(pos[0], pos[1], pos[2]);
		glutSolidSphere(0.5, 5, 5);
		glPopMatrix();
	}
	

}
