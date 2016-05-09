#include "Grid.h"


Grid::Grid(double cell_size, int width, int height) :
	CELL_SIZE(cell_size),
	GRID_WIDTH(width),
	GRID_HEIGHT(height)
{
	init_grid();
}
Grid::Grid()
{

}

void Grid::init_grid() {
	// Initialize the grid
	for (int i = 0; i < GRID_WIDTH*GRID_HEIGHT; i++) {
		vector<Particle> new_list;
		grid.push_back(new_list);
	}
}

void Grid::create_grid(int width, int height, double cell_size) {
	CELL_SIZE = cell_size;
	GRID_WIDTH = width;
	GRID_HEIGHT = height;
	init_grid();
}

void Grid::insert_particle(Particle p) {
	// Get the grid cell coordinates of the particle and insert it into the grid
	vec3 coords = get_grid_coords(p.pos);
	if (coords[0] + GRID_WIDTH * coords[1] < grid.size() && coords[0] + GRID_WIDTH * coords[1] >= 0)
		grid[coords[0] + GRID_WIDTH * coords[1]].push_back(p);
}

void Grid::clear_grid() {
	// Go over all grid cells and clear all vectors for particles
	for (int i = 0; i < grid.size(); i++) {
		grid[i].clear();
	}
}

vector<Particle> Grid::get_particles_in_sphere(Particle p, double radius) {
	// Get the grid coordinates of the particle
	vec3 p_coords = get_grid_coords(p.pos);
	// Get the number of cells to search depending on the radius of the circle
	int num_of_cells = (int)floor((radius / CELL_SIZE));

	vector<Particle> particles;
	// Loop over grid cells around the particle
	for (int i = -num_of_cells; i <= num_of_cells; i++) {
		for (int j = -num_of_cells; j <= num_of_cells; j++) {
			float x = i + p_coords[0];
			float y = j + p_coords[1];
			if (x >= 0 && x < GRID_WIDTH && y >= 0 && y < GRID_HEIGHT) {
				// Loop over all particles in current grid cell
				for each(Particle par in grid[x + GRID_WIDTH * y]) {
					// Calculate distance to current particle and check if its inside the circle
					vec3 v = par.pos - p.pos;
					double distance = v.length();
					if (distance <= radius && par.id != p.id) {
						// Store the particle and distance
						particles.push_back(par);
					}
				}
			}
		}
	}
	return particles;
}

bool Grid::get_closest_particle(Particle p, Particle &closest, double max_dist) {
	// Current_radius is set to two cell sizes
	double current_radius = CELL_SIZE*2.0;
	vector<Particle> close_particles;
	// Search for close particles by expanding the radius for every iteration
	while (current_radius < max_dist) {
		close_particles = get_particles_in_sphere(p, current_radius);
		if (close_particles.size() != 0) {
			break;
		}
		current_radius += CELL_SIZE;
	}
	// No close particles
	if (close_particles.size() == 0)
		return false;
	// Calculate the closest of all found particles
	Particle closest_particle = close_particles[0];
	double dis = max_dist;
	for each(Particle par in close_particles) {
		vec3 v = par.pos - p.pos;
		if (v.length() < dis && p.is_inside) {
			closest_particle = par;
			dis = v.length();
		}
	}
	// Set the input closest to the closest particle and return true
	closest = closest_particle;
	return true;
}

vec3 Grid::get_grid_coords(vec3 pos) {
	// Calculate the grid cell coordinates of the input position
	int x = (int)floor(pos[0] / CELL_SIZE);
	int y = (int)floor(pos[1] / CELL_SIZE);

	return vec3(x, y, 0); // NEEDS Z
}