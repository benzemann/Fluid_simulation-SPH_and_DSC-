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
	grid = new vector<Particle>**[GRID_WIDTH];
	for (int i = 0; i < GRID_WIDTH; ++i) {
		grid[i] = new vector<Particle>*[GRID_LENGTH];

		for (int j = 0; j < GRID_LENGTH; j++) {
			grid[i][j] = new vector<Particle>[GRID_HEIGHT];
		}

	}
	/*
	// Initialize the grid
	for (int i = 0; i < GRID_WIDTH*GRID_WIDTH*GRID_HEIGHT; i++) {
		vector<Particle> new_list;
		vector<vector<Particle>> new_vector_list;
		new_vector_list.push_back(new_list);
		grid.push_back(new_vector_list);
	}
	*/
}

void Grid::create_grid(int width,int length, int height, double cell_size) {
	CELL_SIZE = cell_size;
	GRID_WIDTH = width;
	GRID_LENGTH = length;
	GRID_HEIGHT = height;
	init_grid();
}

void Grid::insert_particle(Particle p) {
	vec3 coords = get_grid_coords(p.pos);
	if (coords[0] < GRID_WIDTH && coords[0] >= 0 &&
		coords[1] < GRID_LENGTH && coords[1] >= 0 &&
		coords[2] < GRID_HEIGHT && coords[2] >= 0) 
	{

		grid[(int)coords[0]][(int)coords[1]][(int)coords[2]].push_back(p);

	}
}

void Grid::clear_grid() {
	for (int x = 0; x < GRID_WIDTH; x++) {
		for (int y = 0; y < GRID_LENGTH; y++) {
			for (int z = 0; z < GRID_HEIGHT; z++) {

				grid[x][y][z].clear();

			}
		}
	}
}

vector<Particle> Grid::get_particles_in_sphere(Particle p, double radius) {
	// Get the grid coordinates of the particle
	vec3 p_coords = get_grid_coords(p.pos);

	int num_of_cells = (int)floor((radius / CELL_SIZE));

	vector<Particle> particles;

	for (int i = -num_of_cells; i <= num_of_cells; i++) {
		for (int j = -num_of_cells; j <= num_of_cells; j++) {
			for (int w = -num_of_cells; w <= num_of_cells; w++) {

				int x = i + p_coords[0];
				int y = j + p_coords[1];
				int z = w + p_coords[2];

				if (x >= 0 && x < GRID_WIDTH &&
					y >= 0 && y < GRID_LENGTH &&
					z >= 0 && z < GRID_HEIGHT)
				{
					for each(Particle par in grid[x][y][z]) {
						vec3 v = par.pos - p.pos;
						double distance = v.length();
						if (distance <= radius && par.id != p.id) {
							particles.push_back(par);
						}
					}
				}

			}
		}
	}
	return particles;
}

bool Grid::get_closest_particle(Particle p, Particle &closest, double max_dist) {
	
	/*
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
	*/
	return false;
}

vec3 Grid::get_grid_coords(vec3 pos) {
	int x = (int)floor(pos[0] / CELL_SIZE);
	int y = (int)floor(pos[1] / CELL_SIZE);
	int z = (int)floor(pos[2] / CELL_SIZE);

	return vec3(x, y, z);
}