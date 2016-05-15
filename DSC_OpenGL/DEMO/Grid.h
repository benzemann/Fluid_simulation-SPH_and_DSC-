#pragma once
#include "Particle_System.h"
#include <vector>
#include <CGLA/Vec3d.h>

typedef CGLA::Vec3d vec3;
using namespace std;
using namespace particlesystem;

class Grid
{
public:
	/*
	Creates a grid at the starting point and out of the x-axis.
	cell_size defines the size of each cell and width and height defines the number
	of cells horisontal and vertical.
	*/
	Grid(double cell_size, int width, int height);
	Grid();
	/*
	Inserts a particle into the grid.
	*/
	void insert_particle(Particle p);
	/*
	Clear the grid.
	*/
	void clear_grid();
	/*
	Returns false if no close particles else true, the input closest will be assigned
	the closest particles. max_distance defines the maximum searching distance.
	*/
	bool get_closest_particle(Particle p, Particle &closest, double max_dist = 10.0);
	/*
	Returns a list of pointers to all particles in the circle with radius around the
	given particle.
	*/
	vector<Particle> get_particles_in_sphere(Particle p, double radius);

	void init_grid();

	void create_grid(int width, int length, int height, double cell_size);

private:
	/*
	Returns a vector with the grid coordinates cooresponding to the input position.
	*/
	vec3 get_grid_coords(vec3 pos);

private:

	double CELL_SIZE;
	int GRID_WIDTH;
	int GRID_LENGTH;
	int GRID_HEIGHT;

	vector<Particle>*** grid;
};
