#include "Particle_System.h"
#include "Grid.h"

using namespace particlesystem;

Particle_System::Particle_System(int max):
	MAX_NUMBER_OF_PARTICLES(max),
	grid(Grid(1.0,10,10))
{
	
}
