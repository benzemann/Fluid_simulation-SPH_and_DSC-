#include "ParticleSystem.h"


DSC2D::ParticleSystem::ParticleSystem()
{
	// Start by adding some water particles in a grid inside the container
	int index = 0;
	for (int i = 0; i < 15; i++){
		for (int j = 0; j < 20; j++){
			particles.push_back(Particle(vec3(80 + (i*15), 80 + (j*15), 0.0), index));
			index++;
		}
	}

	doorCounter = 0;
}

// function to calculate the new velocity and position based on pressure, viscocity, and external forces
void DSC2D::ParticleSystem::update(DSC2D::DeformableSimplicialComplex& dsc){
	
	// User-defined variables
	float k = 350; // gass constant
	float envPressure = 0.00008; // environmental pressure
	float viscosityOfFluid = 25; // viscosity of the water fluid
	float kernelRadius = 50.0; // radius of the smoothing kernels
	float tension = 80; // surface tension coefficient
	float tension_threshold = 0.02; // The threshold for surface tension

	border_floor = 80;
	doorCounter++;
	if (doorCounter < 100){ 
		border_right = 240;
	}
	else { // Container expansion
		border_right = 420;
	}
	// calculate densities and local pressure
	for (int i = 0; i < particles.size(); i++){
		particles[i].is_inside = false;
		for (auto fi = dsc.faces_begin(); fi != dsc.faces_end(); ++fi) {
			if (dsc.get_label(*fi) == 1) {
				auto vpos = dsc.get_pos(*fi);
				DSC2D::vec2 ppos = get_particle_pos(i);
				float area = 0.5f * (-vpos[1][1] * vpos[2][0] + vpos[0][1] * (-vpos[1][0] + vpos[2][0]) + vpos[0][0] * (vpos[1][1] - vpos[2][1]) + vpos[1][0] * vpos[2][1]);
				float s = 1.0 / (2.0 * area) * (vpos[0][1] * vpos[2][0] - vpos[0][0] * vpos[2][1] + (vpos[2][1] - vpos[0][1]) * ppos[0] + (vpos[0][0] - vpos[2][0]) * ppos[1]);
				float t = 1.0 / (2.0 * area) * (vpos[0][0] * vpos[1][1] - vpos[0][1] * vpos[1][0] + (vpos[0][1] - vpos[1][1]) * ppos[0] + (vpos[1][0] - vpos[0][0]) * ppos[1]);
				if (s >= 0.0f && t >= 0.0f && 1.0 - s - t > 0.0f) {
					particles[i].is_inside = true;
				}
			}
			
		}
		particles[i].density = 0;
		float x_temp = particles[i].pos[0];
		float y_temp = particles[i].pos[1];
		vec3 inward_normal = vec3(0);
		float lap_c = 0.0f;
		for (int j = 0; j < particles.size(); j++){
			float diff_x = abs(x_temp - particles[j].pos[0]);
			float diff_y = abs(y_temp - particles[j].pos[1]);
			if (j != particles[i].index && diff_x < kernelRadius && diff_y < kernelRadius){ // Check if particle is the same particle and check if the particle is too far away
				vec3 r = particles[i].pos - particles[j].pos;
				particles[i].density += particles[j].mass * smoothKernel(length(r), kernelRadius);
				inward_normal += (particles[j].mass / particles[j].density) * gradientKernel(length(r), r, kernelRadius);
				lap_c += (particles[j].mass / particles[j].density) * lapKernel(length(r), kernelRadius);
			}
		}
		particles[i].localPressure = k * (particles[i].density - envPressure);
		if (particles[i].localPressure < 0.0f){
			particles[i].localPressure = 0.0f;
		}
		float inward_lenght = CGLA::length(inward_normal);
		particles[i].surface_tension = vec3(0);
		if (inward_lenght > tension_threshold){
			particles[i].surface_tension = -tension * CGLA::normalize(inward_normal) * lap_c;
		}
		//printf("inward normal: [%f , %f]", particles[i].pressure[0], particles[i].surface_tension[1]);
	}
	border_particles_floor.clear();
	border_particles_right.clear();
	border_particles_left.clear();
	for (int i = 0; i < particles.size(); i++){
		int indx = 0;
		if (particles[i].pos[1] < border_floor + 5){
			Particle p = Particle(vec3(particles[i].pos[0], (border_floor - particles[i].pos[1]) + border_floor - 5, 0), indx);
			indx++;
			p.velocity = particles[i].velocity;
			p.density = particles[i].density;
			p.localPressure = particles[i].localPressure;
			border_particles_floor.push_back(p);
		}
	}
	// Calculate the forces on the particles (external, pressure, viscosity)
	for (int i = 0; i < particles.size(); i++){
		particles[i].externalForces = vec3(0.0, -0.479, 0.0); // Gravity
		particles[i].pressure = vec3(0);
		particles[i].viscosity = vec3(0);
		float x_temp = particles[i].pos[0];
		float y_temp = particles[i].pos[1];
		for (int j = 0; j < particles.size(); j++){
			float diff_x = abs(x_temp - particles[j].pos[0]);
			float diff_y = abs(y_temp - particles[j].pos[1]);
			if (j != particles[i].index && diff_x < kernelRadius && diff_y < kernelRadius){ // Check if particle is the same particle and check if the particle is too far away
				if (particles[j].density != 0 && particles[i].density != 0){
					vec3 r = particles[i].pos - particles[j].pos;
					particles[i].pressure += (particles[j].mass * ((particles[i].localPressure / pow(particles[i].density, 2.0)) + (particles[j].localPressure / pow(particles[j].density, 2.0))) * gradientKernel(length(r), r, kernelRadius));
					particles[i].viscosity += particles[j].mass * ((particles[j].velocity - particles[i].velocity) / particles[j].density) * lapKernel(length(r), kernelRadius);
				}
			}
		}
		
		for (int j = 0; j < border_particles_floor.size(); j++){ // Use the "ghost particles"
			float diff_x = abs(x_temp - border_particles_floor[j].pos[0]);
			float diff_y = abs(y_temp - border_particles_floor[j].pos[1]);
			if (diff_x < kernelRadius && diff_y < kernelRadius){ // Check if particle is the same particle and check if the particle is too far away
				if (border_particles_floor[j].density != 0 && particles[i].density != 0){
					vec3 r = particles[i].pos - border_particles_floor[j].pos;
					particles[i].pressure += (border_particles_floor[j].mass * ((particles[i].localPressure / pow(particles[i].density, 2.0)) + (border_particles_floor[j].localPressure / pow(border_particles_floor[j].density, 2.0))) * gradientKernel(length(r), r, kernelRadius));
					particles[i].viscosity += border_particles_floor[j].mass * ((border_particles_floor[j].velocity - particles[i].velocity) / border_particles_floor[j].density) * lapKernel(length(r), kernelRadius);
				}
			}
		}
		particles[i].pressure = -particles[i].pressure;
		particles[i].viscosity = particles[i].viscosity * viscosityOfFluid;
	}

	// Update position and velocity of all particles
	for (int i = 0; i < particles.size(); i++){
		boundingBox(i); // Add external forces to ensure particles are inside the container. 
		vec3 acc = particles[i].pressure + particles[i].externalForces + particles[i].viscosity + particles[i].surface_tension;
		if (length(acc) > 10.0){ // Restrict acceleration to 10.0
			acc = normalize(acc) * 10.0;
		}
		particles[i].velocity += acc;
		if (length(particles[i].velocity) > 10.0){ // Restrict velocity to 10.0
			particles[i].velocity = normalize(particles[i].velocity) * 10.0;
		}
		particles[i].pos = particles[i].pos + particles[i].velocity;
	}
	
}
// Calculate fluid density at the point x with a given radius
float DSC2D::ParticleSystem::fluid_density(vec2 x, float radius){
	float fluid_density = 0.0;
	for (int i = 0; i < get_no_particles(); i++){
		vec2 r = x - vec2(particles[i].pos[0], particles[i].pos[1]);
		if (length(r) > 0.0){
			fluid_density += particles[i].mass * smoothKernel(length(r) / radius, radius);
		}
	}
	return fluid_density;
}
// The three different smoothing kernels: 
float DSC2D::ParticleSystem::smoothKernel(float r, float d){
	if (r < 0 || r > d){
		return 0;
	}
	else {
		return (315 / (64 * M_PI*pow(d, 9)))*pow(pow(d, 2) - pow(r, 2), 3);
	}
}

DSC2D::vec3 DSC2D::ParticleSystem::gradientKernel(float r, DSC2D::vec3 v, float d){
	if (r < 0 || r > d){
		return vec3(0);
	}
	else {
		return (-45/ (M_PI*pow(d, 6)))*pow(d - r, 2) * normalize(v);
	}
}

float DSC2D::ParticleSystem::lapKernel(float r, float d){
	if (r < 0 || r > d){
		return 0;
	}
	else {
		return (45 / (M_PI*pow(d, 6)))*(d-r);
	}
}
// Function to apply external forces to the particles if they are outside the container. 
void DSC2D::ParticleSystem::boundingBox(int i){
	if (particles[i].pos[1] < border_floor){
		particles[i].pos[1] = border_floor;
		float bouncyFactor = 0.2;
		float friction = 0.2;
		vec3 initVel = particles[i].velocity;
		vec3 normal = vec3(0, 1, 0);
		float normalSpeed = dot(normal, initVel);
		float j = -(1 + bouncyFactor)*particles[i].mass*normalSpeed;
		vec3 tmp = vec3(initVel[0] - (normalSpeed*normal[0]), initVel[1] - (normalSpeed*normal[1]), 0);
		vec3 J;
		vec3 dotNV = cross(normal, normalize(initVel));
		if (dotNV[0] == 0 && dotNV[1] == 0 && dotNV[2] == 0){
			J = j*normal;
		}
		else{
			J = j*normal - friction * j * normalize(tmp);
		}

		J[2] = 0;
		particles[i].externalForces += J;
	}
	if (particles[i].pos[0] < border_left){
		particles[i].pos[0] = border_left;
		float bouncyFactor = 0.2;
		float friction = 0.2;
		vec3 initVel = particles[i].velocity;
		vec3 normal = vec3(1, 0, 0);
		float normalSpeed = dot(normal, initVel);
		float j = -(1 + bouncyFactor)*particles[i].mass*normalSpeed;
		vec3 tmp = vec3(initVel[0] - (normalSpeed*normal[0]), initVel[1] - (normalSpeed*normal[1]), 0);
		vec3 J;
		vec3 dotNV = cross(normal, normalize(initVel));
		if (dotNV[0] == 0 && dotNV[1] == 0 && dotNV[2] == 0){
			J = j*normal;
		}
		else{
			J = j*normal - friction * j * normalize(tmp);
		}
		J[2] = 0;
		particles[i].externalForces += J;
	}
	if (particles[i].pos[0] > border_right){ 
		particles[i].pos[0] = border_right;
		float bouncyFactor = 0.2;
		float friction = 0.2;
		vec3 initVel = particles[i].velocity;
		vec3 normal = vec3(-1, 0, 0);
		float normalSpeed = dot(normal, initVel);
		float j = -(1 + bouncyFactor)*particles[i].mass*normalSpeed;
		vec3 tmp = vec3(initVel[0] - (normalSpeed*normal[0]), initVel[1] - (normalSpeed*normal[1]), 0);
		vec3 J;
		vec3 dotNV = cross(normal, normalize(initVel));
		if (dotNV[0] == 0 && dotNV[1] == 0 && dotNV[2] == 0){
			J = j*normal;
		}
		else{
			J = j*normal - friction * j * normalize(tmp);
		}
		J[2] = 0;
		particles[i].externalForces += J;
	}
}


