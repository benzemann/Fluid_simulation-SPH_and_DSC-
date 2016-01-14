#include "SPH.h"

Particle_System particle_system(500);

SPH::SPH(int tmp)
{
	
}

void SPH::init() {
	// Create all particles in a square
	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 10; j++) {
			vec2 particle_pos = vec2(75 + (i * 15), 250 + (j * 15));
			particle_system.create_particle(particle_pos);
		}
	}
	// Create collision boxes, the first is always the movable one!
	vec3 grey = vec3(100, 100, 100);
	vec3 brown = vec3(122, 88, 37);
	create_collision_box(vec2(350.0, 250.0), 100.0, 100.0, brown);
	create_collision_box(vec2(250.0, 50.0), 125.0, 500.0, grey);
	create_collision_box(vec2(250.0, 450.0), 100.0, 500.0, grey);
	create_collision_box(vec2(500.0, 250.0), 500.0, 50.0, grey);
	create_collision_box(vec2(0.0, 250.0), 500.0, 50.0, grey);
}

void SPH::reset() {
	particle_system.clear();
	collision_boxes.clear();
	init();
}

void SPH::update(double delta_time) {
	particle_system.update_grid();
	// First loop calculate density and local pressure
	for (int i = 0; i < particle_system.get_number_of_particles(); i++) {
		Particle* p_ptr = particle_system.get_particle_ptr(i);
		vector<double> distances;
		//vector<Particle> close_particles = particle_system.get_particles_in_circle(*p_ptr, KERNEL_RADIUS, distances);
		//p_ptr->density = calculate_density(*p_ptr, close_particles);
		p_ptr->local_pressure = calculate_local_pressure(*p_ptr);
	}
	// Second loop to calculate Viscocity, pressure, and external forces
	for (int i = 0; i < particle_system.get_number_of_particles(); i++) {
		Particle* p_ptr = particle_system.get_particle_ptr(i);
		vector<double> distances;
		//vector<Particle> close_particles = particle_system.get_particles_in_circle(*p_ptr, KERNEL_RADIUS, distances);
		//p_ptr->pressure = calculate_pressure(*p_ptr, close_particles);
		//p_ptr->viscocity = calculate_viscocity(*p_ptr, close_particles);
		p_ptr->external_forces = calculate_external_forces(*p_ptr);
	}
	update_velocity(delta_time);
	update_position(delta_time);
}

double SPH::calculate_density(Particle p, vector<Particle> close_particles) {
	double density = 0.0;
	for each(Particle close_particle in close_particles) {
		vec2 p_to_p = p.pos - close_particle.pos;
		density += p.mass * poly6_kernel(p_to_p.length(), KERNEL_RADIUS);
	}
	return density;
}

double SPH::calculate_local_pressure(Particle p) {
	return GAS_CONSTANT * (p.density - REST_DENSITY);
}

vec2 SPH::calculate_pressure(Particle p, vector<Particle> close_particles) {

	vec2 pressure = vec2(0.0);
	for each(Particle close_particle in close_particles) {
		vec2 p_to_p = p.pos - close_particle.pos;
		if (close_particle.density != 0.0) {
			pressure += p.mass * ((p.local_pressure + close_particle.local_pressure) / close_particle.density) * gradient_kernel(p_to_p.length(), p_to_p, KERNEL_RADIUS);
		}
	}
	return -pressure;
}

vec2 SPH::calculate_viscocity(Particle p, vector<Particle> close_particles) {

	vec2 viscocity = vec2(0.0);
	for each(Particle close_particle in close_particles) {
		vec2 p_to_p = p.pos - close_particle.pos;
		viscocity += p.mass * ((close_particle.vel - p.vel) / close_particle.local_pressure) * laplacian_kernel(p_to_p.length(), KERNEL_RADIUS);
	}
	return VISCOCITY_TERM * viscocity;

}

vec2 SPH::calculate_external_forces(Particle p) {
	vec2 gravity = vec2(0.0, -0.4);
	return gravity;
}

void SPH::update_velocity(double delta_time) {
	for (int i = 0; i < particle_system.get_number_of_particles(); i++) {
		Particle* p_ptr = particle_system.get_particle_ptr(i);
		vec2 acceleration = p_ptr->pressure + p_ptr->viscocity + p_ptr->external_forces;
		vec2 velocity = p_ptr->vel + acceleration;
		p_ptr->vel = velocity;
		p_ptr->vel += check_collision(*p_ptr, delta_time);
	}
}

void SPH::update_position(double delta_time) {
	
	for (int i = 0; i < particle_system.get_number_of_particles(); i++) {
		Particle* p_ptr = particle_system.get_particle_ptr(i);
		p_ptr->pos += p_ptr->vel * delta_time;	
	}
}

void SPH::draw_collision_boxes() {
	glBegin(GL_QUADS);
	for each (Collision_Box c_b in collision_boxes)
	{
		vec2 pos = c_b.pos;
		float half_height = c_b.height * 0.5;
		float half_width = c_b.width * 0.5;
		vec3 color = c_b.color;
		glColor3f(color[0]/255.0, color[1]/255.0, color[2]/255.0);
		glVertex2f(pos[0] - half_width, pos[1] - half_height);
		glVertex2f(pos[0] + half_width, pos[1] - half_height);
		glVertex2f(pos[0] + half_width, pos[1] + half_height);
		glVertex2f(pos[0] - half_width, pos[1] + half_height);
		
	}
	glEnd();
}

void SPH::draw_particles() {
	particle_system.draw();
}

void SPH::create_collision_box(vec2 pos, float height, float width, vec3 c) {
	Collision_Box new_collision_box(pos, height, width);
	new_collision_box.set_color(c);
	collision_boxes.push_back(new_collision_box);
}

void SPH::create_collision_box_at_mouse_pos() {
	POINT curPos;
	BOOL result = GetCursorPos(&curPos);
	if (result)
	{
		RECT rect;
		HWND hwnd = GetForegroundWindow();

		if (GetWindowRect(hwnd, &rect)) {
			int x = abs(rect.left - curPos.x);
			int y = rect.bottom - curPos.y;
			vec3 grey = vec3(100, 100, 100);
			create_collision_box(vec2(x, y), 100.0, 100.0, grey);
		}
	}
}

vec2 SPH::check_collision(Particle p, double delta_time) {
	vec2 new_pos = p.pos + (p.vel * delta_time);
	vec2 impulse = vec2(0.0);
	for each(Collision_Box c_b in collision_boxes) {
	
		if (c_b.is_inside(new_pos)) {
			// particle is inside collision box
			vec2 wall_normal = vec2(0.0, 0.0);
			// Check distance to closest collision wall
			vector<double> distances = signed_distance_to_walls(new_pos, c_b);
			double closest_distance = distances[0];
			int wall_id = 0;
			int i = 0;
			for each(double distance in distances) {
				if (distance < closest_distance) {
					closest_distance = distance;
					wall_id = i;
				}
				i++;
			}
			// Determine which wall normal to use
			switch (wall_id)
			{
			case 0:
				wall_normal = vec2(0.0, 1.0);
				break;
			case 1:
				wall_normal = vec2(1.0, 0.0);
				break;
			case 2:
				wall_normal = vec2(0.0, -1.0);
				break;
			case 3:
				wall_normal = vec2(-1.0, 0.0);
			default:
				break;
			}
			// teleport to edge of box if 
			if (c_b.is_inside(p.pos)) {
				vec2 tmp = wall_normal * (c_b.height*0.5) + c_b.pos;
				Particle* p_ptr = particle_system.get_particle_ptr(p.id);
				if (wall_id == 1 || wall_id == 3) {
					tmp = wall_normal * (c_b.width*0.5) + c_b.pos;
					p_ptr->pos[0] = tmp[0];
				}
				else {
					p_ptr->pos[1] = tmp[1];
				}
				//impulse += calculate_collision_impulse(p, wall_normal);
			}
			if (closest_distance > 0.01 || wall_id != 0) {
				impulse += calculate_collision_impulse(p, wall_normal);
			}
			else {
				// If particle very close to collider, slow it down
				impulse -= p.vel * 0.05;
			}
		}
	}
	return impulse;
}

vec2 SPH::calculate_collision_impulse(Particle p, vec2 wall_normal) {

	double normal_speed = dot(wall_normal, p.vel);
	double j = -(1.0 + COEFFICIENT_OF_RESTITUTION) * p.mass * normal_speed;
	vec2 tmp = (-p.vel - (normal_speed * wall_normal));
	if (tmp.length() > 0.001) {
		tmp = normalize(tmp);
	}
	vec2 J = j * wall_normal - COEFFICIENT_OF_FRICTION * j * tmp;

	return J;
}

vector<double> SPH::signed_distance_to_walls(vec2 p, Collision_Box c_b) {
	
	vec2 pos = c_b.pos;
	float height = c_b.height;
	float width = c_b.width;

	vector<double> distances;
	// upper wall
	vec2 wall_point(pos[0] - (width*0.5), pos[1] + (height*0.5));
	vec2 wall_dir(1.0, 0.0);
	double distance = signed_distance(p, wall_point, wall_dir);
	distances.push_back(distance);
	// Rightmost wall
	wall_point = vec2(pos[0] + (width*0.5), pos[1] - (height*0.5));
	wall_dir = vec2(0.0, 1.0);
	distance = signed_distance(p, wall_point, wall_dir);
	distances.push_back(distance);
	// Lower wall
	wall_point = vec2(pos[0] - (width*0.5), pos[1] - (height*0.5));
	wall_dir = vec2(1.0, 0.0);
	distance = signed_distance(p, wall_point, wall_dir);
	distances.push_back(distance);
	// Leftmost wall
	wall_point = vec2(pos[0] - (width*0.5), pos[1] - (height*0.5));
	wall_dir = vec2(0.0, 1.0);
	distance = signed_distance(p, wall_point, wall_dir);
	distances.push_back(distance);

	return distances;
}

double SPH::signed_distance(vec2 p, vec2 plane_point, vec2 plane_dir) {

	vec2 w = p - plane_point;
	vec2 u = normalize(plane_dir); 
	vec2 d = w - (dot(w, u) * u);
	return d.length();
}

void SPH::draw_kernel_radius() {
	if (is_drawing_kernel) {
		
		const float DEG2RAD = 3.14159 / 180;
		for (int i = 0; i < particle_system.get_number_of_particles(); i++) {
			glBegin(GL_LINE_LOOP);
			glColor3f(0.0f, 0.0f, 1.0f);
			Particle p = particle_system.get_particle(i);
			for (int j = 0; j < 360; j++)
			{
				float degInRad = j*DEG2RAD;
				glVertex2f(p.pos[0] + cos(degInRad)*KERNEL_RADIUS, p.pos[1] + sin(degInRad)*KERNEL_RADIUS);
			}
			glEnd();
		}



	}
}

void SPH::draw_velocities() {
	if (is_drawing_velocities) {
		glBegin(GL_LINES);
		glColor3f(0.0f, 1.0f, 0.0f);
		for (int i = 0; i < particle_system.get_number_of_particles(); i++) {
			Particle p = particle_system.get_particle(i);

			glVertex2f(p.pos[0], p.pos[1]);
			glVertex2f(p.pos[0] + p.vel[0], p.pos[1] + p.vel[1]);

		}
		glEnd();
	}
}

void SPH::move_collision_box(int id) {
	POINT curPos;
	BOOL result = GetCursorPos(&curPos);
	if (result)
	{
		RECT rect;
		HWND hwnd = GetForegroundWindow();

		if (GetWindowRect(hwnd, &rect)) {
			int x = abs(rect.left - curPos.x);
			int y = rect.bottom - curPos.y;
			collision_boxes[id].pos = vec2(x, y);
		}
	}
}

float SPH::poly6_kernel(float r, float d) {
	if (r < 0 || r > d) {
		return 0;
	}
	else {
		return (315 / (64 * M_PI*pow(d, 9)))*pow(pow(d, 2) - pow(r, 2), 3);
	}
}

vec2 SPH::gradient_kernel(float r, DSC2D::vec2 v, float d) {
	if (r < 0 || r > d) {
		return vec2(0);
	}
	else {
		return (-45 / (M_PI*pow(d, 6)))*pow(d - r, 2) * normalize(v);
	}
}

float SPH::laplacian_kernel(float r, float d) {
	if (r < 0 || r > d) {
		return 0;
	}
	else {
		return (45 / (M_PI*pow(d, 6)))*(d - r);
	}
}