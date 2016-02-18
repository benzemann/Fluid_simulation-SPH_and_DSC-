#include "SPH.h"

Particle_System particle_system(500);

SPH::SPH(int tmp)
{
	
}

void SPH::init() {
	// Create all particles in a square
	int indx = 0;
	/*for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 1; j++) {
			vec2 particle_pos = vec2(75 + (i * 5), 250 + (j * 5));
			particle_system.create_particle(particle_pos);
			indx++;
		}
	}*/
	/*for (int i = 0; i < 100; i++) {
		vec2 particle_pos = vec2(5.0 * i, 200.0);
		particle_system.create_particle(particle_pos);
		Particle* p_ptr = particle_system.get_particle_ptr(indx);
		indx++;
		p_ptr->is_fixed = true;
	}
	for (int i = 0; i < 100; i++) {
		vec2 particle_pos = vec2(100.0, i * 5.0);
		particle_system.create_particle(particle_pos);
		Particle* p_ptr = particle_system.get_particle_ptr(indx);
		indx++;
		p_ptr->is_fixed = true;
	}
	for (int i = 0; i < 100; i++) {
		vec2 particle_pos = vec2(200.0, i * 5.0);
		particle_system.create_particle(particle_pos);
		Particle* p_ptr = particle_system.get_particle_ptr(indx);
		indx++;
		p_ptr->is_fixed = true;
	}*/
	// Create grid (width, height, cell size)

	particle_system.create_grid(34.0, 34.0, 15.0);
	// Create collision boxes, the first is always the movable one!
	if (collision_boxes.size() == 0) {
		vec3 grey = vec3(100, 100, 100);
		vec3 brown = vec3(122, 88, 37);

		create_collision_box(vec2(350.0, 250.0), 200.0, 100.0, brown);
		create_collision_box(vec2(250.0, 50.0), 125.0, 500.0, grey);
		create_collision_box(vec2(250.0, 500.0), 100.0, 500.0, grey);
		create_collision_box(vec2(525.0, 250.0), 500.0, 100.0, grey);
		create_collision_box(vec2(-25.0, 250.0), 500.0, 100.0, grey);
	}
}

void SPH::reset() {
	particle_system.clear();
	collision_boxes.clear();
	init();
}

void SPH::update(double delta_time) {
	//cout << gradient_kernel(10.0, vec2(1.0,0.0), KERNEL_RADIUS) << endl;
	if(enabled_grid)
		particle_system.update_grid();
	// First loop calculate density and local pressure
	for (int i = 0; i < particle_system.get_number_of_particles(); i++) {
		Particle* p_ptr = particle_system.get_particle_ptr(i);
		vector<double> distances;
		vector<Particle> close_particles = get_close_particles(p_ptr);
		p_ptr->density = calculate_density(*p_ptr, close_particles);
		p_ptr->local_pressure = calculate_local_pressure(*p_ptr);
		//cout << p_ptr->density << " " << p_ptr->local_pressure << endl;
	}
	// Second loop to calculate Viscocity, pressure, and external forces
	for (int i = 0; i < particle_system.get_number_of_particles(); i++) {
		Particle* p_ptr = particle_system.get_particle_ptr(i);
		vector<double> distances;
		vector<Particle> close_particles = get_close_particles(p_ptr);
		p_ptr->pressure = calculate_pressure(*p_ptr, close_particles);
		p_ptr->viscocity = calculate_viscocity(*p_ptr, close_particles);
		p_ptr->surface_tension = calculate_surface_tension(*p_ptr, close_particles);
		//cout << p_ptr->viscocity << " " << p_ptr->pressure << endl;
		
		p_ptr->external_forces = calculate_external_forces(*p_ptr);
	}
	update_velocity(delta_time);
	update_position(delta_time);
}

void SPH::update_tmp(double delta_time) {
	//cout << gradient_kernel(10.0, vec2(1.0,0.0), KERNEL_RADIUS) << endl;
	if (enabled_grid)
		particle_system.update_grid();
	// First loop calculate density and local pressure
	for (int i = 0; i < particle_system.get_number_of_particles(); i++) {
		Particle* p_ptr = particle_system.get_particle_ptr(i);
		vector<double> distances;
		vector<Particle> close_particles = get_close_particles(p_ptr);
		p_ptr->density = calculate_density(*p_ptr, close_particles);
		p_ptr->local_pressure = calculate_local_pressure(*p_ptr);
		//cout << p_ptr->density << " " << p_ptr->local_pressure << endl;
	}
	// Second loop to calculate Viscocity, pressure, and external forces
	for (int i = 0; i < particle_system.get_number_of_particles(); i++) {
		Particle* p_ptr = particle_system.get_particle_ptr(i);
		vector<double> distances;
		vector<Particle> close_particles = get_close_particles(p_ptr);
		p_ptr->pressure = calculate_pressure(*p_ptr, close_particles);
		p_ptr->viscocity = calculate_viscocity(*p_ptr, close_particles);
		p_ptr->surface_tension = calculate_surface_tension(*p_ptr, close_particles);
		//cout << p_ptr->viscocity << " " << p_ptr->pressure << endl;

		p_ptr->external_forces = calculate_external_forces(*p_ptr);
	}
}

vector<Particle> SPH::get_close_particles(Particle* p_ptr) {
	
	if (enabled_grid) {
		// Use spatial data grid
		return particle_system.get_particles_in_circle(*p_ptr, KERNEL_RADIUS);
	}
	else {
		// Loop over all particle
		vector<Particle> close_particles;
		for (int i = 0; i < particle_system.get_number_of_particles(); i++) {
			Particle p = particle_system.get_particle(i);
			vec2 p_to_p = p.pos - p_ptr->pos;
			double distance = p_to_p.length();
			if (distance <= KERNEL_RADIUS && p.id != p_ptr->id) {
				close_particles.push_back(p);
			}
		}
		return close_particles;
	}
}

double SPH::calculate_density(Particle p, vector<Particle> close_particles) {
	double density = 0.0;
	for each(Particle close_particle in close_particles) {
		vec2 p_to_p = p.pos - close_particle.pos;
		density += (close_particle.mass * poly6_kernel(p_to_p.length(), KERNEL_RADIUS));
	}
	/*if (!p.is_fixed)
		cout << poly6_kernel(KERNEL_RADIUS - 0.1, KERNEL_RADIUS) << endl;*/
	return density + p.mass * poly6_kernel(0.0, KERNEL_RADIUS);
}

double SPH::calculate_local_pressure(Particle p) {
	return GAS_CONSTANT * (p.density - REST_DENSITY);
}

vec2 SPH::calculate_surface_tension(Particle p, vector<Particle> close_particles) {
	vec2 inward_normal = vec2(0.0);
	double lap_c = 0.0;
	for each(Particle close_particle in close_particles) {
		vec2 p_to_p = p.pos - close_particle.pos;
		inward_normal += (close_particle.mass / close_particle.density) * gradient_kernel(p_to_p.length(), p_to_p, KERNEL_RADIUS);
		lap_c += (close_particle.mass / close_particle.density) * laplacian_kernel(p_to_p.length(), KERNEL_RADIUS);
	}
	vec2 surface_tension = vec2(0.0);
	if (inward_normal.length() > 0.035) {
		inward_normal.normalize();
		surface_tension = -400.0 * inward_normal * lap_c;
	}
	return surface_tension;
}

vec2 SPH::calculate_pressure(Particle p, vector<Particle> close_particles) {

	vec2 pressure = vec2(0.0);
	for each(Particle close_particle in close_particles) {
		vec2 p_to_p = p.pos - close_particle.pos;
		if (close_particle.density != 0.0) {
			pressure += close_particle.mass * ((p.local_pressure + close_particle.local_pressure) / (2.0 * close_particle.density)) * gradient_kernel(p_to_p.length(), p_to_p, KERNEL_RADIUS);		
			//pressure += close_particle.mass * ((p.local_pressure/pow(p.density,2.0)) + (close_particle.local_pressure / pow(close_particle.density, 2.0))) * gradient_kernel(p_to_p.length(), p_to_p, KERNEL_RADIUS);
		}
	}
	return -pressure;
}

vec2 SPH::calculate_viscocity(Particle p, vector<Particle> close_particles) {

	vec2 viscocity = vec2(0.0);
	for each(Particle close_particle in close_particles) {
		vec2 p_to_p = p.pos - close_particle.pos;
		if(close_particle.density != 0.0)
			viscocity += p.mass * ((close_particle.vel - p.vel) / close_particle.density) * laplacian_kernel(p_to_p.length(), KERNEL_RADIUS);
	}
	return VISCOCITY_TERM * viscocity;

}

vec2 SPH::calculate_external_forces(Particle p) {
	vec2 gravity = vec2(0.0, -9.4);
	return gravity;
}

void SPH::update_velocity(double delta_time) {
	for (int i = 0; i < particle_system.get_number_of_particles(); i++) {
		Particle* p_ptr = particle_system.get_particle_ptr(i);
		if (!p_ptr->is_fixed) {
			/*if (p_ptr->pressure.length() > 100.0)
				p_ptr->pressure = vec2(0.0);*/
			
			//vec2 v_minus_half_t = p_ptr->vel - delta_time * 0.5 * p_ptr->a_0;
			vec2 acceleration = p_ptr->pressure + p_ptr->viscocity + p_ptr->external_forces;
			/*if (p_ptr->density != 0.0)
				acceleration = acceleration / p_ptr->density;*/
			//p_ptr->a_0 = acceleration;
			//vec2 v_plus_half_t = v_minus_half_t + delta_time * acceleration;
			
			//cout << p_ptr->density << endl;
			vec2 velocity = p_ptr->vel + acceleration;
			//p_ptr->vel_0 = p_ptr->vel;
			//p_ptr->vel_0 = v_plus_half_t;
			//p_ptr->vel = (v_minus_half_t + v_plus_half_t) / 2.0;
			p_ptr->vel = velocity;
			p_ptr->vel += check_collision(*p_ptr, delta_time);
			/*if (p_ptr->vel.length() > 150.0)
				p_ptr->vel = normalize(p_ptr->vel) * 150.0;*/

			
		}
	}
	/*for (int i = 0; i < particle_system.get_number_of_particles(); i++) {
		Particle* p_ptr = particle_system.get_particle_ptr(i);
		p_ptr->pos_0 = p_ptr->pos;
		if (!p_ptr->is_fixed)
			p_ptr->pos += p_ptr->vel * delta_time;
		
	}
	update_tmp(delta_time);
	for (int i = 0; i < particle_system.get_number_of_particles(); i++) {
		Particle* p_ptr = particle_system.get_particle_ptr(i);
		if (!p_ptr->is_fixed) {
			
			vec2 acceleration = p_ptr->pressure + p_ptr->viscocity + p_ptr->external_forces;
			
			vec2 velocity = p_ptr->vel + acceleration;

			p_ptr->vel = velocity;

			p_ptr->vel += check_collision(*p_ptr, delta_time);
			if (p_ptr->vel.length() > 150.0)
				p_ptr->vel = normalize(p_ptr->vel) * 150.0;
		}
	}
	for (int i = 0; i < particle_system.get_number_of_particles(); i++) {
		Particle* p_ptr = particle_system.get_particle_ptr(i);
		if (!p_ptr->is_fixed)
			p_ptr->pos = p_ptr->pos_0;
	}*/
}

void SPH::update_position(double delta_time) {
	
	for (int i = 0; i < particle_system.get_number_of_particles(); i++) {
		Particle* p_ptr = particle_system.get_particle_ptr(i);
		if(!p_ptr->is_fixed)
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

void SPH::create_particle_at_mouse_pos() {
	POINT curPos;
	BOOL result = GetCursorPos(&curPos);
	static int step = -15;
	if (result)
	{
		RECT rect;
		HWND hwnd = GetForegroundWindow();

		if (GetWindowRect(hwnd, &rect)) {
			int range = 15 - (-15) + 1;
			int num = rand() % range - 15;
			int x = abs(rect.left - curPos.x) + step + (num-10);
			step += 15;
			if (step > 15)
				step = -15;
			int y = rect.bottom - curPos.y + num;
			particle_system.create_particle(vec2(x, y));
		}
	}
}

void SPH::create_collision_box(vec2 pos, float height, float width, vec3 c) {
	Collision_Box new_collision_box(pos, height, width);
	new_collision_box.set_color(c);
	// top layer of ghost particles
	/*for (int i = 0; i < ((width - (KERNEL_RADIUS / 2.0)) / 10.0); i++) {

		vec2 p_pos = vec2(pos[0] - (width / 2.0) + (KERNEL_RADIUS / 4.0) + (i * 10.0) + 15.0,
			pos[1] + (height / 2.0) - 15.0);
		int indx = particle_system.create_particle(p_pos);
		Particle* p_ptr = particle_system.get_particle_ptr(indx);
		p_ptr->is_fixed = true;
		new_collision_box.ghost_particles.push_back(indx);
	}*/

	/*for (int j = 0; j < 1; j++) {
		for (int i = 0; i < ((width - (KERNEL_RADIUS / 2.0)) / 10.0); i++) {

			vec2 p_pos = vec2(pos[0] - (width / 2.0) + (KERNEL_RADIUS / 4.0) + (i * 10.0),
				pos[1] + (height / 2.0));
			int indx = particle_system.create_particle(p_pos);
			Particle* p_ptr = particle_system.get_particle_ptr(indx);
			p_ptr->is_fixed = true;
			new_collision_box.ghost_particles.push_back(indx);
		}
		for (int i = 0; i < ((height - (KERNEL_RADIUS / 2.0)) / 10.0); i++) {
			vec2 p_pos = vec2(pos[0] + (width / 2.0),
				pos[1] - (height / 2.0) + (KERNEL_RADIUS / 4.0) + (i * 10.0));
			int indx = particle_system.create_particle(p_pos);
			Particle* p_ptr = particle_system.get_particle_ptr(indx);
			p_ptr->is_fixed = true;
			new_collision_box.ghost_particles.push_back(indx);
		}
		// left layer of ghost particles
		for (int i = 0; i < ((height - (KERNEL_RADIUS / 2.0)) / 10.0); i++) {
			vec2 p_pos = vec2(pos[0] - (width / 2.0),
				pos[1] - (height / 2.0) + (KERNEL_RADIUS / 4.0) + (i * 10.0));
			int indx = particle_system.create_particle(p_pos);
			Particle* p_ptr = particle_system.get_particle_ptr(indx);
			p_ptr->is_fixed = true;
			new_collision_box.ghost_particles.push_back(indx);
		}
	}*/
	
	// Right layer of ghost particles

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
	Particle* p_ptr = particle_system.get_particle_ptr(p.id);
	/*if (new_pos[0] < 100.0) {
		p_ptr->pos[0] = 101.0;
		p_ptr->pressure = vec2(0.0);
		p_ptr->vel[0] = 0.0;
		p_ptr->vel[1] = 0.0;
	}
	else if (new_pos[0] > 200.0) {
		p_ptr->pos[0] = 199.0;
		p_ptr->pressure = vec2(0.0);
		p_ptr->vel[0] = 0.0;
		p_ptr->vel[1] = 0.0;
	}

	if (new_pos[1] < 200.0) {
		p_ptr->pressure = vec2(0.0);
		p_ptr->pos[1] = 201.0;
		p_ptr->vel[1] = 0.0;
		p_ptr->vel[0] = 0.0;
	}*/

	for each(Collision_Box c_b in collision_boxes) {
	
		if (c_b.is_inside(p.pos)) {
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

			if (c_b.is_inside(p.pos)) {

				vec2 tmp = (wall_normal) * ((c_b.height)*0.5) + c_b.pos;
				Particle* p_ptr = particle_system.get_particle_ptr(p.id);
				//cout << wall_normal << endl;
				if (wall_id == 1 || wall_id == 3) {
					tmp = wall_normal * ((c_b.width)*0.5) + c_b.pos;
					
					p_ptr->pos[0] = tmp[0];
					/*if (p_ptr->pressure[0] > 0 && wall_id == 1) {
						p_ptr->vel = vec2(0.0);
						p_ptr->pressure = vec2(0.0);
					}
					
					if (p_ptr->pressure[0] < 0 && wall_id == 2) {
						p_ptr->vel = vec2(0.0);
						p_ptr->pressure = vec2(0.0);
					}*/
				}
				else {
					//p_ptr->vel[1] = tmp[1] - p_ptr->pos[1];
					//impulse += tmp[1] - p_ptr->pos[1];
					p_ptr->pos[1] = tmp[1];
					/*if (p_ptr->pressure[1] < 0 && wall_id == 0) {
						p_ptr->vel = vec2(0.0);
						p_ptr->pressure = vec2(0.0);
					}
					if (p_ptr->pressure[1] > 0 && wall_id == 3) {
						p_ptr->vel = vec2(0.0);
						p_ptr->pressure = vec2(0.0);
					}*/
				}	
			}

			// teleport to edge of box if 
			if (c_b.is_inside(new_pos)) {
				vec2 tmp = wall_normal * (c_b.height*0.5) + c_b.pos;
				Particle* p_ptr = particle_system.get_particle_ptr(p.id);
				if (wall_id == 1 || wall_id == 3) {
					tmp = wall_normal * (c_b.width*0.5) + c_b.pos;
					//p_ptr->vel[0] = tmp[0] - p_ptr->pos[0];
					//impulse += tmp[0] - p_ptr->pos[0];
					//if (p_ptr->pressure.length() > 30.0) {
						//p_ptr->vel[0] = 0;
					//}
					//p_ptr->pos[0] = tmp[0];
					/*if (p_ptr->pressure[0] > 0 && wall_id == 1)
						p_ptr->vel *= 0.1;
					if (p_ptr->pressure[0] < 0 && wall_id == 2)
						p_ptr->vel *= 0.1;*/
				}
				else {
					//p_ptr->vel[1] = tmp[1] - p_ptr->pos[1];
					//impulse += tmp[1] - p_ptr->pos[1];
					/*if (p_ptr->pressure.length() > 30.0) {
						p_ptr->vel[1] = 1;
					}*/
					//p_ptr->pos[1] = tmp[1];
					/*if (p_ptr->pressure[1] < 0 && wall_id == 0)
						p_ptr->vel *= 0.1;
					if (p_ptr->pressure[1] > 0 && wall_id == 3)
						p_ptr->vel *= 0.1;*/
				}
				impulse += calculate_collision_impulse(p, wall_normal);
			}
			
			if (closest_distance > 5.0) {
				//impulse += calculate_collision_impulse(p, wall_normal);
			}
			else {
				// If particle very close to collider, slow it down

				impulse -= p.vel * 0.1;
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
			Particle p = particle_system.get_particle(i);
			if (!p.is_fixed) {
				glBegin(GL_LINE_LOOP);
				glColor3f(0.0f, 0.0f, 1.0f);
				for (int j = 0; j < 360; j++)
				{
					float degInRad = j*DEG2RAD;
					glVertex2f(p.pos[0] + cos(degInRad)*KERNEL_RADIUS, p.pos[1] + sin(degInRad)*KERNEL_RADIUS);
				}
				glEnd();
			}
		}
	}
}

void SPH::draw_pressure() {
	if (is_drawing_pressure) {
		glBegin(GL_LINES);
		glColor3f(1.0f, 0.0f, 0.0f);
		for (int i = 0; i < particle_system.get_number_of_particles(); i++) {
			Particle p = particle_system.get_particle(i);
			if (!p.is_fixed) {
				glVertex2f(p.pos[0], p.pos[1]);
				glVertex2f(p.pos[0] + p.pressure[0], p.pos[1] + p.pressure[1]);
			}
		}
		glEnd();
	}
}

void SPH::draw_viscocity() {
	if (is_drawing_viscocity) {
		glBegin(GL_LINES);
		
		glColor3f(0.0f, 0.0f, 1.0f);
		for (int i = 0; i < particle_system.get_number_of_particles(); i++) {
			Particle p = particle_system.get_particle(i);
			if (!p.is_fixed) {
				glVertex2f(p.pos[0], p.pos[1]);
				glVertex2f(p.pos[0] + p.viscocity[0], p.pos[1] + p.viscocity[1]);
			}
		}
		glEnd();
	}
}

void SPH::draw_surface_tension() {
	if (is_drawing_surface_tension) {
		glBegin(GL_LINES);

		glColor3f(0.0f, 1.0f, 1.0f);
		for (int i = 0; i < particle_system.get_number_of_particles(); i++) {
			Particle p = particle_system.get_particle(i);
			if (!p.is_fixed) {
				glVertex2f(p.pos[0], p.pos[1]);
				glVertex2f(p.pos[0] + p.surface_tension[0], p.pos[1] + p.surface_tension[1]);
			}
		}
		glEnd();
	}
}

void SPH::draw_velocities() {
	if (is_drawing_velocities) {
		glBegin(GL_LINES);
		glColor3f(0.0f, 1.0f, 0.0f);
		for (int i = 0; i < particle_system.get_number_of_particles(); i++) {
			Particle p = particle_system.get_particle(i);
			if (!p.is_fixed) {
				glVertex2f(p.pos[0], p.pos[1]);
				glVertex2f(p.pos[0] + (p.vel[0] * 0.25), p.pos[1] + (p.vel[1] * 0.25));
			}
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
			vec2 translation = vec2(x, y) - collision_boxes[id].pos;
			collision_boxes[id].pos = vec2(x, y);
			for (int i = 0; i < collision_boxes[id].ghost_particles.size(); i++) {
				Particle* p_ptr = particle_system.get_particle_ptr(collision_boxes[id].ghost_particles[i]);
				p_ptr->pos += translation;
			}
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
		//return -(945 / (32 * M_PI * pow(d, 9))) * v * pow(pow(d, 2) - pow(r, 2), 2);
		return -(45 / (M_PI*pow(d, 6))) * normalize(v) *pow(d - r, 2);
		//return ((15 / (M_PI*pow(d, 6))) * pow(d-r,3)) * normalize(v);
	}
}

float SPH::laplacian_kernel(float r, float d) {
	if (r < 0 || r > d) {
		return 0;
	}
	else {
		//return -(945 / (32 * M_PI * pow(d, 9))) * (pow(d, 2) - pow(r, 2)) * (3 * pow(d, 2) - 7 * pow(r, 2));
		return (45 / (M_PI*pow(d, 6)))*(d - r);
		//return (15 / (2.0*M_PI*pow(d, 3))) - (pow(r, 3) / (2 * pow(d, 3))) + (pow(r, 2) / pow(d, 2)) + (d / (2 * r)) - 1;
	}
}

vector<double*> SPH::get_user_variables_ptr() {
	vector<double*> out_vector;
	out_vector.push_back(&COEFFICIENT_OF_RESTITUTION);
	out_vector.push_back(&COEFFICIENT_OF_FRICTION);
	out_vector.push_back(&KERNEL_RADIUS);
	out_vector.push_back(&GAS_CONSTANT);
	out_vector.push_back(&REST_DENSITY);
	out_vector.push_back(&VISCOCITY_TERM);
	return out_vector;
}