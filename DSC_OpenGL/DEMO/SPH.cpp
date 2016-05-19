#include "SPH.h"

Particle_System ps;

SPH::SPH()
{
}

void SPH::init() {
	ps = Particle_System();
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 5; j++) {
			for (int k = 0; k < 45; k++) {
				int id = ps.create_particle(vec3( 
					100 + (i * 15.0),
					100 + (j * 15.0),
					150 + (k * 15.0)
					));
			}
		}
	}

	create_collision_box(vec3(200.0, 200.0, 0.0), 100.0, 1500, 1500);

	create_collision_box(vec3(200.0, 0.0, 200.0), 10000.0, 1500, 100);
	create_collision_box(vec3(200.0, 200.0, 200.0), 10000.0, 1500, 100);

	create_collision_box(vec3(400.0, 200.0, 200.0), 10000.0, 100, 1500);
	create_collision_box(vec3(0.0, 200.0, 200.0), 10000.0, 100, 1500);

	ps.create_grid(80, 80, 80, 15);



	int iso_dims = 54;
	double iso_cell_size = 7.5;

	iso_surface_nodes = new Node***[iso_dims];
	for (int i = 0; i < iso_dims; i++) {
		iso_surface_nodes[i] = new Node**[iso_dims];
		for (int j = 0; j < iso_dims; j++) {
			iso_surface_nodes[i][j] = new Node*[iso_dims];
		}
	}

	for (int x = 0; x < iso_dims; x++) {
		for (int y = 0; y < iso_dims; y++) {
			for (int z = 0; z < iso_dims; z++) {
				iso_surface_nodes[x][y][z] = new Node(vec3(iso_cell_size*x, iso_cell_size*y, iso_cell_size*z), 0.0);
			}
		}
	}
}

void SPH::reset() {
	ps.clear();
	collision_boxes.clear();
	init();
}
int tmp = 0;
void SPH::update(double delta_time) {
	/*if (tmp < 2000) {
		int rx = rand() % 20;
		int ry = rand() % 20;
		if (tmp % 3 == 0)
			ps.create_particle(vec3(200.0 + rx, 200.0 + ry, 400.0));
	}
	tmp++;*/

	ps.update_grid();
	// First loop calculate density and local pressure
	avg_density = 0.0;
	for (int i = 0; i < ps.get_number_of_particles(); i++) {
		Particle* p_ptr = ps.get_particle_ptr(i);
		vector<double> distances;
		vector<Particle> close_particles = get_close_particles(p_ptr);
		p_ptr->old_density = p_ptr->density;
		p_ptr->density = calculate_density(*p_ptr, close_particles);
		avg_density += p_ptr->density;
		p_ptr->local_pressure = calculate_local_pressure(*p_ptr);
	}
	// Second loop to calculate Viscocity, pressure, and external forces
	for (int i = 0; i < ps.get_number_of_particles(); i++) {

		Particle* p_ptr = ps.get_particle_ptr(i);
		vector<double> distances;
		vector<Particle> close_particles = get_close_particles(p_ptr);
		p_ptr->a = calculate_a(*p_ptr, close_particles);
		p_ptr->pressure = calculate_pressure(*p_ptr, close_particles);
		p_ptr->viscocity = calculate_viscocity(*p_ptr, close_particles);
		p_ptr->external_forces = calculate_external_forces(*p_ptr);
	}

}

double SPH::CFL_correction(double delta_time) {
	double d_t = delta_time;
	double v_max = 0.0;
	for (int i = 0; i < ps.get_number_of_particles(); i++) {
		if (ps.get_particle(i).vel.length() > v_max) {
			v_max = ps.get_particle(i).vel.length();
		}
	}
	double CFL = 0.2 * ((KERNEL_RADIUS*2.0) / v_max);
	if (d_t > CFL) {
		//cout << delta;
		d_t = CFL;
		if (delta > 0.008)
			delta = 0.008;
		//cout << " corrected to " << delta << endl;
	}
	return d_t;
}

double SPH::calculate_density(Particle p, vector<Particle> close_particles) {
	double density = 0.0;
	for each(Particle close_particle in close_particles) {
		vec3 p_to_p = p.pos - close_particle.pos;
		density += (close_particle.mass * poly6_kernel(p_to_p.length(), KERNEL_RADIUS));
	}
	return density + p.mass * poly6_kernel(0.0, KERNEL_RADIUS);
}

double SPH::calculate_local_pressure(Particle p) {
	return GAS_CONSTANT * (p.density - REST_DENSITY);
}

vec3 SPH::calculate_pressure(Particle p, vector<Particle> close_particles) {

	vec3 pressure = vec3(0.0);
	for each(Particle close_particle in close_particles) {
		vec3 p_to_p = p.pos - close_particle.pos;
		if (close_particle.density != 0.0) {
			pressure += close_particle.mass * ((p.local_pressure + close_particle.local_pressure) / (2.0 * close_particle.density)) * gradient_kernel(p_to_p.length(), p_to_p, KERNEL_RADIUS);
		}
	}
	return -pressure;
}

vec3 SPH::calculate_viscocity(Particle p, vector<Particle> close_particles) {

	vec3 viscocity = vec3(0.0);
	for each(Particle close_particle in close_particles) {
		vec3 p_to_p = p.pos - close_particle.pos;
		if (close_particle.density != 0.0)
			viscocity += p.mass * ((close_particle.vel - p.vel) / close_particle.density) * laplacian_kernel(p_to_p.length(), KERNEL_RADIUS);
	}
	return VISCOCITY_TERM * viscocity;

}

vec3 SPH::calculate_external_forces(Particle p) {
	vec3 gravity = vec3(0.0, 0.0, -9.4);
	return gravity;
}

void SPH::update_velocity(double delta_time) {

	for (int i = 0; i < ps.get_number_of_particles(); i++) {

		Particle* p_ptr = ps.get_particle_ptr(i);

		vec3 acceleration = p_ptr->pressure + p_ptr->viscocity + p_ptr->external_forces;

		vec3 velocity = p_ptr->vel + acceleration;

		p_ptr->vel = velocity;

		p_ptr->vel += check_collision(*p_ptr, delta_time);

	}

}

void SPH::update_position(double delta_time) {

	for (int i = 0; i < ps.get_number_of_particles(); i++) {
		Particle* p_ptr = ps.get_particle_ptr(i);
		p_ptr->pos += p_ptr->vel * delta_time;
		p_ptr->pos_scaled = (p_ptr->pos / down_scale);
		p_ptr->pos_scaled = vec3(p_ptr->pos_scaled[0] - 1.0, p_ptr->pos_scaled[1] - 1.0, p_ptr->pos_scaled[2] - 1.0);
	}

}

void SPH::correct_density_error(double delta_time) {

	double density_error = avg_density - REST_DENSITY;
	int max_iterations = 2;
	int iterations = 0;
	while (density_error > 0.00006 && iterations <= max_iterations) {
		double sum_density = 0.0;
		for (int i = 0; i < ps.get_number_of_particles(); i++) {
			Particle* p_ptr = ps.get_particle_ptr(i);
			vec3 next_pos = p_ptr->pos * (p_ptr->vel * delta_time);
			vector<Particle> close_particles = get_close_particles_to_pos(next_pos, KERNEL_RADIUS);
			double d = calculate_density(*p_ptr, close_particles);
			sum_density += d;
			p_ptr->density = d;
		}

		for (int i = 0; i < ps.get_number_of_particles(); i++) {
			Particle* p_ptr = ps.get_particle_ptr(i);
			double k_i = (p_ptr->density - REST_DENSITY) * (1 / (delta_time*delta_time)) * p_ptr->a;
			vector<Particle> close_particles = get_close_particles(p_ptr, KERNEL_RADIUS);
			vec3 tmp = vec3(0.0);
			for each(Particle close_particle in close_particles) {
				vec3 p_to_p = p_ptr->pos - close_particle.pos;
				double k_j = (close_particle.density - REST_DENSITY) * (1 / delta_time * delta_time) * close_particle.a;
				if (p_ptr->density != 0.0 && close_particle.density != 0.0)
					tmp += close_particle.mass * ((k_i / p_ptr->density) + (k_j / close_particle.density)) * gradient_kernel(p_to_p.length(), p_to_p, KERNEL_RADIUS);
			}
			p_ptr->vel = p_ptr->vel - delta_time * tmp;
		}

		avg_density = sum_density / ps.get_number_of_particles();
		density_error = avg_density - REST_DENSITY;
		iterations++;
	}
}

double SPH::calculate_a(Particle p, vector<Particle> close_particles) {
	vec3 v = vec3(0.0);
	double l = 0.0;
	for each(Particle close_particle in close_particles) {
		vec3 p_to_p = p.pos - close_particle.pos;
		vec3 tmp = close_particle.mass * gradient_kernel(p_to_p.length(), p_to_p, KERNEL_RADIUS);
		v += tmp;
		l += tmp.length() * tmp.length();
	}
	double a = ((v.length()*v.length()) + l);
	if (a < 0.00001) {
		a = 0.00001;
	}
	double a_m = p.density / a;
	//cout << "a: " << a << " a_m: " << a_m << endl;
	return a_m;
}

vector<Particle> SPH::get_close_particles(Particle* p_ptr, double radius) {

	if (radius == -1)
		radius = KERNEL_RADIUS;

	return ps.get_particles_in_sphere(*p_ptr, radius);

	vector<Particle> close_particles;
	for (int i = 0; i < ps.get_number_of_particles(); i++) {
		Particle p = ps.get_particle(i);	
		vec3 p_to_p = p.pos - p_ptr->pos;
		double distance = p_to_p.length();
		if (distance <= radius && p.id != p_ptr->id) {
			close_particles.push_back(p);
		}
	}
	return close_particles;
}

vector<Particle> SPH::get_close_particles_to_pos(vec3 pos, double radius) {
	Particle tmp_particle = Particle(pos, -1);
	vector<Particle> close_particles = get_close_particles(&tmp_particle, radius);
	return close_particles;
}

Particle SPH::get_closest_particle(vec3 pos) {
	Particle tmp_particle = Particle(pos, -1);
	Particle closest_particle = Particle(pos, -1);
	ps.get_closest_particle(tmp_particle, closest_particle, 1000.0);
	return closest_particle;
}

void SPH::create_collision_box(vec3 pos, float height, float width, float length) {

	Collision_Box new_collision_box(pos, height, width, length);
	collision_boxes.push_back(new_collision_box);

}

vec3 SPH::check_collision(Particle p, double delta_time) {

	vec3 new_pos = p.pos + (p.vel * delta_time);
	vec3 impulse = vec3(0.0);
	Particle* p_ptr = ps.get_particle_ptr(p.id);

	for each(Collision_Box c_b in collision_boxes) {
		if (c_b.is_inside(p.pos)) {

			vec3 wall_normal = vec3(0.0);

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

			switch (wall_id)
			{
			case 0: // z
				wall_normal = vec3(0.0, 0.0, 1.0);
				break;
			case 1: // -z
				wall_normal = vec3(0.0, 0.0, -1.0);
				break;
			case 2: // x
				wall_normal = vec3(1.0, 0.0, 0.0);
				break;
			case 3: // -x
				wall_normal = vec3(-1.0, 0.0, 0.0);
				break;
			case 4: // y
				wall_normal = vec3(0.0, 1.0, 0.0);
				break;
			case 5: // -y
				wall_normal = vec3(0.0, -1.0, 0.0);
				break;
			default:
				break;
			}

			if (c_b.is_inside(new_pos)) {
				impulse += calculate_collision_impulse(p, wall_normal);
			}
			if (closest_distance <= 5.0) {
				// If particle very close to collider, slow it down
				//impulse -= p.vel * 0.1;
			}
		}
		
	}

	return impulse;

}

vec3 SPH::calculate_collision_impulse(Particle p, vec3 wall_normal) {


	double normal_speed = dot(wall_normal, p.vel);
	double j = -(1.0 + COEFFICIENT_OF_RESTITUTION) * p.mass * normal_speed;
	vec3 tmp = (-p.vel - (normal_speed * wall_normal));
	if (tmp.length() > 0.001) {
		tmp = normalize(tmp);
	}
	vec3 J = j * wall_normal - COEFFICIENT_OF_FRICTION * j * tmp;

	return J;

}

vector<double> SPH::signed_distance_to_walls(vec3 p, Collision_Box c_b) {

	vec3 pos = c_b.pos;
	float height = c_b.height;
	float width = c_b.width;
	float lenght = c_b.length;

	vector<double> distances;

	// z wall 
	vec3 wall_point = vec3(pos[0], pos[1] - (lenght*0.5), pos[2] + (height*0.5));
	vec3 wall_dir = vec3(0.0, 0.0, 1.0);
	double distance = signed_distance(p, wall_point, wall_dir);
	distances.push_back(distance);
	// -z wall
	wall_point = vec3(pos[0], pos[1] - (lenght*0.5), pos[2] - (height*0.5));
	wall_dir = vec3(0.0, 0.0, -1.0);
	distance = signed_distance(p, wall_point, wall_dir);
	distances.push_back(distance);
	// x wall
	wall_point = vec3(pos[0] + (width*0.5), pos[1] - (lenght*0.5), pos[2]);
	wall_dir = vec3(1.0, 0.0, 0.0);
	distance = signed_distance(p, wall_point, wall_dir);
	distances.push_back(distance);
	// -x wall
	wall_point = vec3(pos[0] - (width*0.5), pos[1] - (lenght*0.5), pos[2]);
	wall_dir = vec3(-1.0, 0.0, 0.0);
	distance = signed_distance(p, wall_point, wall_dir);
	distances.push_back(distance);
	// y wall
	wall_point = vec3(pos[0], pos[1] + (lenght*0.5), pos[2] - (height*0.5));
	wall_dir = vec3(0.0, 1.0, 0.0);
	distance = signed_distance(p, wall_point, wall_dir);
	distances.push_back(distance);
	// -y wall
	wall_point = vec3(pos[0], pos[1] - (lenght*0.5), pos[2] - (height*0.5));
	wall_dir = vec3(0.0, -1.0, 0.0);
	distance = signed_distance(p, wall_point, wall_dir);
	distances.push_back(distance);

	return distances;
}

double SPH::signed_distance(vec3 p, vec3 plane_point, vec3 plane_dir) {

	float    sb, sn, sd;
	sn = -dot(plane_dir, (p - plane_point));
	sd = dot(plane_dir, plane_dir);
	sb = sn / sd;
	vec3 base_point = p + sb * plane_dir;
	vec3 base_to_p = p - base_point;
	return base_to_p.length();
}

void SPH::draw_particles(vec3 light_pos) {
	ps.draw(light_pos);
}

void SPH::draw_collision_boxes(vec3 light_pos) {

	GLfloat light_position[] = { -light_pos[0] , -light_pos[1], -light_pos[2] , 1.0 };
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);

	for (int i = 0; i < collision_boxes.size(); i++) {
		
		//glLoadIdentity();                 // Reset the model-view matrix
		vec3 pos_scaled = collision_boxes[i].pos / down_scale;
		pos_scaled = vec3(pos_scaled[0] - 1, pos_scaled[1] - 1, pos_scaled[2] - 1);
		vec3 scale_scaled = vec3(collision_boxes[i].width / (down_scale*2.0), collision_boxes[i].length / (down_scale*2), collision_boxes[i].height / (down_scale*2));
		glPushMatrix();
		glTranslatef(pos_scaled[0], pos_scaled[1], pos_scaled[2]);
		glScalef(scale_scaled[0], scale_scaled[1], scale_scaled[2]);

		glBegin(GL_QUADS);                // Begin drawing the color cube with 6 quads
		glNormal3f(0.0f, 1.0f, 0.0f);
		glVertex3f(1.0f, 1.0f, -1.0f);
		glVertex3f(-1.0f, 1.0f, -1.0f);
		glVertex3f(-1.0f, 1.0f, 1.0f);
		glVertex3f(1.0f, 1.0f, 1.0f);

		glNormal3f(0.0f, -1.0f, 0.0f);
		glVertex3f(1.0f, -1.0f, 1.0f);
		glVertex3f(-1.0f, -1.0f, 1.0f);
		glVertex3f(-1.0f, -1.0f, -1.0f);
		glVertex3f(1.0f, -1.0f, -1.0f);

		glNormal3f(0.0f, 0.0f, 1.0f);
		glVertex3f(1.0f, 1.0f, 1.0f);
		glVertex3f(-1.0f, 1.0f, 1.0f);
		glVertex3f(-1.0f, -1.0f, 1.0f);
		glVertex3f(1.0f, -1.0f, 1.0f);

		glNormal3f(0.0f, 0.0f, -1.0f);
		glVertex3f(1.0f, -1.0f, -1.0f);
		glVertex3f(-1.0f, -1.0f, -1.0f);
		glVertex3f(-1.0f, 1.0f, -1.0f);
		glVertex3f(1.0f, 1.0f, -1.0f);

		glNormal3f(-1.0f, 0.0f, 0.0f);
		glVertex3f(-1.0f, 1.0f, 1.0f);
		glVertex3f(-1.0f, 1.0f, -1.0f);
		glVertex3f(-1.0f, -1.0f, -1.0f);
		glVertex3f(-1.0f, -1.0f, 1.0f);

		glNormal3f(1.0f, 0.0f, 0.0f);
		glVertex3f(1.0f, 1.0f, -1.0f);
		glVertex3f(1.0f, 1.0f, 1.0f);
		glVertex3f(1.0f, -1.0f, 1.0f);
		glVertex3f(1.0f, -1.0f, -1.0f);

		glEnd();  // End of drawing color-cube
		glPopMatrix();
	}

	

}

void SPH::update_iso_surface() {
	
	iso_points.clear();
	


	for (int x = 0; x < 53; x++) {
		for (int y = 0; y < 53; y++) {
			for (int z = 0; z < 53; z++) {

				vector<Particle> close_particles = get_close_particles_to_pos(iso_surface_nodes[x][y][z]->pos);


				double sum = 0.0;
				for each (Particle p in close_particles) {
					vec3 p_to_p = iso_surface_nodes[x][y][z]->pos - p.pos;
					sum += poly6_kernel(p_to_p.length(), 30);
				}
				iso_surface_nodes[x][y][z]->value = sum;
				
			}
		}
	}

	for (int x = 0; x < 53; x++) {
		for (int y = 0; y < 53; y++) {
			for (int z = 0; z < 53; z++) {


				double iso_value = 0.0001;


				Node* n000 = iso_surface_nodes[x][y][z];
				Node* n100 = iso_surface_nodes[x + 1][y][z];
				Node* n010 = iso_surface_nodes[x][y + 1][z];
				Node* n110 = iso_surface_nodes[x + 1][y + 1][z];
				Node* n001 = iso_surface_nodes[x][y][z+1];
				Node* n101 = iso_surface_nodes[x + 1][y][z+1];
				Node* n011 = iso_surface_nodes[x][y + 1][z+1];
				Node* n111 = iso_surface_nodes[x + 1][y + 1][z+1];

				if (n000->value > 0.0 || n100->value > 0.0) {
					if ((n000->value <= iso_value && n100->value > iso_value) ||
						(n000->value > iso_value && n100->value <= iso_value)) {
						//cout << n00->value << endl;
						vec3 interpolated_pos = interpolate_iso_nodes(n000, n100, iso_value);
						interpolated_pos = interpolated_pos / down_scale;
						interpolated_pos = vec3(-1 + interpolated_pos[0], -1 + interpolated_pos[1], -1 + interpolated_pos[2]);
						iso_points.push_back(interpolated_pos);
					}
				}
				if (n100->value > 0.0 || n110->value > 0.0) {
					if ((n100->value <= iso_value && n110->value > iso_value) ||
						(n100->value > iso_value && n110->value <= iso_value)) {

						vec3 interpolated_pos = interpolate_iso_nodes(n100, n110, iso_value);

						interpolated_pos = interpolated_pos / down_scale;
						interpolated_pos = vec3(-1 + interpolated_pos[0], -1 + interpolated_pos[1], -1 + interpolated_pos[2]);
						iso_points.push_back(interpolated_pos);
					}
				}
				if (n000->value > 0.0 || n010->value > 0.0) {
					if ((n000->value <= iso_value && n010->value > iso_value) ||
						(n000->value > iso_value && n010->value <= iso_value)) {

						vec3 interpolated_pos = interpolate_iso_nodes(n000, n010, iso_value);

						interpolated_pos = interpolated_pos / down_scale;
						interpolated_pos = vec3(-1 + interpolated_pos[0], -1 + interpolated_pos[1], -1 + interpolated_pos[2]);
						iso_points.push_back(interpolated_pos);
					}
				}
				if (n010->value > 0.0 || n110->value > 0.0) {
					if ((n010->value <= iso_value && n110->value > iso_value) ||
						(n010->value > iso_value && n110->value <= iso_value)) {

						vec3 interpolated_pos = interpolate_iso_nodes(n010, n110, iso_value);

						interpolated_pos = interpolated_pos / down_scale;
						interpolated_pos = vec3(-1 + interpolated_pos[0], -1 + interpolated_pos[1], -1 + interpolated_pos[2]);
						iso_points.push_back(interpolated_pos);
					}
				}

				if (n001->value > 0.0 || n101->value > 0.0) {
					if ((n001->value <= iso_value && n101->value > iso_value) ||
						(n001->value > iso_value && n101->value <= iso_value)) {
						//cout << n00->value << endl;
						vec3 interpolated_pos = interpolate_iso_nodes(n001, n101, iso_value);
						interpolated_pos = interpolated_pos / down_scale;
						interpolated_pos = vec3(-1 + interpolated_pos[0], -1 + interpolated_pos[1], -1 + interpolated_pos[2]);
						iso_points.push_back(interpolated_pos);
					}
				}
				if (n101->value > 0.0 || n111->value > 0.0) {
					if ((n101->value <= iso_value && n111->value > iso_value) ||
						(n101->value > iso_value && n111->value <= iso_value)) {

						vec3 interpolated_pos = interpolate_iso_nodes(n101, n111, iso_value);

						interpolated_pos = interpolated_pos / down_scale;
						interpolated_pos = vec3(-1 + interpolated_pos[0], -1 + interpolated_pos[1], -1 + interpolated_pos[2]);
						iso_points.push_back(interpolated_pos);
					}
				}
				if (n001->value > 0.0 || n011->value > 0.0) {
					if ((n001->value <= iso_value && n011->value > iso_value) ||
						(n001->value > iso_value && n011->value <= iso_value)) {

						vec3 interpolated_pos = interpolate_iso_nodes(n001, n011, iso_value);

						interpolated_pos = interpolated_pos / down_scale;
						interpolated_pos = vec3(-1 + interpolated_pos[0], -1 + interpolated_pos[1], -1 + interpolated_pos[2]);
						iso_points.push_back(interpolated_pos);
					}
				}
				if (n011->value > 0.0 || n111->value > 0.0) {
					if ((n011->value <= iso_value && n111->value > iso_value) ||
						(n011->value > iso_value && n111->value <= iso_value)) {

						vec3 interpolated_pos = interpolate_iso_nodes(n011, n111, iso_value);

						interpolated_pos = interpolated_pos / down_scale;
						interpolated_pos = vec3(-1 + interpolated_pos[0], -1 + interpolated_pos[1], -1 + interpolated_pos[2]);
						iso_points.push_back(interpolated_pos);
					}
				}
				if (n000->value > 0.0 || n001->value > 0.0) {
					if ((n000->value <= iso_value && n001->value > iso_value) ||
						(n000->value > iso_value && n001->value <= iso_value)) {

						vec3 interpolated_pos = interpolate_iso_nodes(n000, n001, iso_value);

						interpolated_pos = interpolated_pos / down_scale;
						interpolated_pos = vec3(-1 + interpolated_pos[0], -1 + interpolated_pos[1], -1 + interpolated_pos[2]);
						iso_points.push_back(interpolated_pos);
					}
				}
				if (n100->value > 0.0 || n101->value > 0.0) {
					if ((n100->value <= iso_value && n101->value > iso_value) ||
						(n100->value > iso_value && n101->value <= iso_value)) {

						vec3 interpolated_pos = interpolate_iso_nodes(n100, n101, iso_value);

						interpolated_pos = interpolated_pos / down_scale;
						interpolated_pos = vec3(-1 + interpolated_pos[0], -1 + interpolated_pos[1], -1 + interpolated_pos[2]);
						iso_points.push_back(interpolated_pos);
					}
				}
				if (n010->value > 0.0 || n011->value > 0.0) {
					if ((n010->value <= iso_value && n011->value > iso_value) ||
						(n010->value > iso_value && n011->value <= iso_value)) {

						vec3 interpolated_pos = interpolate_iso_nodes(n010, n011, iso_value);

						interpolated_pos = interpolated_pos / down_scale;
						interpolated_pos = vec3(-1 + interpolated_pos[0], -1 + interpolated_pos[1], -1 + interpolated_pos[2]);
						iso_points.push_back(interpolated_pos);
					}
				}
				if (n110->value > 0.0 || n111->value > 0.0) {
					if ((n110->value <= iso_value && n111->value > iso_value) ||
						(n110->value > iso_value && n111->value <= iso_value)) {

						vec3 interpolated_pos = interpolate_iso_nodes(n110, n111, iso_value);

						interpolated_pos = interpolated_pos / down_scale;
						interpolated_pos = vec3(-1 + interpolated_pos[0], -1 + interpolated_pos[1], -1 + interpolated_pos[2]);
						iso_points.push_back(interpolated_pos);
					}
				}
				// DRAW ISO NODES //
				/*if (iso_surface_nodes[x][y][z]->value > 0) {
					vec3 pos = iso_surface_nodes[x][y][z]->pos / down_scale;
					pos = vec3(-1 + pos[0], -1 + pos[1], -1 + pos[2]);
					glPushMatrix();
					glTranslated(pos[0], pos[1], pos[2]);
					glutSolidSphere(0.01, 5, 5);
					glPopMatrix();
				}*/
				
			}
		}
	}
	

	/*for each (vec3 iso_point in iso_points)
	{
		vec3 pos = iso_point;
		glPushMatrix();
		glTranslated(pos[0], pos[1], pos[2]);
		glutSolidSphere(.08, 5, 5);
		glPopMatrix();
	}*/
	
}


vec3 SPH::interpolate_iso_nodes(Node* n_a, Node* n_b, double iso_value) {

	double l_between_nodes = 15.0;
	double t = l_between_nodes * ((iso_value - n_b->value) / (n_a->value - n_b->value));

	vec3 n_a_to_n_b = n_b->pos - n_a->pos;
	n_a_to_n_b = CGLA::normalize(n_a_to_n_b);
	vec3 interpolated_pos = n_a->pos + (n_a_to_n_b * (l_between_nodes - t));
	return interpolated_pos;
}

double SPH::poly6_kernel(float r, float d) {
	if (r < 0 || r > d) {
		return 0;
	}
	else {
		return (315 / (64 * M_PI*pow(d, 9)))*pow(pow(d, 2) - pow(r, 2), 3);
	}
}

vec3 SPH::gradient_kernel(float r, vec3 v, float d) {
	if (r < 0 || r > d) {
		return vec3(0.0);
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