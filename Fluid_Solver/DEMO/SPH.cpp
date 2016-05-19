#include "SPH.h"

Particle_System particle_system(500);

SPH::SPH(int tmp)
{
	
}

void SPH::init() {
	// Create all particles in a square
	int indx = 0;

	/*for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 20; j++) {
			DSC2D::vec2 particle_pos = DSC2D::vec2(500 + (i * 15), 250 + (j * 10));
			particle_system.create_particle(particle_pos);
			indx++;
		}
	}*/
	// create grid for of isosurface nodes
	for (int i = 0; i < 67; i++) {
		vector<Node*> tmp;
		for (int j = 0; j < 67; j++) {
			tmp.push_back(new Node(DSC2D::vec2(i*15.0, j*15.0), 0.0));
		}
		isosurface_nodes.push_back(tmp);
	}

	particle_system.create_grid(67.0, 67.0, 15.0);
	// Create collision boxes, the first is always the movable one!
	if (collision_boxes.size() == 0) {
		DSC2D::vec3 grey = DSC2D::vec3(100, 100, 100);
		DSC2D::vec3 brown = DSC2D::vec3(122, 88, 37);
		create_collision_box(DSC2D::vec2(350.0, 500.0), 1000.0, 100.0, brown);
		create_collision_box(DSC2D::vec2(500.0, 50.0), 125.0, 1000.0, grey);
		create_collision_box(DSC2D::vec2(500.0, 1000.0), 100.0, 1000.0, grey);
		create_collision_box(DSC2D::vec2(1000.0, 500.0), 1000.0, 150.0, grey);
		create_collision_box(DSC2D::vec2(0.0, 500.0), 1000.0, 150.0, grey);
	}
}

void SPH::reset() {
	particle_system.clear();
	collision_boxes.clear();
	init();
}

void SPH::update(double delta_time) {
	if(enabled_grid)
		particle_system.update_grid();
	avg_density = 0.0;
	// First loop calculate density and local pressure
	for (int i = 0; i < particle_system.get_number_of_particles(); i++) {
		Particle* p_ptr = particle_system.get_particle_ptr(i);
		vector<double> distances;
		vector<Particle> close_particles = get_close_particles(p_ptr);
		p_ptr->old_density = p_ptr->density;
		p_ptr->density = calculate_density(*p_ptr, close_particles);
		avg_density += p_ptr->density;
		p_ptr->local_pressure = calculate_local_pressure(*p_ptr);
	}
	// Second loop to calculate Viscocity, pressure, and external forces
	for (int i = 0; i < particle_system.get_number_of_particles(); i++) {
		Particle* p_ptr = particle_system.get_particle_ptr(i);
		vector<double> distances;
		vector<Particle> close_particles = get_close_particles(p_ptr);
		p_ptr->a = calculate_a(*p_ptr, close_particles);
		p_ptr->pressure = calculate_pressure(*p_ptr, close_particles);
		p_ptr->viscocity = calculate_viscocity(*p_ptr, close_particles);
		//p_ptr->surface_tension = calculate_surface_tension(*p_ptr, close_particles);
		p_ptr->surface_tension = calculate_inward_normal(p_ptr->pos) * 2000.0;
		p_ptr->external_forces = calculate_external_forces(*p_ptr);
	}
	//cout << (avg_density / get_no_of_particle()) - REST_DENSITY << endl;
}

vector<Particle> SPH::get_close_particles(Particle* p_ptr, double radius) {
	if (radius == -1)
		radius = KERNEL_RADIUS;
	if (enabled_grid) {
		// Use spatial data grid
		return particle_system.get_particles_in_circle(*p_ptr, radius);
	}
	else {
		// Loop over all particle
		vector<Particle> close_particles;
		for (int i = 0; i < particle_system.get_number_of_particles(); i++) {
			Particle p = particle_system.get_particle(i);
			DSC2D::vec2 p_to_p = p.pos - p_ptr->pos;
			double distance = p_to_p.length();
			if (distance <= radius && p.id != p_ptr->id) {
				close_particles.push_back(p);
			}
		}
		return close_particles;
	}
}

void SPH::CFL_delta_time_update() {
	double v_max = 0.0;
	for (int i = 0; i < particle_system.get_number_of_particles(); i++) {
		if (CGLA::length(particle_system.get_particle(i).vel) > v_max) {
			v_max = CGLA::length(particle_system.get_particle(i).vel);
		}
	}

	double CFL = 0.05 * ((KERNEL_RADIUS*2.0) / v_max);
	if (delta > CFL) {
		//cout << delta;
		delta = CFL;
		if (delta > 0.009)
			delta = 0.009;
		//cout << " corrected to " << delta << endl;
	}
}

void SPH::update_isosurface() {
	if (get_no_of_particle() > 0) {
		for (int i = 0; i < isosurface_nodes.size(); i++) {
			for (int j = 0; j < isosurface_nodes.size(); j++) {

				vector<Particle> close_particles = get_close_particles_to_pos(isosurface_nodes[i][j]->pos);
				double sum = 0.0;
				for each (Particle p in close_particles)
				{


					//if (p.vel.length() < 1300.0) {
					DSC2D::vec2 p_to_p = isosurface_nodes[i][j]->pos - p.pos;
					sum += poly6_kernel(p_to_p.length(), 15);
					//}
					
				}
				isosurface_nodes[i][j]->value = sum;
			}
		}

		double iso_value = 0.000003;
		int no_of_squares = (isosurface_nodes.size() * isosurface_nodes[0].size()) - (isosurface_nodes.size() + (isosurface_nodes[0].size() - 1));
		lines_start.clear();
		lines_end.clear();
		iso_points.clear();
		for (int i = 0; i < isosurface_nodes.size() - 1; i++) {
			for (int j = 0; j < isosurface_nodes.size() - 1; j++) {

				Node* n00 = isosurface_nodes[i][j];
				Node* n10 = isosurface_nodes[i + 1][j];
				Node* n01 = isosurface_nodes[i][j + 1];
				Node* n11 = isosurface_nodes[i + 1][j + 1];

				vector<DSC2D::vec2> points;

				if (n00->value > 0.0 || n10->value > 0.0) {
					if ((n00->value <= iso_value && n10->value > iso_value) ||
						(n00->value > iso_value && n10->value <= iso_value)) {
						//cout << n00->value << endl;
						DSC2D::vec2 interpolated_pos = interpolate_iso_nodes(n00, n10, iso_value);
						
						points.push_back(interpolated_pos); // Middle point
					}
				}
				if (n10->value > 0.0 || n11->value > 0.0) {
					if ((n10->value <= iso_value && n11->value > iso_value) ||
						(n10->value > iso_value && n11->value <= iso_value)) {

						DSC2D::vec2 interpolated_pos = interpolate_iso_nodes(n10, n11, iso_value);

						points.push_back(interpolated_pos); // Middle point
					}
				}
				if (n00->value > 0.0 || n01->value > 0.0) {
					if ((n00->value <= iso_value && n01->value > iso_value) ||
						(n00->value > iso_value && n01->value <= iso_value)) {

						DSC2D::vec2 interpolated_pos = interpolate_iso_nodes(n00, n01, iso_value);

						points.push_back(interpolated_pos); // Middle point
					}
				}
				if (n01->value > 0.0 || n11->value > 0.0) {
					if ((n01->value <= iso_value && n11->value > iso_value) ||
						(n01->value > iso_value && n11->value <= iso_value)) {

						DSC2D::vec2 interpolated_pos = interpolate_iso_nodes(n01, n11, iso_value);

						points.push_back(interpolated_pos); // Middle point
					}
				}
				if (points.size() == 2 || points.size() == 3) {
					lines_start.push_back(points[0]);
					lines_end.push_back(points[1]);
				}
				for each (DSC2D::vec2 p in points)
				{
					iso_points.push_back(p);
				}
			}
		}
	}
	/*Node* n_1 = new Node(DSC2D::vec2(0.0), 0);
	Node* n_2 = new Node(DSC2D::vec2(15.0, 0.0), 0.0000009);
	cout << interpolate_iso_nodes(n_1, n_2, 0.0000003) << endl;
	delete n_1;
	delete n_2;*/
}

DSC2D::vec2 SPH::interpolate_iso_nodes(Node* n_a, Node* n_b, double iso_value) {

	double l_between_nodes = 15.0;
	double t = l_between_nodes * ((iso_value - n_b->value) / (n_a->value - n_b->value));

	DSC2D::vec2 n_a_to_n_b = n_b->pos - n_a->pos;
	n_a_to_n_b = CGLA::normalize(n_a_to_n_b);
	DSC2D::vec2 interpolated_pos = n_a->pos + (n_a_to_n_b * (l_between_nodes - t));
	return interpolated_pos;
}

vector<Particle> SPH::get_close_particles_to_pos(DSC2D::vec2 pos, double radius) {
	Particle tmp_particle = Particle(pos, -1);
	vector<Particle> close_particles = get_close_particles(&tmp_particle, radius);
	return close_particles;
}

Particle SPH::get_particle(int id) {
	return particle_system.get_particle(id);
}

Particle* SPH::get_particle_ptr(int id) {
	return particle_system.get_particle_ptr(id);
}

double SPH::calculate_density(Particle p, vector<Particle> close_particles) {
	double density = 0.0;
	for each(Particle close_particle in close_particles) {
		DSC2D::vec2 p_to_p = p.pos - close_particle.pos;
		density += (close_particle.mass * poly6_kernel(p_to_p.length(), KERNEL_RADIUS));
	}
	return density + p.mass * poly6_kernel(0.0, KERNEL_RADIUS);
}

double SPH::calculate_a(Particle p, vector<Particle> close_particles) {
	DSC2D::vec2 v = DSC2D::vec2(0.0);
	double l = 0.0;
	for each(Particle close_particle in close_particles) {
		DSC2D::vec2 p_to_p = p.pos - close_particle.pos;
		DSC2D::vec2 tmp = close_particle.mass * gradient_kernel(p_to_p.length(), p_to_p, KERNEL_RADIUS);
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

double SPH::calculate_local_pressure(Particle p) {
	return GAS_CONSTANT * (p.density - REST_DENSITY);
}

DSC2D::vec2 SPH::calculate_surface_tension(Particle p, vector<Particle> close_particles) {
	DSC2D::vec2 inward_normal = DSC2D::vec2(0.0);
	double lap_c = 0.0;
	for each(Particle close_particle in close_particles) {
		DSC2D::vec2 p_to_p = p.pos - close_particle.pos;
		inward_normal += (close_particle.mass / close_particle.density) * gradient_kernel(p_to_p.length(), p_to_p, KERNEL_RADIUS);
		lap_c += (close_particle.mass / close_particle.density) * laplacian_kernel(p_to_p.length(), KERNEL_RADIUS);
	}
	DSC2D::vec2 surface_tension = DSC2D::vec2(0.0);
	if (inward_normal.length() > 0.035) {
		inward_normal.normalize();
		surface_tension = -400.0 * inward_normal * lap_c;
	}
	return surface_tension;
}

void SPH::correct_density_error() {

	double density_error = avg_density - REST_DENSITY;
	int max_iterations = 10;
	int iterations = 0;
	while (density_error > 0.00006 && iterations <= max_iterations) {
		double sum_density = 0.0;
		for (int i = 0; i < get_no_of_particle(); i++) {
			Particle* p_ptr = get_particle_ptr(i);
			DSC2D::vec2 next_pos = p_ptr->pos * (p_ptr->vel * get_delta());
			vector<Particle> close_particles = get_close_particles_to_pos(next_pos, KERNEL_RADIUS);
			double d = calculate_density(*p_ptr, close_particles);
			sum_density += d;
			p_ptr->density = d;
		}

		for (int i = 0; i < get_no_of_particle(); i++) {
			Particle* p_ptr = get_particle_ptr(i);
			double k_i = (p_ptr->density - REST_DENSITY) * (1 / (get_delta()*get_delta())) * p_ptr->a;
			vector<Particle> close_particles = get_close_particles(p_ptr, KERNEL_RADIUS);
			DSC2D::vec2 tmp = DSC2D::vec2(0.0);
			for each(Particle close_particle in close_particles) {
				DSC2D::vec2 p_to_p = p_ptr->pos - close_particle.pos;
				double k_j = (close_particle.density - REST_DENSITY) * (1 / get_delta() * get_delta()) * close_particle.a;
				if( p_ptr->density != 0.0 && close_particle.density != 0.0)
					tmp += close_particle.mass * (( k_i / p_ptr->density ) + (k_j / close_particle.density)) * gradient_kernel(p_to_p.length(), p_to_p, KERNEL_RADIUS);
			}
			p_ptr->vel = p_ptr->vel - get_delta() * tmp;
		}

		avg_density = sum_density / get_no_of_particle();
		density_error = avg_density - REST_DENSITY;
		iterations++;
	}
}

void SPH::correct_divergence_error() {
	int max_iterations = 12;
	int iterations = 0;
	DSC2D::vec2 divergence_avg = DSC2D::vec2(0.0);

	for (int i = 0; i < get_no_of_particle(); i++) {
		Particle* p_ptr = get_particle_ptr(i);

		vector<Particle> close_particles = get_close_particles(p_ptr);

		DSC2D::vec2 den_change = DSC2D::vec2(0.0);

		for each (Particle close_particle in close_particles)
		{
			DSC2D::vec2 p_to_p = p_ptr->pos - close_particle.pos;
			den_change += close_particle.mass * (p_ptr->vel - close_particle.vel) * gradient_kernel(p_to_p.length(), p_to_p, KERNEL_RADIUS);
		}

		//double density_change = -abs(p_ptr->density - p_ptr->old_density);
		p_ptr->density_divergence = den_change;

		divergence_avg += den_change;

	}

	if (get_no_of_particle() > 0)
		divergence_avg = divergence_avg / get_no_of_particle();

	while (divergence_avg.length() > 0.00001 && iterations <= max_iterations) {

		for (int i = 0; i < get_no_of_particle(); i++) {
			Particle* p_ptr = get_particle_ptr(i);
			DSC2D::vec2 k_i = (1.0 / get_delta()) * p_ptr->density_divergence * p_ptr->a;
			vector<Particle> close_particles = get_close_particles(p_ptr);
			DSC2D::vec2 tmp = DSC2D::vec2(0.0);
			for each(Particle close_particle in close_particles) {
				DSC2D::vec2 p_to_p = p_ptr->pos - close_particle.pos;

				DSC2D::vec2 k_j = (1.0 / get_delta()) * close_particle.density_divergence * close_particle.a;

				if (p_ptr->density != 0.0 && close_particle.density != 0.0) {
					tmp += close_particle.mass * ((k_i / p_ptr->density) + (k_j / close_particle.density)) * gradient_kernel(p_to_p.length(), p_to_p, KERNEL_RADIUS);
				}

			}
			if (!isnan(tmp[0]) && !isnan(tmp[1])) {
				p_ptr->vel = p_ptr->vel - (get_delta() * tmp);
			}
			else {
				cout << "NAN" << endl;
			}
		}
		divergence_avg = DSC2D::vec2(0.0);
		for (int i = 0; i < get_no_of_particle(); i++) {
			Particle* p_ptr = get_particle_ptr(i);

			vector<Particle> close_particles = get_close_particles(p_ptr);

			DSC2D::vec2 den_change = DSC2D::vec2(0.0);

			for each (Particle close_particle in close_particles)
			{
				DSC2D::vec2 p_to_p = p_ptr->pos - close_particle.pos;
				den_change += close_particle.mass * (p_ptr->vel - close_particle.vel) * gradient_kernel(p_to_p.length(), p_to_p, KERNEL_RADIUS);
			}

			p_ptr->density_divergence = den_change;

			divergence_avg += p_ptr->density_divergence;
		}
		if (get_no_of_particle() > 0)
			divergence_avg = divergence_avg / get_no_of_particle();
		iterations++;
	}

	/*for (int i = 0; i < get_no_of_particle(); i++) {
		Particle* p_ptr = get_particle_ptr(i);

		vector<Particle> close_particles = get_close_particles(p_ptr);

		DSC2D::vec2 den_change = DSC2D::vec2(0.0);

		for each (Particle close_particle in close_particles)
		{
			DSC2D::vec2 p_to_p = p_ptr->pos - close_particle.pos;
			den_change += close_particle.mass * (p_ptr->vel - close_particle.vel) * gradient_kernel(p_to_p.length(), p_to_p, KERNEL_RADIUS);
		}
		
		den_change = den_change;
		p_ptr->density_divergence = den_change;
		if (den_change.length() > 0.00001) {
			p_ptr->is_inside = false;
		}
		else {
			p_ptr->is_inside = true;
		}
	}*/
}

DSC2D::vec2 SPH::calculate_inward_normal(DSC2D::vec2 pos) {
	DSC2D::vec2 inward_normal = DSC2D::vec2(0.0);

	vector<Particle> close_particles = get_close_particles_to_pos(pos, 60.0);
	double lap_c = 0.0;
	for each(Particle close_particle in close_particles) {
		DSC2D::vec2 p_to_p = pos - close_particle.pos;
		if (close_particle.density != 0.0 && p_to_p.length() != 0.0) {
			inward_normal += (close_particle.mass / close_particle.density) * gradient_kernel(p_to_p.length(), p_to_p, 60.0);
			lap_c += (close_particle.mass / close_particle.density) * laplacian_kernel(p_to_p.length(), 60.0);
		}
	}
	return inward_normal;
}

DSC2D::vec2 SPH::calculate_pressure(Particle p, vector<Particle> close_particles) {

	DSC2D::vec2 pressure = DSC2D::vec2(0.0);
	for each(Particle close_particle in close_particles) {
		DSC2D::vec2 p_to_p = p.pos - close_particle.pos;
		if (close_particle.density != 0.0) {
			pressure += close_particle.mass * ((p.local_pressure + close_particle.local_pressure) / (2.0 * close_particle.density)) * gradient_kernel(p_to_p.length(), p_to_p, KERNEL_RADIUS);		
		}
	}
	return -pressure;
}

DSC2D::vec2 SPH::calculate_viscocity(Particle p, vector<Particle> close_particles) {

	DSC2D::vec2 viscocity = DSC2D::vec2(0.0);
	for each(Particle close_particle in close_particles) {
		DSC2D::vec2 p_to_p = p.pos - close_particle.pos;
		if(close_particle.density != 0.0)
			viscocity += p.mass * ((close_particle.vel - p.vel) / close_particle.density) * laplacian_kernel(p_to_p.length(), KERNEL_RADIUS);
	}
	return VISCOCITY_TERM * viscocity;

}

DSC2D::vec2 SPH::calculate_external_forces(Particle p) {
	DSC2D::vec2 gravity = DSC2D::vec2(0.0, -9.4);
	return gravity;
}

void SPH::update_velocity(double delta_time) {
	for (int i = 0; i < particle_system.get_number_of_particles(); i++) {
		Particle* p_ptr = particle_system.get_particle_ptr(i);
		if (!p_ptr->is_fixed) {
			
			DSC2D::vec2 acceleration = p_ptr->pressure + p_ptr->viscocity + p_ptr->external_forces;
			
			DSC2D::vec2 velocity = p_ptr->vel + acceleration;
			
			p_ptr->vel = velocity;
			
			p_ptr->vel += check_collision(*p_ptr, delta_time);
		}
	}
}

void SPH::update_position(double delta_time) {
	
	for (int i = 0; i < particle_system.get_number_of_particles(); i++) {
		Particle* p_ptr = particle_system.get_particle_ptr(i);
		
		if(!p_ptr->is_fixed)
			p_ptr->pos += p_ptr->vel * delta_time;
		p_ptr->pos_0 = p_ptr->pos / 2.0;
	}
}

void SPH::draw_collision_boxes() {
	glBegin(GL_QUADS);
	for each (Collision_Box c_b in collision_boxes)
	{
		DSC2D::vec2 pos = c_b.pos;
		float half_height = c_b.height * 0.5;
		float half_width = c_b.width * 0.5;
		DSC2D::vec3 color = c_b.color;
		glColor3f(color[0]/255.0, color[1]/255.0, color[2]/255.0);
		glVertex2f((pos[0] - half_width)*scale, (pos[1] - half_height)*scale);
		glVertex2f((pos[0] + half_width)*scale, (pos[1] - half_height)*scale);
		glVertex2f((pos[0] + half_width)*scale, (pos[1] + half_height)*scale);
		glVertex2f((pos[0] - half_width)*scale, (pos[1] + half_height)*scale);
		
	}
	glEnd();
}

void SPH::draw_particles() {
	particle_system.draw();
}

void SPH::create_particle_at_mouse_pos() {
	if (get_no_of_particle() >= 900)
		return;
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
			particle_system.create_particle(DSC2D::vec2(x, y) * reverse_scale);
		}
	}
}

int SPH::create_particle(DSC2D::vec2 pos) {

	return particle_system.create_particle(pos);

}

void SPH::create_collision_box(DSC2D::vec2 pos, float height, float width, DSC2D::vec3 c) {
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
			DSC2D::vec3 grey = DSC2D::vec3(100, 100, 100);
			create_collision_box(DSC2D::vec2(x, y) * reverse_scale, 100.0, 100.0, grey);
		}
	}
}

DSC2D::vec2 SPH::check_collision(Particle p, double delta_time) {
	DSC2D::vec2 new_pos = p.pos + (p.vel * delta_time);
	DSC2D::vec2 impulse = DSC2D::vec2(0.0);
	Particle* p_ptr = particle_system.get_particle_ptr(p.id);

	for each(Collision_Box c_b in collision_boxes) {
	
		if (c_b.is_inside(p.pos)) {
			// particle is inside collision box
			DSC2D::vec2 wall_normal = DSC2D::vec2(0.0, 0.0);
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
				wall_normal = DSC2D::vec2(0.0, 1.0);
				break;
			case 1:
				wall_normal = DSC2D::vec2(1.0, 0.0);
				break;
			case 2:
				wall_normal = DSC2D::vec2(0.0, -1.0);
				break;
			case 3:
				wall_normal = DSC2D::vec2(-1.0, 0.0);
			default:
				break;
			}

			if (c_b.is_inside(p.pos)) {

				DSC2D::vec2 tmp = (wall_normal) * ((c_b.height)*0.5) + c_b.pos;
				Particle* p_ptr = particle_system.get_particle_ptr(p.id);

				if (wall_id == 1 || wall_id == 3) {
					tmp = wall_normal * ((c_b.width)*0.5) + c_b.pos;
					
					p_ptr->pos[0] = tmp[0];
					
				}
				else {
				
					p_ptr->pos[1] = tmp[1];
					
				}	
			}

			// teleport to edge of box if 
			if (c_b.is_inside(new_pos)) {
				impulse += calculate_collision_impulse(p, wall_normal);
			}
			
			if (closest_distance <= 5.0) {
				// If particle very close to collider, slow it down
				impulse -= p.vel * 0.1;
			}
		}
	}
	return impulse;
}

DSC2D::vec2 SPH::calculate_collision_impulse(Particle p, DSC2D::vec2 wall_normal) {

	double normal_speed = dot(wall_normal, p.vel);
	double j = -(1.0 + COEFFICIENT_OF_RESTITUTION) * p.mass * normal_speed;
	DSC2D::vec2 tmp = (-p.vel - (normal_speed * wall_normal));
	if (tmp.length() > 0.001) {
		tmp = normalize(tmp);
	}
	DSC2D::vec2 J = j * wall_normal - COEFFICIENT_OF_FRICTION * j * tmp;

	return J;
}

vector<double> SPH::signed_distance_to_walls(DSC2D::vec2 p, Collision_Box c_b) {
	
	DSC2D::vec2 pos = c_b.pos;
	float height = c_b.height;
	float width = c_b.width;

	vector<double> distances;
	// upper wall
	DSC2D::vec2 wall_point(pos[0] - (width*0.5), pos[1] + (height*0.5));
	DSC2D::vec2 wall_dir(1.0, 0.0);
	double distance = signed_distance(p, wall_point, wall_dir);
	distances.push_back(distance);
	// Rightmost wall
	wall_point = DSC2D::vec2(pos[0] + (width*0.5), pos[1] - (height*0.5));
	wall_dir = DSC2D::vec2(0.0, 1.0);
	distance = signed_distance(p, wall_point, wall_dir);
	distances.push_back(distance);
	// Lower wall
	wall_point = DSC2D::vec2(pos[0] - (width*0.5), pos[1] - (height*0.5));
	wall_dir = DSC2D::vec2(1.0, 0.0);
	distance = signed_distance(p, wall_point, wall_dir);
	distances.push_back(distance);
	// Leftmost wall
	wall_point = DSC2D::vec2(pos[0] - (width*0.5), pos[1] - (height*0.5));
	wall_dir = DSC2D::vec2(0.0, 1.0);
	distance = signed_distance(p, wall_point, wall_dir);
	distances.push_back(distance);

	return distances;
}

double SPH::signed_distance(DSC2D::vec2 p, DSC2D::vec2 plane_point, DSC2D::vec2 plane_dir) {

	DSC2D::vec2 w = p - plane_point;
	DSC2D::vec2 u = normalize(plane_dir); 
	DSC2D::vec2 d = w - (dot(w, u) * u);
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
					glVertex2f(p.pos_0[0] + cos(degInRad)*KERNEL_RADIUS, p.pos_0[1] + sin(degInRad)*KERNEL_RADIUS);
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
				DSC2D::vec2 pos = p.pos_0;
				glVertex2f(pos[0], pos[1]);
				glVertex2f(pos[0] + p.pressure[0], pos[1] + p.pressure[1]);
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
				DSC2D::vec2 pos = p.pos_0;
				glVertex2f(pos[0], pos[1]);
				glVertex2f(pos[0] + p.viscocity[0], pos[1] + p.viscocity[1]);
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
				DSC2D::vec2 pos = p.pos_0;
				glVertex2f(pos[0], pos[1]);
				glVertex2f(pos[0] + p.surface_tension[0], pos[1] + p.surface_tension[1]);
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
				DSC2D::vec2 pos = p.pos_0;
				glVertex2f(pos[0], pos[1]);
				glVertex2f(pos[0] + (p.vel[0] * 0.15), pos[1] + (p.vel[1] * 0.15));
			}
		}
		glEnd();
	}
}

void SPH::draw_isosurface() {
	/*glPointSize(5.0);
	glBegin(GL_POINTS);
	glColor3f(0.0f, 1.0f, 0.0f);
	for (int i = 0; i < isosurface_nodes.size(); i++) {
		for (int j = 0; j < isosurface_nodes[i].size(); j++) {
			if (isosurface_nodes[i][j]->value > 0.0) {
				DSC2D::vec2 pos = isosurface_nodes[i][j]->pos * scale;
				glVertex2f(pos[0], pos[1]);
			}
		}
	}
	glEnd();*/
	/*if (is_drawing_isosurfacee) {
		glBegin(GL_LINES);
		glColor3f(0.0f, 1.0f, 0.0f);
		for (int i = 0; i < lines_start.size(); i++) {
			glVertex2f(lines_start[i][0] * scale, lines_start[i][1] * scale);
			glVertex2f(lines_end[i][0] * scale, lines_end[i][1] * scale);
		}

		glEnd();
	}*/
	if (is_drawing_isosurfacee) {
		glBegin(GL_POINTS);
		glColor3f(0.0f, 1.0f, 0.0f);
		for (int i = 0; i < iso_points.size(); i++) {
			glVertex2f(iso_points[i][0] * scale, iso_points[i][1] * scale);
		}

		glEnd();
	}
	
}

void SPH::draw_dsc_velocities(DSC2D::DeformableSimplicialComplex& dsc) {
	if (is_drawing_dsc_vel) {
		DSC2D::DeformableSimplicialComplex* dsc_ptr = &dsc;
		glBegin(GL_LINES);

		glColor3f(1.0f, 0.0f, 0.0f);
		for (auto vi = dsc_ptr->vertices_begin(); vi != dsc_ptr->vertices_end(); ++vi) {

			DSC2D::vec2 pos = dsc_ptr->get_pos(*vi);
			DSC2D::vec2 vel = dsc_vert_velocities[*vi];

			glVertex2f(pos[0], pos[1]);
			glVertex2f(pos[0] + (vel[0] * 0.25), pos[1] + (vel[1] * 0.25));
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
			DSC2D::vec2 translation = (DSC2D::vec2(x, y)*reverse_scale) - collision_boxes[id].pos;
			collision_boxes[id].pos = DSC2D::vec2(x, y)*reverse_scale;
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

DSC2D::vec2 SPH::gradient_kernel(float r, DSC2D::vec2 v, float d) {
	if (r < 0 || r > d) {
		return DSC2D::vec2(0);
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

int SPH::get_no_of_particle() {
	return particle_system.get_number_of_particles();
}

DSC2D::vec2 SPH::get_particle_pos(int id) {
	return particle_system.get_particle(id).pos;
}

void SPH::project_velocities(DSC2D::DeformableSimplicialComplex& dsc) {
	
	DSC2D::DeformableSimplicialComplex* dsc_ptr = &dsc;

	/*for (auto vi = dsc_ptr->vertices_begin(); vi != dsc_ptr->vertices_end(); ++vi) {
		dsc_vert_velocities[*vi] = DSC2D::vec2(0.0);
	}

	for (auto fi = dsc_ptr->faces_begin(); fi != dsc_ptr->faces_end(); ++fi) {

		vector<HMesh::VertexID> verts = dsc_ptr->get_verts(*fi);
		DSC2D::vec2 p0 = dsc_ptr->get_pos(verts[0]);
		DSC2D::vec2 p1 = dsc_ptr->get_pos(verts[1]);
		DSC2D::vec2 p2 = dsc_ptr->get_pos(verts[2]);

		double area = 0.5 * (-p1[1] * p2[0] + p0[1] * (-p1[0] + p2[0]) + p0[0] * (p1[1] - p2[1]) + p1[0] * p2[1]);

		
		for (int j = 0; j < 3;j++) {
			DSC2D::vec2 pos = dsc_ptr->get_pos(*fi)[j];
			vector<Particle> close_particles = get_close_particles_to_pos(pos);
			int no_of_close_particles = close_particles.size();
			int count = 0;
			double sum_weight = 0.0;


			for (int i = 0; i < no_of_close_particles; i++) {
				DSC2D::vec2 p = close_particles[i].pos;
				if (in_triangle(p, p0, p1, p2)) {
					DSC2D::vec2 r = pos - p;
					sum_weight += poly6_kernel(r.length(),15.0);
				}
			
			}

			for (int i = 0; i < no_of_close_particles; i++) {
				DSC2D::vec2 p = close_particles[i].pos;
				if (in_triangle(p, p0, p1, p2)) {
					DSC2D::vec2 r = pos - close_particles[i].pos;
					Particle* p_ptr = particle_system.get_particle_ptr(close_particles[i].id);
					dsc_vert_velocities[verts[j]] += (p_ptr->vel * (poly6_kernel(r.length(),15.0) / sum_weight));
				}
			}
		}
	}*/

	for (auto vi = dsc_ptr->vertices_begin(); vi != dsc_ptr->vertices_end(); ++vi) {

		dsc_vert_velocities[*vi] = DSC2D::vec2(0.0);

		DSC2D::vec2 velocity = DSC2D::vec2(0.0);
		DSC2D::vec2 pos = dsc.get_pos(*vi);


		vector<Particle> close_particles = get_close_particles_to_pos(pos, 30.0);
		int no_of_close_particles = close_particles.size();
		double sum_weight = 0.0;
		for (int i = 0; i < no_of_close_particles; i++) {
			DSC2D::vec2 r = pos - close_particles[i].pos;
			sum_weight += poly6_kernel(r.length(), 30.0);

		}
	
		for (int i = 0; i < no_of_close_particles; i++) {
			DSC2D::vec2 r = pos - close_particles[i].pos;
			Particle* p_ptr = particle_system.get_particle_ptr(close_particles[i].id);

			velocity += p_ptr->vel * (poly6_kernel(r.length(), 30.0) / sum_weight);
		}
		if (!isnan(velocity[0])) {
			dsc_vert_velocities[*vi] = velocity;
		}
	}
}

void SPH::project_velocities_to_particles(DSC2D::DeformableSimplicialComplex& dsc) {

	DSC2D::DeformableSimplicialComplex* dsc_ptr = &dsc;

	for (auto fi = dsc_ptr->faces_begin(); fi != dsc_ptr->faces_end(); ++fi) {
		
		vector<HMesh::VertexID> verts = dsc_ptr->get_verts(*fi);
		DSC2D::vec2 p0 = dsc_ptr->get_pos(verts[0]);
		DSC2D::vec2 p1 = dsc_ptr->get_pos(verts[1]);
		DSC2D::vec2 p2 = dsc_ptr->get_pos(verts[2]);

		double area = 0.5 * (-p1[1] * p2[0] + p0[1] * (-p1[0] + p2[0]) + p0[0] * (p1[1] - p2[1]) + p1[0] * p2[1]);

		DSC2D::vec2 pos = dsc_ptr->get_pos(*fi)[0];
		vector<Particle> close_particles = get_close_particles_to_pos(pos);

		for (int i = 0; i < close_particles.size(); i++) {
			DSC2D::vec2 p = close_particles[i].pos;
			if (in_triangle(p, p0, p1, p2)) {
				double area_0 = 0.5 * (-p1[1] * p2[0] + p[1] * (-p1[0] + p2[0]) + p[0] * (p1[1] - p2[1]) + p1[0] * p2[1]);
				double area_1 = 0.5 * (-p[1] * p2[0] + p0[1] * (-p[0] + p2[0]) + p0[0] * (p[1] - p2[1]) + p[0] * p2[1]);
				double area_2 = 0.5 * (-p1[1] * p[0] + p0[1] * (-p1[0] + p[0]) + p0[0] * (p1[1] - p[1]) + p1[0] * p[1]);
				Particle* p_ptr = particle_system.get_particle_ptr(close_particles[i].id);
				DSC2D::vec2 velocity = dsc_vert_velocities[verts[0]] * (area_0 / area) + dsc_vert_velocities[verts[1]] * (area_1 / area) + dsc_vert_velocities[verts[2]] * (area_2 / area);
				
				if (!isnan(velocity[0]))
					p_ptr->vel = velocity;
			}
		}
	}

}

Particle SPH::get_closest_particle(DSC2D::vec2 pos) {
	Particle tmp_particle = Particle(pos, -1);
	Particle closest_particle = Particle(pos, -1);
	particle_system.get_closest_particle(tmp_particle, closest_particle, 1000.0);
	return closest_particle;
}

bool SPH::in_triangle(DSC2D::vec2 p, DSC2D::vec2 p0, DSC2D::vec2 p1, DSC2D::vec2 p2) {
	double area = 0.5 * (-p1[1] * p2[0] + p0[1] * (-p1[0] + p2[0]) + p0[0] * (p1[1] - p2[1]) + p1[0] * p2[1]);
	double s = 1 / (2 * area)*(p0[1] * p2[0] - p0[0] * p2[1] + (p2[1] - p0[1])*p[0] + (p0[0] - p2[0])*p[1]);
	double t = 1 / (2 * area)*(p0[0] * p1[1] - p0[1] * p1[0] + (p0[1] - p1[1])*p[0] + (p1[0] - p0[0])*p[1]);
	double u = 1 - s - t;
	if (s > 0 && t > 0 && u > 0) {
		// is inside triangle
		return true;
	}
	return false;
}