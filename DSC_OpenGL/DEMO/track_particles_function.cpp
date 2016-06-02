#include "track_particles_function.h"


void track_particles_function::deform(DSC::DeformableSimplicialComplex<>& dsc) {
	for (auto n = dsc.nodes_begin(); n != dsc.nodes_end(); n++) {
		if (dsc.is_movable(n.key())) {
			double radius = 100.0;
			vec3 vert_pos = dsc.get_pos(n.key());
			vert_pos = vec3(vert_pos[0] + 1, vert_pos[1] + 1, vert_pos[2] + 1);
			vert_pos = vert_pos * sph->get_down_scale();
			vector<Particle> close_particles = sph->get_close_particles_to_pos(vert_pos, radius);
			vec3 vel = vec3(0.0);

			if (close_particles.size() > 0) {

				double iso_value = 0.00004;
				Particle p_tmp = Particle(vert_pos, -1);
				double density = sph->calculate_density(p_tmp, close_particles, radius);
				//cout << density << endl;
				vec3 den_grad = sph->calculate_density_gradient(p_tmp, close_particles, radius);
				double len = den_grad.length();
				vel = -den_grad * ((density - iso_value) / (len * len));
				vel = vel / sph->get_down_scale();
				//vel = vec3(-1 + vel[0], -1 + vel[1], -1 + vel[2]);
				//cout << density << endl;
				if (vel.length() > 0.5) {
					vel.normalize();
					vel = vel * 0.5;
				}
			}
			else {

				Particle closest_p = sph->get_closest_particle(vert_pos);
				vel = closest_p.pos - vert_pos;
				vel.normalize();

			}

			dsc.set_destination(n.key(), dsc.get_pos(n.key()) + vel);
		}
	}
	
	// TEST // 
	/*for (auto n = dsc.nodes_begin(); n != dsc.nodes_end(); n++) {

		if (dsc.is_movable(n.key())) {
			vec3 to_center = vec3(0.0, 0.0, -1.0);
			to_center.normalize();
			if(dsc.get_pos(n.key())[2] > -0.9)
				dsc.set_destination(n.key(), dsc.get_pos(n.key()) + (to_center * 0.01));

		}
	}*/

	/*for (int i = 0; i < 1; i++) {
		for (auto nit = dsc.nodes_begin(); nit != dsc.nodes_end(); nit++) {
			if (dsc.is_movable(nit.key()))
				dsc.laplacian(nit.key());
		}
	}*/
	

	dsc.deform();
	
	/*for (int i = 0; i < 1; i++) {
		vec3 new_pos, p;
		for (auto nit = dsc.nodes_begin(); nit != dsc.nodes_end(); nit++)
		{
			if (dsc.is_movable(nit.key()))
			{
				//dsc.laplacian(nit.key());
				p = nit->get_pos();
				new_pos = p + 0.9*(dsc.get_barycenter(nit.key(), true) - p);
				dsc.set_destination(nit.key(), new_pos);

			}
		}

		dsc.deform();
	}*/

	/*
	for (auto nit = dsc.nodes_begin(); nit != dsc.nodes_end(); nit++)
	{
		if (dsc.is_movable(nit.key()))
		{
			p = nit->get_pos();
			new_pos = p + 0.5 * (dsc.get_barycenter(nit.key(), true) - p);
			dsc.set_destination(nit.key(), new_pos);

		}
	}

	dsc.deform();

	
	for (auto nit = dsc.nodes_begin(); nit != dsc.nodes_end(); nit++)
	{
		if (dsc.is_movable(nit.key()))
		{
			p = nit->get_pos();
			new_pos = p + 0.5 * (dsc.get_barycenter(nit.key(), true) - p);
			dsc.set_destination(nit.key(), new_pos);

		}
	}
	
	dsc.deform();
	*/
}