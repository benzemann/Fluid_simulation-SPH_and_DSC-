#include "track_particles_function.h"


void track_particles_function::deform(DSC::DeformableSimplicialComplex<>& dsc) {
	
	
	vector<vec3> iso_points = sph->get_iso_points();
	for (auto n = dsc.nodes_begin(); n != dsc.nodes_end(); n++) {
		if (dsc.is_movable(n.key())) {

			vec3 closest_point = vec3(1000.0);
			double closest_distance = 1000.0;

			for each (vec3 point in iso_points)
			{
				vec3 n_to_point = point - dsc.get_pos(n.key());
				double distance = n_to_point.length();
				if (distance < closest_distance) {
					closest_distance = distance;
					closest_point = point;
				}
			}


			vec3 vel = closest_point - dsc.get_pos(n.key());
			//cout << vel.length() << endl;
			//if (vel.length() > 0.01) {
			//	vel.normalize();
			//	vel *= 0.01;
			//}
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
	
	for (int i = 0; i < 1; i++) {
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
	}

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