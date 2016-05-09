#include "track_particles_function.h"


void track_particles_function::deform(DSC::DeformableSimplicialComplex<>& dsc) {
	// TEST // 
	for (auto n = dsc.nodes_begin(); n != dsc.nodes_end(); n++) {

		if (dsc.is_movable(n.key())) {
			vec3 to_center = vec3(0.0, 0.0, -1.0);
			to_center.normalize();
			if(dsc.get_pos(n.key())[2] > -0.9)
				dsc.set_destination(n.key(), dsc.get_pos(n.key()) + (to_center * 0.01));

		}
	}

	dsc.deform();

}