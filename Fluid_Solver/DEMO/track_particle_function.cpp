#include "track_particle_function.h"
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "object_generator.h"
#include <iostream>
#include <fstream>

using namespace std;
using namespace Eigen;
using namespace DSC2D;

using EigMat = MatrixXd;
using EigVec = VectorXd;

HMesh::VertexAttributeVector<int> particle_vertices;

void track_particle_function::init(DSC2D::DeformableSimplicialComplex& dsc) {
	cout << "init" << endl;

	
	// initialize particles at all vertices marked "inside"
	for (auto vi = dsc.vertices_begin(); vi != dsc.vertices_end(); vi++) {
		//Particle* tmp_par = new Particle();
		//tmp_par->id = -1;
		//particle_vertices[*vi] = tmp_par;
		if (!dsc.is_outside(*vi)) {
			int i = sph.create_particle(dsc.get_pos(*vi));
			//particle_vertices[*vi] = sph.get_particle_ptr(i);
			p_verts.push_back(i);
			Particle* p_ptr = sph.get_particle_ptr(i);
			
			p_ptr->pos = dsc.get_pos(*vi);
			particle_vertices[*vi] = i;
			dsc.set_destination(*vi, dsc.get_pos(*vi));
			
		}
	}

}

void track_particle_function::deform(DSC2D::DeformableSimplicialComplex& dsc) {

	HMesh::VertexAttributeVector<vec2> vels;
	for (auto vi = dsc.vertices_begin(); vi != dsc.vertices_end(); vi++) {
		vels[*vi] = vec2(0.0);
		if (dsc.is_interface(*vi)) {
			double radius = 100.0;
			vec2 vert_pos = dsc.get_pos(*vi) * 2.0;
			vector<Particle> close_particles = sph.get_close_particles_to_pos(vert_pos, radius);
			vec2 vel = vec2(0.0);
			if (close_particles.size() > 0) {
				double iso_value = 0.000016;
				Particle p_tmp = Particle(vert_pos, -1);
				double density = sph.calculate_density(p_tmp, close_particles, radius);
				//cout << density << endl;
				vec2 den_grad = sph.calculate_gradient_density(p_tmp, close_particles, radius);
				double len = den_grad.length();
				vel = -den_grad * ((density - iso_value) / (len * len));
				vel = vel * 0.5;
				
				if (vel.length() > 5.0) {
					vel.normalize();
					vel = vel * 5.0;
				}


				vels[*vi] = vel * 100.0;
			}
			else {

				Particle closest_p = sph.get_closest_particle(vert_pos);
				vel = closest_p.pos - vert_pos;
				vel.normalize();

			}
			sph.set_dsc_vert_velocities(vels);
			dsc.set_destination(*vi, dsc.get_pos(*vi) + vel);

			/*vec2 closest_point = vec2(-100.0);
			double closest_dis = 10000000.0;
			vec2 pos = dsc.get_pos(*vi);

			for each (vec2 point in iso_points)
			{
				vec2 p_to_point = (point*sph.get_scale()) - pos;
				double distance = p_to_point.length();
				if (distance < closest_dis) {
					closest_dis = distance;
					closest_point = (point*sph.get_scale());
				}
			}
			
			if (closest_point != vec2(-100.0)) {
				vec2 vel = closest_point - dsc.get_pos(*vi);
				//if (vel.length() > 1.5) {
				//vel.normalize();
				//vel = vel * 1.5;
				//}
				if (vel.length() > 3.0) {
					vel.normalize();
					vel = vel * 1.5;
				}
				else {
					vel.normalize();
					vel = vel * 0.1;
				}
				dsc.set_destination(*vi, dsc.get_pos(*vi) + vel);
				
			}*/
			/*for (int i = 0; i < isosurface_start.size(); i++) {

				// Find closest point to each line
				vec2 start_p = pos - (isosurface_start[i] * sph.get_scale());
				vec2 line = (isosurface_end[i]* sph.get_scale()) - (isosurface_start[i]* sph.get_scale());

				double line_line = CGLA::dot(line, line);
				double sp_line = CGLA::dot(start_p, line);

				double t = sp_line / line_line;

				if (t < 0.0)
					t = 0.0;
				else if (t > 1.0)
					t = 1.0;
				
				vec2 closest = (isosurface_start[i]* sph.get_scale()) + line * t;
				vec2 v = closest - pos;
				double l = v.length();
				if (l < closest_dis) {
					closest_dis = l;
					closest_point = closest;
				}
			}
			//cout << closest_dis << endl;
			if (closest_point != vec2(-100.0)) {
				vec2 vel = closest_point - dsc.get_pos(*vi);
				//if (vel.length() > 1.5) {
					//vel.normalize();
					//vel = vel * 1.5;
				//}
				if (vel.length() > 1.5) {
					vel.normalize();
					vel = vel * 1.5;
					dsc.set_destination(*vi, dsc.get_pos(*vi) + vel);
				}
			}*/
		}
	}

	dsc.deform();
	/*
	for (auto vi = dsc.vertices_begin(); vi != dsc.vertices_end(); vi++) {
		if (dsc.is_interface(*vi)) {
			HMesh::Walker w = dsc.walker(*vi);

			vector<DSC2D::vec2> neighbor_pos;

			while (!w.full_circle()) {
				HMesh::VertexID v = w.next().vertex();

				HMesh::HalfEdgeID he = w.halfedge();

				if (dsc.is_interface(he)) {
					neighbor_pos.push_back(dsc.get_pos(w.next().vertex()));
				}

				w = w.circulate_vertex_ccw();
			}

			DSC2D::vec2 smoothed_pos = dsc.get_pos(*vi) * 0.9;
			for each (DSC2D::vec2 p in neighbor_pos)
			{
				smoothed_pos += p * 0.05;
			}
			//smoothed_pos = smoothed_pos / 3.0;
			dsc.set_destination(*vi, smoothed_pos);
		}
	}

	dsc.deform();*/
}

void track_particle_function::deform_old_old(DSC2D::DeformableSimplicialComplex& dsc) {

	HMesh::VertexAttributeVector<DSC2D::vec2> velocities;
	HMesh::VertexAttributeVector<DSC2D::vec2> velocities_2;
	std::vector<HMesh::FaceID> face_ids;
	std::vector<HMesh::VertexID> verts_ids;
	for (auto fi = dsc.faces_begin(); fi != dsc.faces_end(); ++fi) {
		auto vis = dsc.get_verts(*fi);
		if (dsc.get_label(*fi) == 1) {
			face_ids.push_back(*fi);
		}
	}
	for (auto vi = dsc.vertices_begin(); vi != dsc.vertices_end(); vi++) {
		if (!dsc.is_outside(*vi)) {
			verts_ids.push_back(*vi);
			velocities[*vi] = vec2(0.0, 0.0);
			velocities_2[*vi] = vec2(0.0, 0.0);
		}
	}


	// ReMapping!
	/*int j = 0;
	for (auto vi = dsc.vertices_begin(); vi != dsc.vertices_end(); vi++) {
	if (!dsc.is_outside(*vi)) {
	if (particle_vertices[*vi] != NULL) {

	}
	else {
	//cout << "need map" << endl;
	Particle p = sph.get_closest_particle(dsc.get_pos(*vi));
	//particle_vertices[*vi] = p.id;
	}

	j++;
	}
	}*/

	//int tmp = 0;
	//int i = 0;
	/*for (int i = 0; i < sph.get_no_of_particle(); i++) {
	Particle* p_ptr = sph.get_particle_ptr(i);

	if (p_ptr->vel.length() < 150.0) {
	p_ptr->is_inside = true;
	}
	else {
	p_ptr->is_inside = false;
	}
	}*/
	//cout << get_curent_volume(dsc) / m_V0 << endl;
	for (int i = 0; i < sph.get_no_of_particle(); i++) {
		Particle* p_ptr = sph.get_particle_ptr(i);
		p_ptr->is_inside = false;
		/*
		if ( volume_d > 1.0 ) {
			vec2 n = p_ptr->surface_tension;
			n.normalize();

			p_ptr->vel += (n * 100.0);
		}
		else if (volume_d < 1.0) {
			vec2 n = p_ptr->surface_tension;
			n.normalize();
			p_ptr->vel += (-n * 100.0);
		}*/
	}
	/*for (auto fi = dsc.faces_begin(); fi != dsc.faces_end(); fi++) {
	for (int i = 0; i < sph.get_no_of_particle(); i++) {
	if (!dsc.is_outside(*fi)) {
	Particle* p_ptr = sph.get_particle_ptr(i);
	vector<DSC2D::vec2> verts = dsc.get_pos(*fi);
	if (sph.in_triangle(p_ptr->pos, verts[0], verts[1], verts[2]) && p_ptr->vel.length() < 150.0) {
	p_ptr->is_inside = true;
	}
	}
	}
	}*/
	for (auto vi = dsc.vertices_begin(); vi != dsc.vertices_end(); vi++) {
		particle_vertices[*vi] = -1;
		velocities[*vi] = vec2(0.0, 0.0);
		if (!dsc.is_outside(*vi)) {

			//if (particle_vertices[*vi] != NULL) {
			/*Particle closest_particle = Particle();
			closest_particle.pos = DSC2D::vec2(1000.0, 1000.0);
			double closest_dis = 10000.0;
			vector<Particle> close_particles = sph.get_close_particles_to_pos(dsc.get_pos(*vi));
			for each (Particle par in close_particles)
			{
			DSC2D::vec2 par_to_vert = par.pos - dsc.get_pos(*vi);
			double dis = par_to_vert.length();
			if (dis < closest_dis) {
			closest_dis = dis;
			closest_particle = par;
			}
			}*/
			Particle p = sph.get_closest_particle(dsc.get_pos(*vi));
			Particle* p_ptr = sph.get_particle_ptr(p.id);
			vec2 vert_to_p = p_ptr->pos - dsc.get_pos(*vi);
			velocities[*vi] = vec2(0.0);

			//if ((p_ptr->vel.length() < 300.0) && p_ptr->is_inside == false) {
			velocities[*vi] = p_ptr->vel;
			particle_vertices[*vi] = p_ptr->id;
			velocities_2[*vi] = p_ptr->pos - dsc.get_pos(*vi);
			p_ptr->is_inside = true;
				//velocities[*vi] = sph.calculate_inward_normal(p_ptr->pos) * 100.0;
			//}
			//cout << p_ptr->vel.length() << endl;
			//}
			//i++;
			//if (tmp == 0) {
			/*Particle* p_ptr = sph.get_particle_ptr(p_verts[i]);
			//cout << p_ptr->id << endl;
			if (p_ptr->id != -1 ) {
			if (!isnan(p_ptr->pos[0])) {
			dsc.set_destination(*vi, p_ptr->pos);
			}
			}

			//}
			i++;
			tmp++;*/
			//cout << p_ptr->pos << endl;
		}

	}
	//cout << c << " " << sph.get_no_of_particle() << endl;
	if (verts_ids.size() > 0) {
		//solve_for_incompressibility(face_ids, verts_ids, velocities, dsc);
	}

	for (auto vi = dsc.vertices_begin(); vi != dsc.vertices_end(); vi++) {
		if (!dsc.is_outside(*vi)) {
			//dsc.set_destination(*vi, dsc.get_pos(*vi) + velocities[*vi]);
			if (particle_vertices[*vi] != -1) {
				Particle* p_ptr = sph.get_particle_ptr(particle_vertices[*vi]);
				//p_ptr->vel = velocities[*vi];
				dsc.set_destination(*vi, p_ptr->pos + (velocities[*vi] * sph.get_delta()));
			}
		}
	}
	dsc.deform();
	//sph.set_dsc_vert_velocities(velocities);
	if (tmp == 0) {
		m_V0 = get_curent_volume(dsc);

		sph.set_delta(0.004);
	}
	double volume_d = get_curent_volume(dsc) / m_V0;
	if (tmp > 0 && ( volume_d > 1.03 || volume_d < 0.97)) {
		//solve_for_incompressibility_test(velocities_2, dsc);
		//dsc.deform();
	}
	//cout << get_curent_volume(dsc) / m_V0 << endl;
	/*for (auto vi = dsc.vertices_begin(); vi != dsc.vertices_end(); vi++) {
		if (!dsc.is_outside(*vi)) {
			//dsc.set_destination(*vi, dsc.get_pos(*vi) + velocities[*vi]);
			if (particle_vertices[*vi] != -1) {
				Particle* p_ptr = sph.get_particle_ptr(particle_vertices[*vi]);
			//	p_ptr->pos = dsc.get_pos(*vi);
				p_ptr->vel += velocities_2[*vi];
				//dsc.set_destination(*vi, dsc.get_pos(*vi) + (p_ptr->vel*sph.get_delta()));
			}
		}
	}*/
	
	//sph.set_volume_diff(m_V0);
	done(dsc);
	tmp++;
}


void track_particle_function::deform_old(DSC2D::DeformableSimplicialComplex& dsc) {
	DeformableSimplicialComplex* dsc_ptr = &dsc;
	SPH* sph_ptr = &sph;

	HMesh::VertexAttributeVector<DSC2D::vec2> velocities = sph_ptr->get_dsc_vert_velocities();
	std::vector<HMesh::FaceID> face_ids;
	std::vector<HMesh::VertexID> verts_ids;
	for (auto fi = dsc.faces_begin(); fi != dsc.faces_end(); ++fi) {
		auto vis = dsc.get_verts(*fi);
		if (dsc.get_label(*fi) == 1) {
			face_ids.push_back(*fi);
		}
	}
	for (auto vi = dsc.vertices_begin(); vi != dsc.vertices_end(); vi++) {
		if (!dsc_ptr->is_outside(*vi)) {
			verts_ids.push_back(*vi);
		}
	}


	for (auto vi = dsc_ptr->vertices_begin(); vi != dsc_ptr->vertices_end(); ++vi) {
		if (dsc_ptr->is_interface(*vi)) {
			vec2 pos = dsc_ptr->get_pos(*vi);

			//vector<Particle> close_particles = sph_ptr->get_close_particles_to_pos(pos);
			vec2 closest_pos = sph_ptr->get_closest_particle(pos).pos;
			//double closest_distance = 1000000.0;
			/*for (int i = 0; i < sph_ptr->get_no_of_particle(); i++) {
			vec2 r = pos - sph_ptr->get_particle_pos(i);
			if (r.length() < closest_distance) {
			closest_distance = r.length();
			closest_pos = sph_ptr->get_particle_pos(i);
			}
			}*/
			vec2 v = closest_pos - pos;
			double dist = v.length();
			if (pos[0] < 40.0) {
				vec2 vec_to_border = vec2(40.0, pos[1]) - pos;
				velocities[*vi] = vec_to_border * 10.0;
			}
			double r = sph_ptr->get_right_wall();
			if (pos[0] > r + 10.0) {
				vec2 vec_to_border = vec2(r + 10.0, pos[1]) - pos;
				velocities[*vi] = vec_to_border * 10.0;
			}
			if (pos[1] < 102.5) {
				vec2 vec_to_border = vec2(pos[0], 102.5) - pos;
				velocities[*vi] = vec_to_border * 10.0;
				if (pos[0] > r + 10.0) {
					vec2 vec_to_border = vec2(r + 10.0, pos[1]) - pos;
					velocities[*vi] += vec_to_border * 10.0;
				}
				if (pos[0] < 40.0) {
					vec2 vec_to_border = vec2(40.0, pos[1]) - pos;
					velocities[*vi] += vec_to_border * 10.0;
				}
			}
			if (dist > 15.0 && velocities[*vi].length() < 5.0) {
				// Particles far away move towards closest particle
				if (!dsc_ptr->is_outside(*vi)) {
					vec2 vec_to_par = closest_pos - pos;
					velocities[*vi] = vec_to_par;
				}
			}
			else {
				/*Particle tmp_particle = Particle(pos, -1);
				tmp_particle.mass = 0.0;
				vector<Particle> close_particles = sph_ptr->get_close_particles_to_pos(pos);
				float mass_at_point = sph_ptr->calculate_density(tmp_particle, close_particles) * 1000000.0;
				if (mass_at_point > 50.0) {
				//vec2 outward_normal = -sph_ptr->calculate_inward_normal(pos);
				auto w = dsc_ptr->walker(*vi);
				while (true) {
				if (dsc_ptr->is_interface(w.next().vertex())) {
				break;
				}
				w = w.circulate_face_cw();
				}
				vec2 tip = dsc_ptr->get_pos(w.vertex());
				vec2 edge = pos - tip;
				vec2 outward_normal = vec2(edge[1], -edge[0]);
				velocities[*vi] = outward_normal * 10.0;
				}*/
			}
		}
		velocities[*vi] = velocities[*vi];

	}
	if (verts_ids.size() > 0) {
		solve_for_incompressibility(face_ids, verts_ids, velocities, dsc);
		sph_ptr->set_dsc_vert_velocities(velocities);
	}
	for (auto vi = dsc_ptr->vertices_begin(); vi != dsc_ptr->vertices_end(); ++vi) {
		vec2 pos = dsc_ptr->get_pos(*vi);
		dsc_ptr->set_destination(*vi, pos + (velocities[*vi] * 0.015));
	}
	if (tmp == 0)
		m_V0 = get_curent_volume(dsc);
	done(dsc);
	double volume = get_curent_volume(dsc);
	//cout << volume / m_V0 << endl;
	//if( tmp > 0)
	//solve_for_incompressibility_test(velocities, dsc);
	//sph_ptr->set_dsc_vert_velocities(velocities);
	//done(dsc);

	tmp++;
}

void track_particle_function::done(DSC2D::DeformableSimplicialComplex& dsc) {
	double volume = get_curent_volume(dsc);
	sph.set_volume_diff(m_V0);
	if (tmp > 0)
		vol_log.push_back(volume);
	if (tmp == 500) {
		ofstream file;
		file.open("without.txt");
		for (double v : vol_log) {
			file << v << endl;
		}
		file.close();
	}

}

void track_particle_function::solve_for_incompressibility(std::vector<HMesh::FaceID> face_ids,
	std::vector<HMesh::VertexID> verts_ids,
	HMesh::VertexAttributeVector<DSC2D::vec2> &velocities,
	DSC2D::DeformableSimplicialComplex& dsc)
{

	int ROWS = face_ids.size();
	int COLS = verts_ids.size() * 2;
	
	EigVec X = EigVec::Zero(COLS);
	HMesh::VertexAttributeVector<int> num;
	
	int i = 0;
	for (auto v : verts_ids) {
		num[v] = i++;
		X[2 * num[v]] = velocities[v][0];
		X[2 * num[v] + 1] = velocities[v][1];
	}

	EigMat DivOp = EigMat::Zero(ROWS, COLS);

	int R = 0;
	for (auto f : face_ids) {
		
		auto w = dsc.walker(f);

		while (!w.full_circle()) {
			HMesh::VertexID v = w.next().vertex();
			DSC2D::vec2 dir = dsc.get_pos(w.vertex()) - dsc.get_pos(w.opp().vertex());
			DivOp(R, 2 * num[v]) = -dir[1];
			DivOp(R, 2 * num[v] + 1) = dir[0];
			++i;
			w = w.circulate_face_ccw();
		}
		++R;
	}

	EigMat GradOp = DivOp.transpose();
	EigMat LapOp = DivOp * GradOp;
	LDLT<EigMat> solver(LapOp);
	EigVec Xrecon = X - GradOp * solver.solve(DivOp * X);

	for (auto v : verts_ids) {
		DSC2D::vec2 des = DSC2D::vec2(Xrecon[2 * num[v]], Xrecon[2 * num[v] + 1]);
		//cout << des << endl;
		velocities[v] = des;
	}
}

double track_particle_function::get_curent_volume(DSC2D::DeformableSimplicialComplex& dsc) {
	DeformableSimplicialComplex* dsc_ptr = &dsc;
	double V = 0;
	for (auto fit = dsc_ptr->faces_begin(); fit != dsc_ptr->faces_end(); fit++)
	{
		auto fkey = *fit;
		if (dsc_ptr->get_label(fkey) == 1)
		{
			V += dsc_ptr->area(fkey);
		}
	}

	return V;
}

void track_particle_function::solve_for_incompressibility_test(HMesh::VertexAttributeVector<DSC2D::vec2> &velocities, DSC2D::DeformableSimplicialComplex& dsc) {
	DeformableSimplicialComplex* dsc_ptr = &dsc;
	/*
	Volume lost compensation
	*/
	static int t = 0;
	/*if (t < 10)
		m_V0 = get_curent_volume(dsc);*/
	t++;
	double volumeLost = get_curent_volume(dsc) - m_V0;

	//cout << volumeLost << endl;
	// Index the interface veritces
	HMesh::VertexAttributeVector<int> vIdxs(dsc_ptr->get_no_vertices(), -1);
	int num = 0;
	for (auto vit = dsc_ptr->vertices_begin(); vit != dsc_ptr->vertices_end(); vit++)
	{
		auto vkey = *vit;
		if (dsc_ptr->is_interface(vkey) || dsc_ptr->is_crossing(vkey))
		{
			vIdxs[vkey] = num++;
		}
	}

	// Build the equation
	//	using EigMat = SparseMatrix<double>;
	//	using EigVec = VectorXd;

	EigMat M(1 + num, 2 * num);
	for (int i = 0; i < M.rows(); i++) {
		for (int j = 0; j < M.cols(); j++) {
			M(i, j) = 0.0;
		}
	}
	for (auto eIt = dsc_ptr->halfedges_begin(); eIt != dsc_ptr->halfedges_end();
	eIt++)
	{
		auto ekey = *eIt;
		if (dsc_ptr->is_interface(ekey))
		{
			auto hew = dsc_ptr->walker(ekey);
			if (dsc_ptr->get_label(hew.face()) != 1)
			{
				hew = hew.opp();
			}

			auto pt_tip = dsc_ptr->get_pos(hew.vertex());
			auto pt_root = dsc_ptr->get_pos(hew.opp().vertex());
			auto line = pt_tip - pt_root;
			DSC2D::vec2 norm = DSC2D::Util::normalize(DSC2D::vec2(line[1], -line[0]));
			double length = line.length();

			// Matrix M
			int idx1 = vIdxs[hew.vertex()];
			int idx2 = vIdxs[hew.opp().vertex()];
			M.coeffRef(0, 2 * idx1) += norm[0] * length / 2.0;
			M.coeffRef(0, 2 * idx1 + 1) += norm[1] * length / 2.0;

			M.coeffRef(0, 2 * idx2) += norm[0] * length / 2.0;
			M.coeffRef(0, 2 * idx2 + 1) += norm[1] * length / 2.0;
		}
	}

	EigVec B(1 + num);
	for (int i = 0; i < 1 + num; i++)
	{
		B[i] = 0;
	}
	B[0] = -volumeLost;

	// Boundary
	int bound_count = 1;
	for (auto vit = dsc_ptr->vertices_begin(); vit != dsc_ptr->vertices_end(); vit++)
	{
		auto vkey = *vit;
		if (dsc_ptr->is_interface(vkey) || dsc_ptr->is_crossing(vkey))
		{
			DSC2D::vec2 norm;
			
			
			auto n = dsc_ptr->get_normal(vkey);
			norm = DSC2D::vec2(n[1], -n[0]);
			

			int idx = vIdxs[vkey];
			M.coeffRef(bound_count, 2 * idx) = norm[0];
			M.coeffRef(bound_count, 2 * idx + 1) = norm[1];
			bound_count++;
		}
	}

	// Solve
	EigMat M_t = M.transpose();
	EigMat MtM = M_t*M;
	EigVec MtB = M_t*B;

	//cout << " MtB: " << MtB[0] << endl;
	ConjugateGradient<EigMat> cg(MtM);
	EigVec X = cg.solve(MtB);


	// Apply displacement
	for (auto vit = dsc_ptr->vertices_begin(); vit != dsc_ptr->vertices_end(); vit++)
	{
		auto vkey = *vit;

		int idx = vIdxs[vkey];
		if (idx != -1)
		{
			DSC2D::vec2 dis(X[2 * idx], X[2 * idx + 1]);
			velocities[*vit] = dis;
			DSC2D::vec2 v = dsc_ptr->get_pos(vkey) + dis;
			v = v - dsc.get_pos(vkey);

			if (v.length() > 5.0) {
				v.normalize();
				v = v * 5.0;
			}
			
			dsc_ptr->set_destination(vkey, dsc_ptr->get_pos(vkey) + v);
		}
	}

	//dsc_ptr->deform();

	/*
	Project back to boundary
	*/
	/*for (auto vit = dsc_ptr->vertices_begin(); vit != dsc_ptr->vertices_end(); vit++)
	{
	auto vkey = *vit;
	if (dsc_ptr->is_interface(vkey) || dsc_ptr->is_crossing(vkey))
	{


	if (dsc.get_pos(vkey)[1] <= particle_system.get_border_floor() && dsc.get_pos(vkey)[0] >= particle_system.get_border_right()) {
	dsc.set_destination(vkey, DSC2D::vec2(particle_system.get_border_right() + 20, particle_system.get_border_floor() - 20));
	}
	else if (dsc.get_pos(vkey)[1] <= particle_system.get_border_floor() && dsc.get_pos(vkey)[0] <= particle_system.get_border_left()) {
	dsc.set_destination(vkey, DSC2D::vec2(particle_system.get_border_left() - 20, particle_system.get_border_floor() - 20));
	}
	else if (dsc.get_pos(vkey)[0] >= particle_system.get_border_right()) {
	dsc.set_destination(vkey, DSC2D::vec2(particle_system.get_border_right() + 20, dsc.get_pos(vkey)[1]));
	}
	else if (dsc.get_pos(vkey)[1] <= particle_system.get_border_floor()) {
	dsc.set_destination(vkey, DSC2D::vec2(dsc.get_pos(vkey)[0], particle_system.get_border_floor() - 20));
	}
	else if (dsc.get_pos(vkey)[0] <= particle_system.get_border_left()) {
	dsc.set_destination(vkey, DSC2D::vec2(particle_system.get_border_left() - 20, dsc.get_pos(vkey)[1]));
	}
	}
	}
	dsc.deform();*/
}