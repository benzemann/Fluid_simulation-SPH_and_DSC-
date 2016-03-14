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

void track_particle_function::deform(DSC2D::DeformableSimplicialComplex& dsc){
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
				if (pos[0] > r+10.0) {
					vec2 vec_to_border = vec2(r+10.0, pos[1]) - pos;
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
	if( tmp > 0)
		solve_for_incompressibility_test(velocities, dsc);
	//sph_ptr->set_dsc_vert_velocities(velocities);
	//done(dsc);
	
	vol_log.push_back(volume);
	if (tmp == 500) {
		ofstream file;
		file.open("par_vel_v0.txt");
		for (double v : vol_log) {
			file << v << endl;
		}
		file.close();
	}
	tmp++;
}

void track_particle_function::done(DSC2D::DeformableSimplicialComplex& dsc) {
	dsc.deform();
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
	if (t < 10)
		m_V0 = get_curent_volume(dsc);
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
			dsc_ptr->set_destination(vkey, dsc_ptr->get_pos(vkey) + dis);
		}
	}

	dsc_ptr->deform();

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