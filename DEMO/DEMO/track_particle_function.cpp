#include "track_particle_function.h"
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "object_generator.h"

using namespace Eigen;

using EigMat = MatrixXd;
using EigVec = VectorXd;

void track_particle_function::deform(DSC2D::DeformableSimplicialComplex& dsc)
{


	dsc_ptr = &dsc;

	/*
	Displace DSC
	*/
	// Gravity
	HMesh::VertexAttributeVector<DSC2D::vec2> vels(dsc_ptr->get_no_vertices(), DSC2D::vec2(0.0));
	for (auto vi = dsc.vertices_begin(); vi != dsc.vertices_end(); ++vi)
	{
		if (dsc.is_interface(*vi) || dsc.is_crossing(*vi))
		{
			vels[*vi] += DSC2D::vec2(0,-3);
		}
	}
	// curvature
	for (auto eIt = dsc_ptr->halfedges_begin(); eIt != dsc_ptr->halfedges_end();
		eIt++)
	{
		auto hekey = *eIt;

		if (dsc_ptr->is_interface(hekey))
		{
			auto hew = dsc_ptr->walker(hekey);
			auto l = dsc_ptr->get_pos(hew.opp().vertex()) - dsc_ptr->get_pos(hew.vertex());
			l.normalize();
			vels[hew.vertex()] += l*1.1;
		}
	}


	for (auto vi = dsc.vertices_begin(); vi != dsc.vertices_end(); ++vi)
	{
		if (dsc.is_movable(*vi))
		{
			if (dsc.is_interface(*vi) || dsc.is_crossing(*vi))
			{
				auto ff = vels[*vi];
				dsc.set_destination(*vi, dsc.get_pos(*vi) + vels[*vi]);
			}
		}
	}

	dsc.deform();


	/*
	Volume lost compensation
	*/
	double volumeLost = get_curent_volume() - m_V0;

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
			if (!is_on_boundary(dsc_ptr->get_pos(vkey), norm))
			{
				auto n = dsc_ptr->get_normal(vkey);
				norm = DSC2D::vec2(n[1], -n[0]);
			}

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
			dsc_ptr->set_destination(vkey, dsc_ptr->get_pos(vkey) + dis);
		}
	}

	dsc_ptr->deform();

	/*
	Project back to boundary
	*/
	for (auto vit = dsc_ptr->vertices_begin(); vit != dsc_ptr->vertices_end(); vit++)
	{
		auto vkey = *vit;
		if (dsc_ptr->is_interface(vkey) || dsc_ptr->is_crossing(vkey))
		{
			DSC2D::vec2 p_v;
			if (is_outside(dsc_ptr->get_pos(vkey), p_v))
			{
				dsc_ptr->set_destination(vkey, p_v);
			}
		}
	}
	dsc.deform();
}

void  track_particle_function::init()
{
	using namespace DSC2D;
	static bool inited = false;
	if (!inited)
	{
		inited = true;

		// Init inside the DSC, half size
		auto corners = dsc_ptr->get_design_domain()->get_corners();
		auto center = dsc_ptr->get_center();
		vec2 diag = corners[2] - corners[0];
		diag = diag / 4.0;

		auto ld_ = center - diag;
		auto ru_ = center + diag;

		DSC2D::vec2 gap(10, 10);

		DSC2D::ObjectGenerator::create_square(*dsc_ptr, ld_ - gap, ru_ - ld_ + gap, 1);


		m_V0 = get_curent_volume();

		// Boundary
		auto pts = dsc_ptr->get_design_domain()->get_corners();

		center_bound = (pts[0] + pts[2]) / 2;
		r_bound = (pts[1] - pts[0]).length() * 0.5 * 0.8;
	}
}


bool track_particle_function::is_on_boundary(DSC2D::vec2 pt, DSC2D::vec2 &norm)
{
	double l = (pt - center_bound).length();
	if (l > r_bound - 1e-2)
	{
		norm = DSC2D::Util::normalize(pt - center_bound);
		return true;
	}

	return false;
}

bool track_particle_function::is_outside(DSC2D::vec2 pt, DSC2D::vec2 &projectionPoint)
{
	double l = (pt - center_bound).length();
	if (l > r_bound)
	{
		projectionPoint = center_bound + DSC2D::Util::normalize(pt - center_bound)*r_bound;
		return true;
	}

	return false;
}

double track_particle_function::get_curent_volume(){
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
void track_particle_function::deform_old(DSC2D::DeformableSimplicialComplex& dsc){
	r =  15;
	alpha = r * 0.3;
	no_segments = 3;
	
	auto init_time = std::chrono::system_clock::now();
	std::vector<DSC2D::vec2> added_forces;
	dsc.ready_destination();
	// loop over all triangles
	float volume = 0.0;
	float mass = 0.0;
	// Q? Volume of entire DSC domain?
	for (auto fi = dsc.faces_begin(); fi != dsc.faces_end(); ++fi){
		volume += dsc.area(*fi);
	}
	for (int i = 0; i < particle_system.get_no_particles(); i++){
		mass += particle_system.get_particle_mass(i);
	}
	mean_density = mass / volume;

	std::vector<HMesh::FaceID> face_ids;
	HMesh::VertexAttributeVector<int> is_vert_inside;

	HMesh::VertexAttributeVector<DSC2D::vec2> velocities;

	for (auto vi = dsc.vertices_begin(); vi != dsc.vertices_end(); vi++) {
			velocities[*vi] = DSC2D::vec2(0.0f);
	}
	
	for (auto fi = dsc.faces_begin(); fi != dsc.faces_end(); ++fi){
		auto vis = dsc.get_verts(*fi);
		DSC2D::vec2 force_v0 = DSC2D::vec2(0);
		DSC2D::vec2 force_v1 = DSC2D::vec2(0);
		if (dsc.get_label(*fi) == 1) {
			face_ids.push_back(*fi);
			for (int i = 0; i < 3; i++)
				is_vert_inside[vis[i]] = 1; // Can use dsc.is_outside() instead
		}


		// Check every edge (vertex pair) if they are on the interface and then apply force

		/*
			Should make a loop on all interface edges instead
			--IMPORTANT--
		*/
		std::vector<int> labels = dsc.get_interface_labels(vis[0]);
		std::vector<int> labels2 = dsc.get_interface_labels(vis[1]);
		if (dsc.is_interface(vis[0]) && dsc.is_interface(vis[1]) && dsc.get_label(*fi) == 1){ // edge v0 to v1
			calculate_forces(dsc.get_pos(vis[0]), dsc.get_pos(vis[1]), &force_v1, &force_v0, dsc, vis[0], vis[1]);
			/*if (dsc.get_pos(vis[0])[1] > 40.0 && dsc.get_pos(vis[1])[1] > 40.0) {
				dsc.add_destination(vis[0], force_v0*alpha);
				dsc.add_destination(vis[1], force_v1*alpha);
			}*/

			//velocities[vis[0]] += force_v0*alpha;
			//velocities[vis[1]] += force_v1*alpha;
			if (dsc.get_pos(vis[1])[1] > 50.0f)
				velocities[vis[1]][1] = -1.0f;

			if (dsc.get_pos(vis[0])[1] > 50.0f)
				velocities[vis[0]][1] = -1.0f;
			//dsc.add_destination(vis[0], force_v0*alpha);
			//dsc.add_destination(vis[1], force_v1*alpha);
			force_v0 = DSC2D::vec2(0);
			force_v1 = DSC2D::vec2(0);
		}
		if (dsc.is_interface(vis[1]) && dsc.is_interface(vis[2]) && dsc.get_label(*fi) == 1){ // edge v1 to v2
			calculate_forces(dsc.get_pos(vis[1]), dsc.get_pos(vis[2]), &force_v1, &force_v0, dsc, vis[1], vis[2]);
			/*if (dsc.get_pos(vis[1])[1] > 40.0 && dsc.get_pos(vis[2])[1] > 40.0) {
				dsc.add_destination(vis[1], force_v0*alpha);
				dsc.add_destination(vis[2], force_v1*alpha);
			}*/
			
			//velocities[vis[1]] += force_v0*alpha;
			//velocities[vis[2]] += force_v1*alpha;
			if (dsc.get_pos(vis[1])[1] > 50.0f)
				velocities[vis[1]][1] = -1.0f;

			if (dsc.get_pos(vis[2])[1] > 50.0f)
				velocities[vis[2]][1] = -1.0f;
			//dsc.add_destination(vis[1], force_v0*alpha);
			//dsc.add_destination(vis[2], force_v1*alpha);
			force_v0 = DSC2D::vec2(0);
			force_v1 = DSC2D::vec2(0);
		}
		if (dsc.is_interface(vis[2]) && dsc.is_interface(vis[0]) && dsc.get_label(*fi) == 1){ // edge v2 to v0
			calculate_forces(dsc.get_pos(vis[2]), dsc.get_pos(vis[0]), &force_v1, &force_v0, dsc, vis[2], vis[0]);
			/*if (dsc.get_pos(vis[2])[1] > 40.0 && dsc.get_pos(vis[0])[1] > 40.0) {
				dsc.add_destination(vis[2], force_v0*alpha);
				dsc.add_destination(vis[0], force_v1*alpha);
			}*/

			//velocities[vis[2]] += force_v0*alpha;
			//velocities[vis[0]] += force_v1*alpha;
			if (dsc.get_pos(vis[2])[1] > 50.0f)
				velocities[vis[2]][1] = -1.0f;

			if (dsc.get_pos(vis[0])[1] > 50.0f)
				velocities[vis[0]][1] = -1.0f;
			//dsc.add_destination(vis[2], force_v0*alpha);
			//dsc.add_destination(vis[0], force_v1*alpha);
			force_v0 = DSC2D::vec2(0);
			force_v1 = DSC2D::vec2(0);
		}
	}

	std::vector<HMesh::VertexID> verts_ids;
	for (auto vi = dsc.vertices_begin(); vi != dsc.vertices_end(); vi++) {
		if (is_vert_inside[*vi] == 1) {
			verts_ids.push_back(*vi);
		}
	}
	float number_of_verts = verts_ids.size();
	
// 	// Solve volume preservation
// 	if (solver_delay > 0) {
// 		solver_delay -= 1;
// 	}
// 	else if (solver_delay == 0) {
// 		solver_delay -= 1;
// 		max_diff = 0;
// 		for (auto fi = dsc.faces_begin(); fi != dsc.faces_end(); ++fi) {
// 			
// 			if (dsc.get_label(*fi) == 1) {
// 				start_area += dsc.area(*fi);
// 			}
// 
// 		}
// 	}
// 	else 
	{
		//printf("number of vertices: %i", velocities.size());
		
		solve_for_incompressibility(face_ids, verts_ids, velocities, dsc);

	}
	
	for (auto vi = dsc.vertices_begin(); vi != dsc.vertices_end(); vi++) {
		dsc.set_destination_list_element(*vi, (dsc.get_pos(*vi) + velocities[*vi]));
		//dsc.set_destination(*vi, dsc.get_pos(*vi) + velocities[*vi]);
		/*
		if (!dsc.is_interface(*vi) && is_vert_inside[*vi] == 1) {
			//printf("* X: %f Y: %f *", des[0], des[1]);
			dsc.set_destination_list_element(*vi, (dsc.get_pos(*vi) + velocities[*vi]));
			//dsc.move_vertex(*vi);
		}*/
	}

	//dsc.sum_destination();
	
	update_compute_time(init_time);
	init_time = std::chrono::system_clock::now();

	dsc.deform();
	float sum_area = 0.0f;
	for (auto fi = dsc.faces_begin(); fi != dsc.faces_end(); ++fi) {
		if (dsc.get_label(*fi) == 1) {
			sum_area += dsc.area(*fi);
		}

	}
	float percentage_changed = (sum_area / last_area) - 1.0;
	if (max_diff < abs(percentage_changed)) {
		max_diff = abs(percentage_changed);
	}
	//printf("AREA: %f", sum_area);
	printf("AREA DIFF: %f, MAX DIFF: %f, START_NOW_DIFF: %f \n", sum_area - last_area, max_diff, start_area - sum_area);
	last_area = sum_area;
	update_deform_time(init_time);
	// Makes sure the border interface gets clamped at the border of the container.
	/*for (auto fi = dsc.faces_begin(); fi != dsc.faces_end(); ++fi){
		auto vis = dsc.get_verts(*fi);
		if (dsc.is_interface(vis[0])){
			for (int i = 0; i < vis.size(); i++){
				if (dsc.get_pos(vis[i])[1] <= particle_system.get_border_floor() && dsc.get_pos(vis[i])[0] >= particle_system.get_border_right()){
					dsc.set_destination(vis[i], DSC2D::vec2(particle_system.get_border_right() + 15, particle_system.get_border_floor() - 15));
				}
				else if (dsc.get_pos(vis[i])[1] <= particle_system.get_border_floor() && dsc.get_pos(vis[i])[0] <= particle_system.get_border_left()){
					dsc.set_destination(vis[i], DSC2D::vec2(particle_system.get_border_left() - 15, particle_system.get_border_floor() - 15));
				}
				else if (dsc.get_pos(vis[i])[0] >= particle_system.get_border_right()){
					dsc.set_destination(vis[i], DSC2D::vec2(particle_system.get_border_right() + 15, dsc.get_pos(vis[i])[1]));
				}
				else if (dsc.get_pos(vis[i])[1] <= particle_system.get_border_floor()){
					dsc.set_destination(vis[i], DSC2D::vec2(dsc.get_pos(vis[i])[0], particle_system.get_border_floor() - 15));
				}
				else if (dsc.get_pos(vis[i])[0] <= particle_system.get_border_left()){
					dsc.set_destination(vis[i], DSC2D::vec2(particle_system.get_border_left() - 15, dsc.get_pos(vis[i])[1]));
				}
			}

		}
	}*/
	//dsc.deform();
}
// Calculate the forces and distribute them to the vertices of the edge
void track_particle_function::calculate_forces(DSC2D::vec2 pos1, DSC2D::vec2 pos2, DSC2D::vec2 *force1, DSC2D::vec2 *force2, DSC2D::DeformableSimplicialComplex& dsc, HMesh::VertexID vert1, HMesh::VertexID vert2){
	DSC2D::vec2 v0v1 = pos2 - pos1; // the edge as a vector
	float edge_length = CGLA::length(v0v1); // lenght of the edge
	float segment_length = edge_length / (float)no_segments; // length of one segment
	for (int i = 0; i < 3; i++){ // loop over all edge-segments
		// Calculate the center of the edge segment
		DSC2D::vec2 segment_center = pos1 + (CGLA::normalize(v0v1) * (((segment_length)* (1 + i)) - (segment_length*0.5)));
		// Loop over all particles and find the closes one
		DSC2D::vec2 closest = DSC2D::vec2(100000);
		DSC2D::vec2 cen_par = closest - segment_center;
		float closest_distance = 100000;
		for (int i = 0; i < particle_system.get_no_particles(); i++){
			if (particle_system.is_particle_inside(i) == true) {
				cen_par = particle_system.get_particle_pos(i) - segment_center;
				float distance = CGLA::length(cen_par);
				if (distance < closest_distance){
					bool is_inside = true;
					
					if (is_inside) {
						//printf("h");
						closest = particle_system.get_particle_pos(i);
						closest_distance = distance;
					}
					else {
						closest = particle_system.get_particle_pos(i);
					}
				}
			}
		}
		
		//printf("Dis: [%f, %f]", closest_distance, closest_distance);
		// Calculate outward normal
		DSC2D::vec2 par_cen = segment_center - closest;
		DSC2D::vec2 outward_normal = (segment_center + par_cen) - segment_center;
		outward_normal = CGLA::normalize(outward_normal);
		// Calculate the forces on the two vertices
		DSC2D::vec2 force_particle_v0 = DSC2D::vec2(0);
		DSC2D::vec2 force_particle_v1 = DSC2D::vec2(0);
		if (closest_distance > r * 1.6){
			force_particle_v0 = -outward_normal * 0.5f;
			force_particle_v1 = -outward_normal * 0.5f;
		}
		else if (closest_distance < r){
			
			force_particle_v0 = (2.0 * particle_system.fluid_density(pos1, 50.0) - mean_density) * outward_normal * 300.0;
			force_particle_v1 = (2.0 * particle_system.fluid_density(pos2, 50.0) - mean_density) * outward_normal * 300.0;
			
		}
		DSC2D::vec2 vec_l1 = segment_center - pos1;
		float l1 = CGLA::length(vec_l1);
		DSC2D::vec2 vec_l2 = segment_center - pos2;
		float l2 = CGLA::length(vec_l1);
		*force1 += force_particle_v0 * (l2 / (l1 + l2));
		*force2 += force_particle_v1 * (l1 / (l1 + l2));
	}
}

void track_particle_function::solve_for_incompressibility(std::vector<HMesh::FaceID> face_ids, 
	std::vector<HMesh::VertexID> verts_ids, 
	HMesh::VertexAttributeVector<DSC2D::vec2> &velocities, 
	DSC2D::DeformableSimplicialComplex& dsc) 
{
	/*int ROWS = 0;
	int COLS = 0;
	/*for (auto f : face_ids) {
		HMesh::Walker w = dsc.walker(dsc.get_verts(f)[0]);
		bool is_on_interface = false;
		while (!w.full_circle()) {
			if (dsc.get_label(w.vertex()) == 1) {
				COLS += 1;
				is_on_interface = true;
			}
			w = w.circulate_face_ccw();
		}
		if (is_on_interface) {
			ROWS += 1;
		}
	}*/
	
	int ROWS = face_ids.size();
	int COLS = verts_ids.size() * 2;
	//COLS *= 2;
	EigVec X = EigVec::Zero(COLS);
	HMesh::VertexAttributeVector<int> num;
	//printf("ROWS: %i, COLS: %i", ROWS, COLS);
	int i = 0;
	for (auto v : verts_ids) {
		num[v] = i++;
		X[2 * num[v]] = velocities[v][0];
		X[2 * num[v] + 1] = velocities[v][1];
	}

	EigMat DivOp = EigMat::Zero(ROWS, COLS);

	int R = 0;
	for (auto f : face_ids) {
		//HMesh::Walker w = dsc.walker(dsc.get_verts(f)[0]);
		// Q? Should be walker(f). Otherwise, it may be random
		// IMPORTANT
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
		//dsc.set_destination_list_element(v, des);
		velocities[v] = des;
		/*if (!dsc.is_interface(v)) {
			//printf("* X: %f Y: %f *", des[0], des[1]);
			dsc.set_destination_list_element(v, (dsc.get_pos(v) + des));
			dsc.move_vertex(v);
		}*/
	}
}