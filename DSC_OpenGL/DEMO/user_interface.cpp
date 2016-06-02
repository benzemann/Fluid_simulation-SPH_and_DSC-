//
//  Deformabel Simplicial Complex (DSC) method
//  Copyright (C) 2013  Technical University of Denmark
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  See licence.txt for a copy of the GNU General Public License.

#include "user_interface.h"
#include "rotate_function.h"
#include "average_function.h"
#include "normal_function.h"
#include "draw_helper.h"
#include "track_particles_function.h"

#include <iostream>
#include <iomanip>
#include <ctime>
#include <chrono>

using namespace DSC;

void display_(){
    UI::get_instance()->display();
}

void keyboard_(unsigned char key, int x, int y){
    UI::get_instance()->keyboard(key, x, y);
}

void reshape_(int width, int height){
    UI::get_instance()->reshape(width, height);
}

void visible_(int v){
    UI::get_instance()->visible(v);
}

void animate_(){
    UI::get_instance()->animate();
}

UI* UI::instance = NULL;


void UI::setup_light()
{
	vec3 center = vec3(2, 2, 2);
	vec3 eye = center + vec3(gl_dis_max*3.0*cos(angle)*cos(angle2),
		gl_dis_max*3.0*cos(angle)*sin(angle2),
		gl_dis_max*3.0*sin(angle));
	eye.normalize();
	eye = eye * 30.0;
	GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat mat_shininess[] = { 50.0 };
	lightPos = vec3(-(GLfloat)eye[0], -(GLfloat)eye[1], -(GLfloat)eye[2]);
	GLfloat light_position[] = { -(GLfloat)eye[0], -(GLfloat)eye[1], -(GLfloat)eye[2], 0.0 };
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glShadeModel(GL_SMOOTH);

	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_DEPTH_TEST);
}
inline vec3 min_vec(vec3 v1, vec3 v2)
{
	vec3 out;
	out[0] = std::min(v1[0], v2[0]);
	out[1] = std::min(v1[1], v2[1]);
	out[2] = std::min(v1[2], v2[2]);
	return out;
}

inline vec3 max_vec(vec3 v1, vec3 v2)
{
	vec3 out;
	out[0] = std::max(v1[0], v2[0]);
	out[1] = std::max(v1[1], v2[1]);
	out[2] = std::max(v1[2], v2[2]);
	return out;
}

UI::UI(int &argc, char** argv)
{
	instance = this;

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_STENCIL | GLUT_MULTISAMPLE);

	glutCreateWindow("Shadowy Leapin' Lizards");

	glutDisplayFunc(display_);
	glutKeyboardFunc(keyboard_);
	glutSetKeyRepeat(GLUT_KEY_REPEAT_ON);
	glutVisibilityFunc(visible_);
	glutReshapeFunc(reshape_);


	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	glEnable(GL_DEPTH_TEST);
	glLineWidth(1.0);

	setup_light();

	glutReshapeWindow(WIN_SIZE_X, WIN_SIZE_Y);
	check_gl_error();
	
	load_model("C:/Users/Jeppe/Desktop/DSC_OpenGL/DSC_OpenGL/DSC.visualstudio/cube_simple.dsc", 2.5);
	
	real velocity = 5.;
	real accuracy = 0.25;
	vel_fun = std::unique_ptr<VelocityFunc<>>(new VelocityFunc<>(velocity, accuracy, 500));
	start("");

	//
	vec3 minc(INFINITY);
	vec3 maxc(-INFINITY);
	for (auto vit = dsc->nodes_begin(); vit != dsc->nodes_end(); vit++)
	{
		vec3 pt = dsc->get_pos(vit.key());

		minc = min_vec(minc, pt);
		maxc = max_vec(maxc, pt);
	}
	_obj_dim = maxc - minc;
	gl_dis_max = std::max(std::max(_obj_dim[0], _obj_dim[1]), _obj_dim[2])*2;
	for (int i = 0; i < 100; i++) {
		real discretization = std::max(dsc->get_avg_edge_length() - 0.5, 1.);
		dsc->set_avg_edge_length(discretization);
		update_title();
	}


	dsc->scale(vec3(8.0, 8.0, 8.0));

	int i = 0;
	vector<is_mesh::TetrahedronKey> keys;
	for (auto te = dsc->tetrahedra_begin(); te != dsc->tetrahedra_end(); te++) {
		keys.push_back(te.key());
	}
	for each (is_mesh::TetrahedronKey key in keys)
	{
		dsc->split(key);
	}
	cout << "SPLIT" << endl;
	keys.clear();
	for (auto te = dsc->tetrahedra_begin(); te != dsc->tetrahedra_end(); te++) {
		keys.push_back(te.key());
	}
	for each (is_mesh::TetrahedronKey key in keys)
	{
		dsc->split(key);
	}
	cout << "SPLIT" << endl;
	keys.clear();
	for (auto te = dsc->tetrahedra_begin(); te != dsc->tetrahedra_end(); te++) {
		keys.push_back(te.key());
	}
	for each (is_mesh::TetrahedronKey key in keys)
	{
		dsc->split(key);
	}
	cout << "SPLIT" << endl;
	keys.clear();
	for (auto te = dsc->tetrahedra_begin(); te != dsc->tetrahedra_end(); te++) {
		vel_fun->set_label(*dsc, te.key(), 0);
		/*if (i < 75000) {
			
			is_mesh::SimplexSet<is_mesh::NodeKey> set = dsc->get_nodes(te.key());
			bool m = false;
			for (auto n = set.begin(); n != set.end(); n++) {
				//std::cout << dsc->get_pos(*n) << std::endl;
				if (dsc->get_pos(*n)[1] < 0.8 && dsc->get_pos(*n)[1] > -0.8 && dsc->get_pos(*n)[0] < 0.8 && dsc->get_pos(*n)[0] > -0.8 && dsc->get_pos(*n)[2] < 0.8 && dsc->get_pos(*n)[2] > -0.8) {
					m = true;
				}
			}
			if (m) {
				vel_fun->set_label(*dsc, te.key(), 1);
			}
			else {
				vel_fun->set_label(*dsc, te.key(), 0);
			}
				

		}
		else {

			
		}*/

		i++;
	}
	std::cout << "Number of tets: " << i << std::endl;
	sph = new SPH();
	sph->init();
}

void UI::load_model(const std::string& file_name, real discretization)
{
	std::cout << file_name << std::endl;
    dsc = nullptr;
    std::vector<vec3> points;
    std::vector<int>  tets;
    std::vector<int>  tet_labels;
	is_mesh::import_tet_mesh(file_name, points, tets, tet_labels);
	
    dsc = std::unique_ptr<DeformableSimplicialComplex<>>(new DeformableSimplicialComplex<>(points, tets, tet_labels));
    
    std::cout << "Loading done" << std::endl << std::endl;
}

void UI::update_title()
{
    std::ostringstream oss;
    oss << "3D DSC\t" << vel_fun->get_name() << ", Time step " << vel_fun->get_time_step();
    oss << " (Nu = " << vel_fun->get_velocity() << ", Delta = " << dsc->get_avg_edge_length() << ", Alpha = " << vel_fun->get_accuracy() << ")";
    std::string str(oss.str());
    glutSetWindowTitle(str.c_str());
}

void UI::update_gl()
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective( /* field of view in degree */ 40.0,
		/* aspect ratio */ 1.0,
		/* Z near */ 1.0, /* Z far */ gl_dis_max*150.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	vec3 center = vec3(2, 2, 2);
	vec3 eye = center + vec3(gl_dis_max*3.0*cos(angle)*cos(angle2),
		gl_dis_max*3.0*cos(angle)*sin(angle2),
		gl_dis_max*3.0*sin(angle));
	vec3 head = vec3(-sin(angle)*cos(angle2),
		-sin(angle)*sin(angle2),
		cos(angle));
	eye.normalize();
	eye = eye * 30.0;
	gluLookAt(eye[0], eye[1], eye[2], /* eye is at (0,8,60) */
		center[0], center[1], center[2],      /* center is at (0,8,0) */
		head[0], head[1], head[2]);      /* up is in postivie Y direction */
	
	int size = std::min(WIN_SIZE_Y, WIN_SIZE_X);
	glViewport((WIN_SIZE_X - size) / 2.0, (WIN_SIZE_Y - size) / 2.0, size, size);

	glClearColor(0.1, 0.1, 0.1, 1.0);
}

bool screenshot(char *fileName) {

	int Xres = 640*2;
	int Yres = 360*2;

	static unsigned char header[54] = {
		0x42, 0x4D, 0x36, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x36, 0x00, 0x00, 0x00, 0x28, 0x00,
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x18, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x00, 0x03, 0x00, 0xC4, 0x0E, 0x00, 0x00, 0xC4, 0x0E, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

	unsigned char *pixels = (unsigned char *)malloc(Xres * Yres * 3);
	((unsigned __int16 *)header)[9] = Xres;
	((unsigned __int16 *)header)[11] = Yres;

	glReadPixels(0, 0, Xres, Yres, GL_RGB, GL_UNSIGNED_BYTE, pixels);

	unsigned char temp;
	for (unsigned int i = 0; i < Xres * Yres * 3; i += 3) {
		temp = pixels[i];
		pixels[i] = pixels[i + 2];
		pixels[i + 2] = temp;
	}

	HANDLE FileHandle;
	unsigned long Size;

	if (fileName == NULL) {
		char file[256];
		int i = 0;
		do {
			sprintf(file, "Screenshot%d.bmp", i);
			FileHandle = CreateFile(file, GENERIC_WRITE, 0, NULL, CREATE_NEW, FILE_ATTRIBUTE_NORMAL, NULL);
			i++;
		} while (FileHandle == INVALID_HANDLE_VALUE);
	}
	else {
		FileHandle = CreateFile(fileName, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
		if (FileHandle == INVALID_HANDLE_VALUE)	return false;
	}

	WriteFile(FileHandle, header, sizeof(header), &Size, NULL);
	WriteFile(FileHandle, pixels, Xres * Yres * 3, &Size, NULL);

	CloseHandle(FileHandle);

	free(pixels);
	return true;
}
int c_tmp = 0;
void UI::display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	update_gl();
	setup_light();

	static double total_time = 100;
	static auto init_time = std::chrono::system_clock::now();
	std::chrono::duration<real> t = std::chrono::system_clock::now() - init_time;
	total_time += t.count();
	init_time = std::chrono::system_clock::now();
	
	draw_helper::dsc_draw_domain(*dsc);
	//draw_helper::dsc_draw_edge(*dsc);
	draw_helper::dsc_draw_interface(*dsc);


	if (CONTINUOUS) {
		//cout << vel_fun->get_time_step() % 10 << endl;

		//if (vel_fun->get_time_step() % 2 == 0 && vel_fun->get_time_step() < 160000) {
		
		string name = "C:/Users/Jeppe/Desktop/Screenshots/" + std::to_string(vel_fun->get_time_step());
		name = name + ".BMP";
		screenshot((char*)name.c_str());
		//}
	}
	else {
		//sph->draw_particles(lightPos);
	}
	sph->draw_particles(lightPos);

	double d_t = sph->CFL_correction(0.008);
	//cout << d_t << endl;
	sph->update(d_t);
	
	sph->update_velocity(d_t);

	sph->correct_density_error(d_t);

	sph->update_position(d_t);

	//sph->draw_particles(lightPos);
	//sph->draw_collision_boxes(lightPos);

	//sph->update_iso_surface();


	string name = "C:/Users/Jeppe/Desktop/Screenshots/" + std::to_string(c_tmp);
	name = name + ".BMP";
	screenshot((char*)name.c_str());
	c_tmp++;
	glutSwapBuffers();
}

void UI::reshape(int width, int height)
{
	WIN_SIZE_X = width;
	WIN_SIZE_Y = height;

	update_gl();

	glutReshapeWindow(WIN_SIZE_X, WIN_SIZE_Y);
}



void UI::animate()
{

	if (CONTINUOUS)
	{
		std::cout << "\n***************TIME STEP " << vel_fun->get_time_step() + 1 << " START*************\n" << std::endl;

		identify_fluid(220);
		vel_fun->take_time_step(*dsc);
	}
    glutPostRedisplay();
}

void UI::keyboard(unsigned char key, int x, int y) {
    switch(key) {
        case '\033':
            stop();
            exit(0);
            break;
        case '0':
            stop();
            QUIT_ON_COMPLETION = false;
            RECORD = false;
            vel_fun = std::unique_ptr<VelocityFunc<>>(new VelocityFunc<>(vel_fun->get_velocity(), vel_fun->get_accuracy(), 500));
            start("");
            break;
        case '1':
            stop();
            QUIT_ON_COMPLETION = true;
            RECORD = true;
            vel_fun = std::unique_ptr<VelocityFunc<>>(new track_particles_function(vel_fun->get_velocity(), vel_fun->get_accuracy(), sph));
            start("rotate");
            break;
        case '2':
            stop();
            QUIT_ON_COMPLETION = true;
            RECORD = true;
            vel_fun = std::unique_ptr<VelocityFunc<>>(new AverageFunc(vel_fun->get_velocity(), vel_fun->get_accuracy()));
            start("smooth");
            break;
        case '3':
            stop();
            QUIT_ON_COMPLETION = true;
            RECORD = true;
            vel_fun = std::unique_ptr<VelocityFunc<>>(new NormalFunc(vel_fun->get_velocity(), vel_fun->get_accuracy()));
            start("expand");
            break;
		case '4':
			sph->reset();
			break;
        case ' ':
            if(!CONTINUOUS)
            {
                std::cout << "MOTION STARTED" << std::endl;
                if(RECORD && basic_log)
                {
                    painter->set_view_position(camera_pos);
                    painter->save_painting(basic_log->get_path(), vel_fun->get_time_step());
                }
            }
            else {
                std::cout << "MOTION PAUSED" << std::endl;
            }
            CONTINUOUS = !CONTINUOUS;
            break;
        case 'm':
            std::cout << "MOVE" << std::endl;
            vel_fun->take_time_step(*dsc);
            painter->update(*dsc);
            break;
        case 'r':
            std::cout << "RELOAD MODEL" << std::endl;
            load_model(model_file_name, dsc->get_avg_edge_length());
            break;
        case 't':
            std::cout << "TEST VELOCITY FUNCTION" << std::endl;
            vel_fun->test(*dsc);
            painter->update(*dsc);
            break;
        case '\t':
            painter->switch_display_type();
            painter->update(*dsc);
            break;
        case 's':
            std::cout << "TAKING SCREEN SHOT" << std::endl;
            //painter->set_view_position(camera_pos);
            //painter->save_painting("LOG");
			screenshot("Test.bmp");
            break;
        case 'e':
        {
            std::cout << "EXPORTING MESH" << std::endl;
            std::string filename("data/mesh.dsc");
            std::vector<vec3> points;
            std::vector<int> tets;
            std::vector<int> tet_labels;
            dsc->extract_tet_mesh(points, tets, tet_labels);
            is_mesh::export_tet_mesh(filename, points, tets, tet_labels);
        }
            break;
        case 'i':
        {
            std::cout << "EXPORTING SURFACE MESH" << std::endl;
            std::string filename("data/mesh.obj");
            std::vector<vec3> points;
            std::vector<int> faces;
            dsc->extract_surface_mesh(points, faces);
            is_mesh::export_surface_mesh(filename, points, faces);
        }
            break;
        case '+':
        {
            real velocity = std::min(vel_fun->get_velocity() + 1., 100.);
            vel_fun->set_velocity(velocity);
            update_title();
        }
            break;
        case '-':
        {
            real velocity = std::max(vel_fun->get_velocity() - 1., 0.);
            vel_fun->set_velocity(velocity);
            update_title();
        }
            break;
        case '.':
        {
            real discretization = std::min(dsc->get_avg_edge_length() + 0.5, 100.);
            dsc->set_avg_edge_length(discretization);
            update_title();
        }
            break;
        case ',':
        {
            real discretization = std::max(dsc->get_avg_edge_length() - 0.5, 1.);
            dsc->set_avg_edge_length(discretization);
            update_title();
        }
            break;
        case '<':
        {
            real accuracy = std::min(vel_fun->get_accuracy() + 1., 100.);
            vel_fun->set_accuracy(accuracy);
            update_title();
        }
            break;
        case '>':
        {
            real accuracy = std::max(vel_fun->get_accuracy() - 1., 1.);
            vel_fun->set_accuracy(accuracy);
            update_title();
        }
            break;
		case 'b':
		{
			angle -= 0.1;
		}
			break;
		case 'h':
		{
			angle += 0.1;
		}
			break;
		case 'f':
		{
			angle2 -= 0.1;
		}
			break;
		case 'g':
		{
			angle2 += 0.1;
		}
			break;
		case 'k':
		{
			sph->reset();
		}
			break;
		case 'w':
		{
			sph->move_wall();
		}
			break;
		case 'j':
		{
			identify_fluid(50);
		}
			break;
		case 'u':
		{
			identify_fluid(220);
		}
		break;
    }
}

void UI::visible(int v)
{
    if(v==GLUT_VISIBLE)
        glutIdleFunc(animate_);
    else
        glutIdleFunc(0);
}

void UI::stop()
{
    CONTINUOUS = false;
    update_title();
    glutPostRedisplay();
}

void UI::start(const std::string& log_folder_name)
{
    update_title();
    glutPostRedisplay();
}

void UI::identify_fluid(double threshold) 
{
	
	for (auto te = dsc->tetrahedra_begin(); te != dsc->tetrahedra_end(); te++) {

		vec3 center = dsc->get_barycenter(dsc->get_nodes(te.key()));
		center = vec3(center[0] + 1, center[1] + 1, center[2] + 1);
		center = center * sph->get_down_scale();
		vector<Particle> close_particles = sph->get_close_particles_to_pos(center);
		if (close_particles.size() > 0) {
			
			Particle tmp_particle = Particle(center, -1);
			float mass_at_point = sph->calculate_density(tmp_particle, close_particles) * 1000000.0;
			//cout << mass_at_point << endl;
			if (mass_at_point > threshold) {
				vel_fun->set_label(*dsc, te.key(), 1);
			}
		}
		else {
			//vel_fun->set_label(*dsc, te.key(), 0);
		}

	}

}
