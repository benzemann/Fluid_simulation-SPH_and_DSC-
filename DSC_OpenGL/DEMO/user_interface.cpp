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
	vec3 center = _obj_dim / 2.0;
	vec3 eye = center + vec3(gl_dis_max*2.0*cos(angle)*cos(angle2),
		gl_dis_max*2.0*cos(angle)*sin(angle2),
		gl_dis_max*2.0*sin(angle));

	GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat mat_shininess[] = { 50.0 };
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



	load_model("C:/Users/Jeppe/Desktop/DSC_OpenGL/DSC_OpenGL/DSC.visualstudio/armadillo.dsc", 2.5);

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
		/* Z near */ 10.0, /* Z far */ gl_dis_max*5.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	vec3 center = _obj_dim / 2.0;
	vec3 eye = center + vec3(gl_dis_max*3.0*cos(angle)*cos(angle2),
		gl_dis_max*3.0*cos(angle)*sin(angle2),
		gl_dis_max*3.0*sin(angle));
	vec3 head = vec3(-sin(angle)*cos(angle2),
		-sin(angle)*sin(angle2),
		cos(angle));
	gluLookAt(eye[0], eye[1], eye[2], /* eye is at (0,8,60) */
		center[0], center[1], center[2],      /* center is at (0,8,0) */
		head[0], head[1], head[2]);      /* up is in postivie Y direction */

	int size = std::min(WIN_SIZE_Y, WIN_SIZE_X);
	glViewport((WIN_SIZE_X - size) / 2.0, (WIN_SIZE_Y - size) / 2.0, size, size);

	glClearColor(0.1, 0.1, 0.1, 1.0);
}

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
//	draw_helper::dsc_draw_edge(*dsc);
	draw_helper::dsc_draw_interface(*dsc);

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
            vel_fun = std::unique_ptr<VelocityFunc<>>(new RotateFunc(vel_fun->get_velocity(), vel_fun->get_accuracy()));
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
            painter->set_view_position(camera_pos);
            painter->save_painting("LOG");
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
