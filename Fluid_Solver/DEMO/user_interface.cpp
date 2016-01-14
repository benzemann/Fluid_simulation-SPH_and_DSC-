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

#include "track_particle_function.h"
#include "Particle_System.h"
#include "SPH.h"
#include "trializer.h"
#include "object_generator.h"
#include "draw.h"
#include <iostream>
#include <chrono>

using namespace std;

void _check_gl_error(const char *file, int line)
{
    GLenum err (glGetError());
    
    while(err!=GL_NO_ERROR) {
        std::string error;
        
        switch(err) {
            case GL_INVALID_OPERATION:      error="INVALID_OPERATION";      break;
            case GL_INVALID_ENUM:           error="INVALID_ENUM";           break;
            case GL_INVALID_VALUE:          error="INVALID_VALUE";          break;
            case GL_OUT_OF_MEMORY:          error="OUT_OF_MEMORY";          break;
            case GL_INVALID_FRAMEBUFFER_OPERATION:  error="INVALID_FRAMEBUFFER_OPERATION";  break;
        }
        
        std::cerr << "GL_" << error.c_str() <<" - "<<file<<":"<<line<<std::endl;
        err=glGetError();
    }
}

#define check_gl_error() _check_gl_error(__FILE__,__LINE__)


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

UI::UI(int &argc, char** argv)
{
    instance = this;
	WIN_SIZE_X = 500;
    WIN_SIZE_Y = 500;

    glutInit(&argc, argv);
    glutInitWindowSize(WIN_SIZE_X,WIN_SIZE_Y);
#ifdef WIN32
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_MULTISAMPLE);
#else
    glutInitDisplayString("rgba double samples=16");
#endif
    glutCreateWindow("");
    
    glEnable(GL_MULTISAMPLE);
	
	glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_POINT_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glutDisplayFunc(display_);
    glutKeyboardFunc(keyboard_);
    glutVisibilityFunc(visible_);
    glutReshapeFunc(reshape_);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    dsc = nullptr;
    vel_fun = nullptr;
    
    if(argc > 1)
    {
        QUIT_ON_COMPLETION = true;
        CONTINUOUS = true;
        RECORD = true;
        
        ::Util::ArgExtracter ext(argc, argv);
        ext.extract("nu", VELOCITY);
        ext.extract("delta", DISCRETIZATION);
        ext.extract("alpha", ACCURACY);
    }
    else {
        VELOCITY = 10.;
        DISCRETIZATION = 25.;
        ACCURACY = 5.;
        
        CONTINUOUS = false;
        RECORD = true;
        QUIT_ON_COMPLETION = false;
    }
    update_title();

    check_gl_error();

	track_particle_function* instance = static_cast<track_particle_function*>(vel_fun.get());
	
}

void UI::update_title()
{
    std::ostringstream oss;
    oss << "2D DSC\t";
    if(vel_fun)
    {
        oss << vel_fun->get_name();
        oss << ", Time step " << vel_fun->get_time_step();
    }
    oss << " (Nu = " << VELOCITY << ", Delta = " << DISCRETIZATION << ", Alpha = " << ACCURACY << ")";
    std::string str(oss.str());
    glutSetWindowTitle(str.c_str());
}
double sum_time = 0.0;

void UI::display()
{
    if (glutGet(GLUT_WINDOW_WIDTH) != WIN_SIZE_X || glutGet(GLUT_WINDOW_HEIGHT) != WIN_SIZE_Y) {
        return;
    }
	
    draw();
    update_title();
	
    if(vel_fun && CONTINUOUS)
    {
		double time = vel_fun->get_compute_time() + vel_fun->get_deform_time();
		
		if (vel_fun->get_time_step() >= 100) {
			sum_time += time;
			cout << "avg time: " << sum_time/ (vel_fun->get_time_step()-100) << endl;
			cout << "total time: " << sum_time << endl;
		}
		std::chrono::duration<real> delta_time = std::chrono::system_clock::now() - last_time;
		last_time = std::chrono::system_clock::now();
		sph->update(delta_time.count());
        basic_log->write_timestep(*vel_fun);
        if (vel_fun->is_motion_finished(*dsc))
        {
            //stop();
            if (QUIT_ON_COMPLETION) {
                exit(0);
            }
        }
    }
    check_gl_error();
}

void UI::reshape(int width, int height)
{
    if(dsc)
    {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluOrtho2D(0, WIN_SIZE_X, 0, WIN_SIZE_Y);
    }
    
    glViewport(0, 0, WIN_SIZE_X, WIN_SIZE_Y);
    glutReshapeWindow(WIN_SIZE_X, WIN_SIZE_Y);
}

void UI::animate()
{
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
            break;
		case '1':
			create_fluid_domain();
			break;
        case ' ':
            if(!CONTINUOUS)
            {
				last_time = std::chrono::system_clock::now();
                std::cout << "MOTION STARTED" << std::endl;
            }
            else {
                std::cout << "MOTION PAUSED" << std::endl;
            }
            CONTINUOUS = !CONTINUOUS;
            break;
		case 'b':
			sph->move_collision_box(0);
			break;
		case 'c':
			enter_command();
			break;
        case 'm':
            if(vel_fun)
            {
                std::cout << "MOVE" << std::endl;
                vel_fun->take_time_step(*dsc);
            }
            break;
		case 'n':
			sph->create_collision_box_at_mouse_pos();
			break;
        case 't':
            if(vel_fun)
            {
                std::cout << "TEST" << std::endl;
                vel_fun->test(*dsc);
            }
            break;
        case '\t':
            if(dsc)
            {
                switch_display_type();
            }
            break;
		case 'r':
			sph->reset();
			break;
        case 's':
            if(dsc)
            {
                std::cout << "TAKING SCREEN SHOT" << std::endl;
                Painter::save_painting(WIN_SIZE_X, WIN_SIZE_Y, "LOG");
            }
            break;
		case 'z':
			sph->subtract_costum();
			break;
		case 'x':
			sph->add_costum();
			break;
        case '+':
            if(!vel_fun)
            {
                VELOCITY = std::min(VELOCITY + 1., 100.);
                update_title();
            }
            break;
        case '-':
            if(!vel_fun)
            {
                VELOCITY = std::max(VELOCITY - 1., 0.);
                update_title();
            }
            break;
        case '.':
            if(!dsc)
            {
                DISCRETIZATION = std::min(DISCRETIZATION + 0.5, 100.);
                update_title();
            }
            break;
        case ',':
            if(!dsc)
            {
                DISCRETIZATION = std::max(DISCRETIZATION - 0.5, 1.);
                update_title();
            }
            break;
        case '<':
            if(!vel_fun)
            {
                ACCURACY = std::min(ACCURACY + 1., 100.);
                update_title();
            }
            break;
        case '>':
            if(!vel_fun)
            {
                ACCURACY = std::max(ACCURACY - 1., 1.);
                update_title();
            }
            break;
    }
}

void UI::enter_command() {
	cout << "Enter command:" << endl;
	string command;
	double input = 0.0;
	cin >> command;
	if (command == "help") {
		cout << "Command list: " << endl;
		cout << "    - help : Show the command list" << endl;
		cout << "    - kernel_radius radius : Change the kernel radius to input value" << endl;
		cout << "    - show_kernel_radius true/false : Show/hide kernel radius, input: 1=show, 0=hide" << endl;
	}
	else if (command == "kernel_radius") {
		cin >> input;
		cout << "Kernel radius changed to " << input << endl;
		sph->set_kernel_radius(input);
	}
	else if (command == "show_kernel_radius") {
		sph->change_draw_kernel_radius();
	}
	else if (command == "show_velocities") {
		sph->change_draw_velocities();
	}
	else if(command != "done") {
		cout << "There is no command for: " << command << endl;
	}
	if(command != "done")
		enter_command();
}

void UI::visible(int v)
{
    if(v==GLUT_VISIBLE)
        glutIdleFunc(animate_);
    else
        glutIdleFunc(0);
}

void UI::draw()
{
    Painter::begin();

    if (dsc)
    {
        Painter::draw_complex(*dsc);
		//Painter::draw_vertices_index(*dsc);
		sph->draw_particles();
		sph->draw_collision_boxes();
		sph->draw_kernel_radius();
		sph->draw_velocities();
        if(RECORD && CONTINUOUS)
        {
            Painter::save_painting(WIN_SIZE_X, WIN_SIZE_Y, basic_log->get_path(), vel_fun->get_time_step());
        }
    }
    Painter::end();
	
}

void UI::stop()
{
    draw();
    
    if(basic_log)
    {
        basic_log->write_message("MOTION STOPPED");
        basic_log->write_log(*dsc);
        basic_log->write_log(*vel_fun);
        basic_log->write_timings(*vel_fun);
    }
    
    basic_log = nullptr;
    vel_fun = nullptr;
    dsc = nullptr;
}

void UI::start()
{
    basic_log = std::unique_ptr<Log>(new Log(create_log_path()));
    basic_log->write_message(vel_fun->get_name().c_str());
    basic_log->write_log(*dsc);
    basic_log->write_log(*vel_fun);
    
    update_title();
}

void UI::create_fluid_domain()
{
	stop();

	DISCRETIZATION = 18;
	int width = WIN_SIZE_X - (2 * DISCRETIZATION);
	int height = WIN_SIZE_Y - (2 * DISCRETIZATION);

	std::vector<real> points;
	std::vector<int> faces;

	Trializer::trialize(width, height, DISCRETIZATION, points, faces);

	DesignDomain *domain = new DesignDomain(DesignDomain::RECTANGLE, width, height, DISCRETIZATION);
	
	sph = std::unique_ptr<SPH>(new SPH(100));
	sph->init();
	dsc = std::unique_ptr<DeformableSimplicialComplex>(new DeformableSimplicialComplex(DISCRETIZATION, points, faces, domain));
	
	vel_fun = std::unique_ptr<VelocityFunc<>>(new track_particle_function(VELOCITY, ACCURACY));

	reshape(width + 2 * DISCRETIZATION, height + 2 * DISCRETIZATION);

	start();
}

using namespace DSC2D;



