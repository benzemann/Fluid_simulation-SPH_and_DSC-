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
#include "Debug_GUI.h"
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
int tmp = 0;
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
		
		std::chrono::duration<real> delta_time = std::chrono::system_clock::now() - last_time;
		last_time = std::chrono::system_clock::now();
		if (show_compute_time) {
			frames_compute_time++;
			dt_compute_time += delta_time.count();
			if (frames_compute_time > 100) {
				double ms_per_frame = (dt_compute_time * 1000.0) / frames_compute_time;
				cout << ms_per_frame << " ms per frame";
				dt_compute_time = 0.0;
				frames_compute_time = 0;
				double FPS = 1000.0 / ms_per_frame;
				int particle_count = sph->get_no_of_particle();
				cout << " [" << FPS << " FPS] [" << particle_count << " particles]" << endl;
			}
		}

		double d_t = delta_time.count();

		//if (d_t > 0.008)
			//d_t = 0.008;
		//d_t = 0.008;

		sph->set_delta(d_t);
		
		sph->CFL_delta_time_update();

		if (sph->get_delta() > 0.009) {
			sph->set_delta(0.009);
		}

		d_t = sph->get_delta();

		sph->update(d_t);
		
		sph->update_velocity(d_t);
		
		sph->update_isosurface();

		if (sph->is_it_projecting()) {
			sph->project_velocities(*dsc);
		}
		//sph->draw_dsc_velocities(*dsc);
		if (sph->is_fluid_detection())
			identify_fluid(100.0);
		if (sph->is_dsc_tracking()) {
			//if(tmp % 2 == 0)
			vel_fun->take_time_step(*dsc);
		}
		if (sph->is_density_correction()) {
			//sph->correct_divergence_error();
			sph->correct_density_error();
		}
		if (tmp % 100 == 0)
			cout << "Time step: " << tmp << endl;
		if (tmp == 1500) {
			sph->write_volume_file();
		}
		else if (tmp > 500) {
			sph->log_volume();
		}

		sph->update_position(d_t);


		tmp++;
		

		if (sph->is_it_projecting()) {
			//sph->project_velocities_to_particles(*dsc);
		}

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
        case 'g':
            stop();
            break;
		case 'h':
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
			if (dsc)
				show_debug_ui();
			break;
        case 'm':
            if(vel_fun)
            {
                std::cout << "MOVE" << std::endl;
				sph->update(0.1);
               // vel_fun->take_time_step(*dsc);
            }
            break;
		case 'n':
			sph->create_collision_box_at_mouse_pos();
			break;
		case 'p':
			sph->create_particle_at_mouse_pos();
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
			selection++;
			if (selection > (number_of_user_flags+number_of_user_variables)) {
				selection = 0;
			}
			if(dsc)
				show_debug_ui();
            /*if(dsc)
            {
                std::cout << "TAKING SCREEN SHOT" << std::endl;
                Painter::save_painting(WIN_SIZE_X, WIN_SIZE_Y, "LOG");
            }*/
            break;
		case 'w':
			selection--;
			if (selection < 0) {
				selection = (number_of_user_flags + number_of_user_variables);
			}
			if (dsc)
				show_debug_ui();
			break;
		case 'd':
			if (dsc) {
				change_selection(1);
				show_debug_ui();
			}	
			break;
		case 'a':
			if (dsc) {
				change_selection(-1);
				show_debug_ui();
			}
			break;
		case 'e':
			if (dsc) {
				set_selection();
				show_debug_ui();
			}
			break;
		case 'f':
			identify_fluid(50.0);
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
#define SHOW(a) std::cout << #a << ": " << (a) << std::endl
void UI::show_debug_ui() {
	system("cls");
	cout << "##### DEBUG SCREEN #####" << endl;
	cout << endl;
	cout << "CHANGE USER VARIABLES:" << endl;
	cout << " - Use 'w' and 's' to select variable" << endl;
	cout << " - Use 'd' to increase value and 'a' to decrease" << endl;
	cout << endl;
	string var_names[6] = {
		"COEFFICIENT_OF_RESTITUTION",
		"COEFFICIENT_OF_FRICTION",
		"KERNEL_RADIUS",
		"GAS_CONSTANT",
		"REST_DENSITY",
		"VISCOCITY_TERM"
	};
	
	for (int i = 0; i < number_of_user_variables; i++) {
		cout << var_names[i] << ": " << *user_variables_ptr[i];
		if (selection == i) {
			cout << "      <<";
		}
		cout << endl;
	}

	string flag_names[12] = {
		"Draw kernel",
		"Draw velocities",
		"Draw viscocity",
		"Draw pressure",
		"Draw surface tension",
		"Draw dsc velocity",
		"Draw isosurface",
		"Using spatial data grid",
		"Project velocities to dsc",
		"Enabled dsc tracking",
		"Enable fluid detection",
		"Enable density correction"
	};
	for (int i = 0; i < number_of_user_flags; i++) {
		cout << flag_names[i] << ": " << *user_flags_ptr[i];
		if (selection == i + number_of_user_variables) {
			cout << "      <<";
		}
		cout << endl;
	}
	cout << "Show ms per frame and FPS: " << show_compute_time;
	if (selection == number_of_user_flags+ number_of_user_variables) {
		cout << "    <<";
	}
	cout << endl;
	cout << endl;
	cout << "OTHER COMMANDS: " << endl;
	cout << " - SPACE : Start/Stop simulation" << endl;
	cout << " - r : Reset simulation" << endl;
	cout << " - b : Move brown collision box to mouse position" << endl;
	cout << " - n : Create new collision box at mouse position" << endl;
	cout << " - f : Identy fluid" << endl;
	cout << " - ESC : Quit program" << endl;
	cout << " - c : Show this information" << endl;
}

void UI::change_selection(double value) {
	if (selection < number_of_user_variables) {
		*user_variables_ptr[selection] += (*user_variables_ptr[selection]*0.1)*value;
	}
	else if (selection >= number_of_user_variables && selection < (number_of_user_variables+ number_of_user_flags)) {
		*user_flags_ptr[selection - number_of_user_variables] = !*user_flags_ptr[selection-number_of_user_variables];
	}
	else {
		show_compute_time = !show_compute_time;
	}
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
		sph->draw_pressure();
		sph->draw_viscocity();
		sph->draw_surface_tension();
		sph->draw_dsc_velocities(*dsc);
		sph->draw_isosurface();
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

	
	DISCRETIZATION = 16; // super high res: 5, high res 10, normal 18, coarse 25

	int width = WIN_SIZE_X - (2 * DISCRETIZATION);
	int height = WIN_SIZE_Y - (2 * DISCRETIZATION);

	std::vector<real> points;
	std::vector<int> faces;

	Trializer::trialize(width, height, DISCRETIZATION, points, faces);

	DesignDomain *domain = new DesignDomain(DesignDomain::RECTANGLE, width, height, DISCRETIZATION);
	
	sph = std::unique_ptr<SPH>(new SPH(100));
	sph->init();
	dsc = std::unique_ptr<DSC2D::DeformableSimplicialComplex>(new DSC2D::DeformableSimplicialComplex(DISCRETIZATION, points, faces, domain));
	user_variables_ptr = sph->get_user_variables_ptr();
	user_flags_ptr = sph->get_user_flags_ptr();
	vel_fun = std::unique_ptr<VelocityFunc<>>(new track_particle_function(VELOCITY, ACCURACY, *sph));
	identify_fluid(20.0);
	reshape(width + 2 * DISCRETIZATION, height + 2 * DISCRETIZATION);

	start();
}

struct face_density {
	HMesh::FaceID face_id;
	double density;
	int par_count;
};

struct by_density {
	bool operator()(face_density const &a, face_density const &b) {
		return a.density > b.density;
	}
};

struct by_count {
	bool operator()(face_density const &a, face_density const &b) {
		return a.par_count > b.par_count;
	}
};
void UI::identify_fluid(float treshold) {
	
	// IDENTIFY FLUID BY PARTICLE DENSITY //
	for (auto fi = dsc->faces_begin(); fi != dsc->faces_end(); ++fi) {
		//if (dsc->is_outside(*fi)) {
		std::vector<DSC2D::vec2> vertices = dsc->get_pos(*fi);
		DSC2D::vec2 center = DSC2D::vec2((vertices[0][0] + vertices[1][0] + vertices[2][0]) / 3.0, (vertices[0][1] + vertices[1][1] + vertices[2][1]) / 3.0);
		Particle tmp_particle = Particle(center*2.0, -1);
		tmp_particle.mass = 0.0;
		vector<Particle> close_particles = sph->get_close_particles_to_pos(tmp_particle.pos);
		if (close_particles.size() > 0) {
			float mass_at_point = sph->calculate_density(tmp_particle, close_particles) * 1000000.0;
			// if mass is greater than the threshold value it is marked as inside the fluid
			if (mass_at_point > treshold) {
				//Particle p = sph->get_closest_particle(center);
				//DSC2D::vec2 p_to_center = center - (p.pos * sph->get_scale());
				//p_to_center.normalize();
				//DSC2D::vec2 n = p.surface_tension;
				//n.normalize();

				//double d = CGLA::dot(n, p_to_center);

				//if (d >= 0.0 || p.surface_tension.length() < 20.0) {
				ObjectGenerator::set_water(*dsc, *fi, 1); // Sets the face with the face id to inside the water
				//}
			}
		}
		else {
			ObjectGenerator::set_water(*dsc, *fi, 0);
		}
	}
}


using namespace DSC2D;



