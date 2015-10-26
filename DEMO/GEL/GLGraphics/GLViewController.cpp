/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "gel_glu.h"
#include "GLViewController.h"

using namespace std;
using namespace CGLA;

namespace GLGraphics
{
    
  void GLViewController::reset_projection()
  {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(FOV_DEG, aspect, znear, zfar);
    glMatrixMode(GL_MODELVIEW);
  }

  GLViewController::GLViewController(int _WINX, int _WINY, const CGLA::Vec3f& centre, float rad)
    : FOV_DEG(53),
      WINX(_WINX), WINY(_WINY), 
      aspect(WINX/float(WINY)),
      button_down(false),
      spin(false),
      ball(centre, rad, WINX, WINY)
  {
    znear = 0.01f*rad;
    zfar  = 3*rad;
    reset_projection();
  }

  void GLViewController::grab_ball(TrackBallAction action, const CGLA::Vec2i& pos)
  {
    ball.grab_ball(action,pos);
    if(action==ZOOM_ACTION)
      set_near_and_far();

    spin = false;
    button_down = true;
    last_action = action;
  }

  void GLViewController::roll_ball(const CGLA::Vec2i& pos)
  {
    static Vec2i old_pos = pos;
    Vec2f dir = Vec2f(pos-old_pos);
    float len = dir.length();
    if (len < TINY)
      return;
    
    ball.roll_ball(pos);
    if(last_action==ZOOM_ACTION)
      set_near_and_far();
    
    spin = len>=1.1f;
    old_pos = pos;  
  }


  void GLViewController::release_ball()
  {
    ball.release_ball();
    if(last_action==ZOOM_ACTION)
      set_near_and_far();
  }

  bool GLViewController::try_spin()
  {
    if(spin && !ball.is_grabbed()) 
    {
      ball.do_spin();
      return true;
    }
    return false;
  }
  
  void GLViewController::set_gl_modelview()
  {
    ball.set_gl_modelview();
  }


  void GLViewController::reshape(int W, int H)
  {
    WINX = W;
    WINY = H;
    aspect = WINX/static_cast<float>(WINY);
    glViewport(0,0,WINX,WINY);
    reset_projection();
    ball.set_screen_window(WINX, WINY);
  }  

  void GLViewController::set_near_and_far()
  {  
    float rad = ball.get_eye_dist();
    znear = 0.01f*rad;
    zfar = 3*rad;
    reset_projection();
  }

  bool GLViewController::load(std::ifstream& ifs)
  {
    if(ifs)
    {
      ifs.read(reinterpret_cast<char*>(this), sizeof(GLViewController));
      reset_projection();
      ball.set_screen_window(WINX, WINY);
      return true;
    }
    return false;
  }
  bool GLViewController::save(std::ofstream& ofs) const
  {
    if(ofs)
    {
      ofs.write(reinterpret_cast<const char*>(this), sizeof(GLViewController));
      return true;
    }
    return false;
   }
}
