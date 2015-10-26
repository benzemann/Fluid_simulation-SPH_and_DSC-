/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "Mat4x4f.h"

namespace CGLA 
{

	Mat4x4f rotation_Mat4x4f(Axis axis, float angle)
	{
		Mat4x4f m(0.0f);

		switch(axis)
			{
			case XAXIS:
				m[0][0] = 1.0f;
				m[1][1] = cos(angle);
				m[1][2] = sin(angle);
				m[2][1] = -sin(angle);
				m[2][2] = cos(angle);
				m[3][3] = 1.0f;
				break;
			case YAXIS:
				m[0][0] = cos(angle);
				m[0][2] = -sin(angle);
				m[2][0] = sin(angle);
				m[2][2] = cos(angle);
				m[1][1] = 1.0f;
				m[3][3] = 1.0f;
				break;
			case ZAXIS:
				m[0][0] = cos(angle);
				m[0][1] = sin(angle);
				m[1][0] = -sin(angle);
				m[1][1] = cos(angle);
				m[2][2] = 1.0f;
				m[3][3] = 1.0f;
				break;
			}
		return m;
	}

	Mat4x4f translation_Mat4x4f(const Vec3f& v)
	{
		Mat4x4f m(0.0f);

		m[0][0] = 1.0f;
		m[1][1] = 1.0f;
		m[2][2] = 1.0f;
		m[3][3] = 1.0f;
  
		m[0][3] = v[0];
		m[1][3] = v[1];
		m[2][3] = v[2];
  
		return m;
	}

	Mat4x4f scaling_Mat4x4f(const Vec3f& v)
	{
		Mat4x4f m(0.0f);

		m[0][0] = v[0];
		m[1][1] = v[1];
		m[2][2] = v[2];
		m[3][3] = 1.0f;
   
		return m;
	}
}
