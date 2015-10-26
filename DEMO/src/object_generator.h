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

#pragma once

#include "util.h"
#include "DSC.h"

namespace DSC2D {
    
    class ObjectGenerator {
        
        static void fit_mesh_to_object(DeformableSimplicialComplex& dsc, const std::vector<vec2>& corners);
        
        static void label_faces(DeformableSimplicialComplex& dsc, const std::vector<vec2>& corners, int label);
        
        static void create_object(DeformableSimplicialComplex& dsc, const std::vector<vec2>& corners, int label);
        
    public:
        
        static void create_blob(DeformableSimplicialComplex& dsc, const vec2& center, const real& radius, int label)
        {
            std::vector<vec2> corners;
            for (real a = 0; a < 2*M_PI; a += (1./32.)*2*M_PI)
            {
                corners.push_back(radius*vec2(std::cos(a), -std::sin(a)) + center);
            }
            create_object(dsc, corners, label);
        }
        
        static void create_square(DeformableSimplicialComplex& dsc, const vec2& origin, const vec2& size, int label)
        {
            std::vector<vec2> corners;
            corners.push_back(origin);
            corners.push_back(origin + vec2(0., size[1]));
            corners.push_back(origin + size);
            corners.push_back(origin + vec2(size[0], 0.));
            
            create_object(dsc, corners, label);
        }

		static void set_water(DeformableSimplicialComplex& dsc, HMesh::FaceID fi){
			dsc.update_attributes(fi, 1);
		}

		static void create_water(DeformableSimplicialComplex& dsc, std::vector<vec2> particle_pos){ // not used
			//dsc.init_attributes(); // Maybe not necessary??

			/*Label the faces with particles inside*/
			/*int tmp = 0;
			for (auto fi = dsc.faces_begin(); fi != dsc.faces_end(); fi++)
			{
				// Calculate if particle is inside the face
				for (int i = 0; i < particle_pos.size(); i++){
					auto verts = dsc.get_pos(*fi);
					vec2 particle = particle_pos[i];
					vec2 v1 = vec2(verts[1] - verts[0]);
					vec2 v2 = vec2(verts[2] - verts[0]);
					double a = (CGLA::cross(particle, v2) - CGLA::cross(verts[0], v2)) / CGLA::cross(v1, v2);
					double b = -((CGLA::cross(particle, v1) - CGLA::cross(verts[0], v1)) / CGLA::cross(v1, v2));

					if (a > 0.0 && b > 0 && (a + b) < 1.0){
						dsc.update_attributes(*fi, 1); // Label face as inside 
					}
				}
			}*/

			dsc.update_attributes(); // Maybe not necessary??

			dsc.fix_complex();
			
			for (auto vi = dsc.vertices_begin(); vi != dsc.vertices_end(); ++vi)
			{
				if (dsc.is_movable(*vi))
				{
					std::vector<int> labels = dsc.get_interface_labels(*vi);
					if (labels[0] == 0)
					{
						
						dsc.set_destination(*vi, dsc.get_pos(*vi) + 5.0 *  dsc.get_normal(*vi));
					}
				}
			}
			
			dsc.deform();

		}

    };
    
    
}