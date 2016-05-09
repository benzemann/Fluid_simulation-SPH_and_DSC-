#pragma once

#include "DSC.h"

namespace DSC {
	template <typename DeformableSimplicialComplex = DeformableSimplicialComplex<>>
	class ObjectGenerator {

	public:
		//template <typename DeformableSimplicialComplex = DeformableSimplicialComplex<>>
		static void set_label(DSC::DeformableSimplicialComplex &dsc, is_mesh::TetrahedronKey te, int label) {
			dsc->set_label(te, label);
			
		}
		static void test() {

		}
	};
}




