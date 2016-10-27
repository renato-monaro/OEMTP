/**
 *  Copyright 2015 by Renato Monaro and Silvio Giuseppe
 *
 * This file is part of Open Electromagnetic Transient Program - OEMTP.
 * 
 * OEMTP is free software: you can redistribute 
 * it and/or modify it under the terms of the GNU General Public 
 * License as published by the Free Software Foundation, either 
 * version 3 of the License, or (at your option) any later version.
 * 
 * Some open source application is distributed in the hope that it will 
 * be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @license GPL-3.0+ <http://spdx.org/licenses/GPL-3.0+>
 */

#ifndef COMPONENTS_H
#define COMPONENTS_H


#include <iomanip>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <vector>
#include <math.h>
#include <complex.h>      // std::complex
#include <iostream>
#include <fstream>
#include <string>

#include "constants.h"




using namespace std;

namespace oemtp{

	class Transformer:public Component{
		public:
			Transformer(){};
			Transformer(string N1a, string N1b, string N2a, string N2b,double L11, double L21,double L21, double L22, double dt);
			bool Compute_Ih(bool e);
			bool Set_Value(double Val);
		protected:
		};
		
	
}
#endif
