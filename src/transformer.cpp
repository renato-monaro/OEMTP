/**
 *  Copyright 2016 by Renato Monaro
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

#include "components.h"

using namespace std;
using namespace oemtp;


			
Transformer::Transformer(string N1a, string N1b, string N2a, string N2b,double Lkk, double Lkm,double Lmk, double Lmm, double dT);
	Alias.push_back(N1a);
	Alias.push_back(N1b);
	Alias.push_back(N2a);
	Alias.push_back(N2b);
	Reff=gsl_matrix_alloc(2,2);
	Ic=gsl_vector_alloc(2);
	double lup=2*(Lkk*Lmm-Lkm*Lmk);
	dt=dT;
	gsl_matrix_set(Reff,0,0,lup/(dt*(Lmm-Lmk))+lup/(dt*(Lmk)));
	gsl_matrix_set(Reff,1,1,lup/(dt*(Lkk-Lmk))+lup/(dt*(Lmk)));
	gsl_matrix_set(Reff,0,1,-lup/(dt*(Lmk)));
	gsl_matrix_set(Reff,1,0,-lup/(dt*(Lmk)));


	Changed=true;
	}
	
bool Transformer::Compute_Ih(bool e){
	if(e){
		gsl_vector_memcpy(&View_I_Hist_Pri.vector,Ic);
		}
	else{
		gsl_blas_dgemv (CblasNoTrans, 1.0, &View_Gpr.matrix, &View_V_Pri.vector, 0.0, &View_I_Hist_Pri.vector);
		gsl_vector_add(&View_I_Hist_Pri.vector,Ic); //I=I'+I_h
		}
	return true;
	}
	


