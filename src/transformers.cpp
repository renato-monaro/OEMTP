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

#include "transformers.h"

using namespace std;
using namespace oemtp;


			
Transformer::Transformer(string N1a, string N1b, string N2a, string N2b,double Lkk, double Lkm,double Lmk, double Lmm, double dT){
	cout<<"Lkk:"<<Lkk<<" Lkm:"<<Lkm<<" Lmk:"<<Lmk<<" Lmm:"<<Lmm<<endl;
	Alias.push_back(N1a);
	Alias.push_back(N1b);
	Alias.push_back(N2a);
	Alias.push_back(N2b);
	Reff=gsl_matrix_alloc(2,2);
	Geff=gsl_matrix_alloc(2,2);
	Ic=gsl_vector_alloc(2);
	double lup=2*(Lkk*Lmm-Lkm*Lmk);
	dt=dT;
	gsl_matrix_set(Geff,0,0,(dt*(Lmm-Lmk))/lup+(dt*(Lmk))/lup);
	gsl_matrix_set(Geff,1,1,(dt*(Lkk-Lmk))/lup+(dt*(Lmk))/lup);
	gsl_matrix_set(Geff,0,1,-(dt*(Lmk))/lup);
	gsl_matrix_set(Geff,1,0,-(dt*(Lmk))/lup);

	int s;
	gsl_permutation *ppp;
	gsl_matrix *A;

	A=gsl_matrix_alloc(2,2);
	ppp = gsl_permutation_alloc (2);

	gsl_matrix_memcpy (A, Geff);
	gsl_linalg_LU_decomp (A, ppp, &s);
    	gsl_linalg_LU_invert (A, ppp,Reff);	

	gsl_matrix_free(A);
	gsl_permutation_free(ppp);

	Changed=true;
	}
	
bool Transformer::Compute_Ih(bool e){


	gsl_blas_dgemv (CblasNoTrans, 1.0, &View_Gpr.matrix, &View_V_Pri.vector, 0.0, &View_I_Hist_Pri.vector);
	gsl_vector_add(&View_I_Hist_Pri.vector,Ic); //I=I'+I_h
	return true;
	}
	


