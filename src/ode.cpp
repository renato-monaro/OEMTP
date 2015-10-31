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

#include "ode.h"

using namespace oemt;


//ODE::ODE(){
//}
//int Inductor_ODE::Function (double t, const double y[], double f[], void *params){
//  double mu = *(double *)params;
//  f[0] = 1;
//  return GSL_SUCCESS;
//  }
//int Inductor_ODE::Jacobian (double t, const double y[], double *dfdy, double dfdt[], void *params){
//  return  GSL_SUCCESS;
//	}
//	
//bool  ODE::Compute_Ih(bool e){
//	return true;
//	}
	
Inductor_ODE::Inductor_ODE(string N1,string N2,double Ind, double dT){
	L=Ind;
	gsl_odeiv2_system sys_temp = {Function,NULL,1,&L};
	sys=sys_temp;
	h=dT;
	Alias.push_back(N1);
	Alias.push_back(N2);
	Reff=gsl_matrix_alloc(1,1);
	I=gsl_vector_alloc(1);
	gsl_matrix_set(Reff,0,0,1E9);
	d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45,1e-6, 1e-6, 0.0);
	}

bool Inductor_ODE::Compute_Ih(bool e){
	
	return true;
	}
	
int Inductor_ODE::Function (double t, const double y[], double f[], void *params){
  double ind = *(double *)params;
  f[0] = gsl_vector_get(&View_V_Pri.vector,0)/ind;
  return GSL_SUCCESS;
  }
int Inductor_ODE::Jacobian (double t, const double y[], double *dfdy, double dfdt[], void *params){
  return  GSL_SUCCESS;
	}
		
DFIG_AB0::DFIG_AB0(string N1,string N2,string N3,,string Nn, double rs, double rr, double lm, double ls, double lr, unsigned PP, double j, double d=0.0){

	Alias.push_back(N1);
	Alias.push_back(Nn);
	Alias.push_back(N2);
	Alias.push_back(Nn);
	Alias.push_back(N3);
	Alias.push_back(Nn);
	
	Reff=gsl_matrix_alloc(3,3);
	I=gsl_vector_alloc(3);
	
	gsl_matrix_set(Reff,0,0,1E9);
	gsl_matrix_set(Reff,1,1,1E9);
	gsl_matrix_set(Reff,2,2,1E9);
	
	Data.Rs=rs;
	Data.Rr=rr;
	Data.Lm=lm;
	Data.Ls=ls;
	Data.pp=PP;
	Data.J=j;
	Data.D=d;
	gsl_odeiv2_system sys_temp = {Function,Jacobian,6,&Data};
	}

bool DFIG_AB0::Compute_Ih(bool e){

	return true;
	}

	
	
	
	
	
	Capacitor_ODE::Capacitor_ODE(double Ind){
	L=Ind;
	gsl_odeiv2_system sys_temp = {Function,Jacobian,1,NULL};
	sys=sys_temp;
	}

bool Capacitor_ODE::Compute_Ih(bool e){
	return true;
	}
	
int Capacitor_ODE::Function (double t, const double y[], double f[], void *params){
  double mu = *(double *)params;
  f[0] = 1;
  return GSL_SUCCESS;
  }
int Capacitor_ODE::Jacobian (double t, const double y[], double *dfdy, double dfdt[], void *params){
  return  GSL_SUCCESS;
	}
