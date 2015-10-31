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

#ifndef ODE_H
#define ODE_H

#include "components.h"
#include <gsl/gsl_odeiv2.h>


using namespace std;

namespace oemtp{

struct ODE_Data{
double *parametros;
double *fonte;
unsigned n_param;
unsigned n_fonte;
};

struct InductionMachineData{
	gsl_vector *V;
	gsl_vector *C;
	gsl_matrix *A;
	gsl_matrix *B;
	double *aux;
	};

class ODE:public Current_Source{
	public:
		ODE(){};
		~ODE();
		bool Compute_Ih(bool e);
	protected:
		gsl_odeiv2_system sys;
		struct ODE_Data Parameters;
		gsl_odeiv2_driver *d;
		double dT;
		double *y;
		double T;
	private:

	};
class Inductor_ODE:public ODE{
	public:
		Inductor_ODE(string N1, string N2, double L, double Dt);	
	private:
		static int func(double t, const double y[], double f[],void *params);
	};
	
	
class DFIG_AB0:public Current_Source{
	public:
		DFIG_AB0(string N1,string N2,string N3,string Nn,string Nr1,string Nr2,string Nr3,string Nrn,double Rs, double Rr, double Lm, double Ls, double Lr, unsigned pp, double J, double D, double dT);
		DFIG_AB0(string N1,string N2,string N3,string Nn,double Rs, double Rr, double Lm, double Ls, double Lr, unsigned pp, double J, double D, double dT);
		bool Compute_Ih(bool);
		bool Set_Views(gsl_matrix_view GprV,gsl_vector_view Vv,gsl_vector_view IhV);
		double Get_Speed();
		void Set_Torque(double);
		double Get_Torque();
		double Get_Angle();
		void Set_Speed(double);
	protected:
		static int func(double t, const double y[], double f[],void *params);
		struct InductionMachineData Data;
		gsl_matrix *Clarke,*ClarkeI,*Phi_I,*M,*MI;
		gsl_vector *States,*Iab0,*Temp;
		gsl_vector_view Vab0_s,Vab0_r,Torque,Phi,W,Theta,Iab0_r,Iab0_s;
		gsl_vector_view Vabc_s,Vabc_r,Iabc_r,Iabc_s;
		gsl_vector *Vabc_ant;
		double dT,T;
		gsl_odeiv2_driver *d;
		gsl_odeiv2_system sys;
		unsigned pp;
	};
	
class Synchronous_Machine:public ODE{
	public:
		Synchronous_Machine(string Na, string Nb, string Nc, string Nn,string NF1,string NF2, double L, double Dt);	
	protected:
		gsl_matrix *L;
		gsl_matrix *R;
		
	private:
		static int func(double t, const double y[], double f[],void *params);
	};
}
#endif
