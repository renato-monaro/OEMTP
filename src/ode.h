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

#include <gsl/gsl_odeiv2.h>
#include "analysis.h"
using namespace std;

namespace oemt{

struct param {
   double Carac;
   double u;
};

struct InductionMachineData{
	double Rs;
	double Rr;
	double Lm;
	double Ls;
	double Lr;
	double J;
	double D;
	unsigned pp;
	};

//class ODE:public Current_Source{
//	public:
//		ODE();
//		bool Compute_Ih(bool e);
//	protected:
//		virtual static int Function(double t, const double y[], double f[],void *params);
//		virtual static int Jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);

//	gsl_odeiv2_system sys;
//	double mu;
//	};
class Inductor_ODE:public Current_Source{
	public:
		~Inductor_ODE(){};
		Inductor_ODE(string N1,string N2,double L, double dT);
		bool Compute_Ih(bool e);
	protected:
		double L;
		static int Function(double t, const double y[], double f[],void *params);
		static int Jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
		gsl_odeiv2_system sys;
		double mu;
		double h;
		gsl_odeiv2_driver *d;
	};

class Capacitor_ODE:public Current_Source{
	public:
		~Capacitor_ODE(){};
		Capacitor_ODE(double L);
		bool Compute_Ih(bool e);
	protected:
		double L;
		static int Function(double t, const double y[], double f[],void *params);
		static int Jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
		gsl_odeiv2_system sys;
		double mu;
	};
	
class DFIG_AB0:public Current_Source{
	public:
		~DFIG_AB0(){};
		DFIG_AB0(string N1,string N2,string N3,double Rs, double Rr, double Lm, double Ls, double Lr, unsigned pp, double J, double D=0.0);
		bool Compute_Ih(bool e);
	protected:
		bool V_ABC_AB0();
		bool I_AB0_ABC();
		struct InductionMachineData Data;
		static int Function(double t, const double y[], double f[],void *params);
		static int Jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
		gsl_odeiv2_system sys;
	};
}
#endif
