/**
 *  Copyright 2015 by Renato Monaro and Silvio Giuseppe
 *  Copyright 2016 by Renato Monaro, Silvio Giuseppe and Heitor Kenzo Koga
 *  Copyright 2017 by Renato Monaro, Silvio Giuseppe and Heitor Kenzo Koga
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
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <vector>
#include <math.h>
#include <complex.h>      // std::complex
#include <iostream>
#include <fstream>
#include <string>

#include "constants.h"
#include "cable.h"



using namespace std;

namespace oemtp{
	class Component{
		public:
			Component();
			~Component();
			virtual bool Compute_Ih(bool e);
			bool Compute_Gpr();
			bool Compute_I();
			double Get_I(unsigned N);
			double Get_V(unsigned N);
			unsigned Get_N_Branches();
			unsigned Get_N_Nodes();
			bool Set_Views(gsl_matrix_view GprV,gsl_vector_view Vv,gsl_vector_view IhV);
			string Get_Alias(unsigned Al);
			bool G_Changed();
			bool Set_Value(unsigned k, unsigned l,double);
			virtual bool Reset();
		protected:
			gsl_matrix_view View_Gpr;
			gsl_vector_view View_V_Pri;
			gsl_vector *Ic;
			gsl_vector_view View_I_Hist_Pri;
			vector<string> Alias;
			gsl_matrix *Reff;
			bool Changed;
			double dt;
			gsl_matrix *Refflu;
			gsl_permutation *p;
		};
	class Resistor:public Component{
		public:
			Resistor(){};
			Resistor(string N1,string N2,double Res);
			Resistor(vector<string> N,vector<vector<double> > Res);
			bool Compute_Ih(bool e);
			bool Set_Value(double Val);
		protected:
		};
		
	class Switch:public Resistor{
		public:
			Switch(){};
			Switch(string N1,string N2,bool IniStatus,double ResON, double ResOFF);
			Switch(string N1,string N2,bool IniStatus);
			Switch(string N1,string N2);
			bool Set_Ron(double R);
			bool Set_Roff(double R);
			void Open();
			void Close();
			bool Get_Status();
			void Set_Status(bool);
		protected:
			bool Status;
			double Ron, Roff;
		};
		
	class Switch2:public Switch{
		public:
			Switch2(string N1,string N2,bool IniStatus,double ResON, double ResOFF);
			bool Compute_Ih(bool e);
		};
		
	class Diode:public Switch{
		public:
			Diode(string N1,string N2):Switch(N1,N2){};
			Diode(string N1,string N2,double ResON, double ResOFF):Switch(N1,N2,false, ResON, ResOFF){};
			bool Compute_Ih(bool e);
		};		
	class Voltmeter:public Resistor{
		public:
			Voltmeter(string N1,string N2);
			Voltmeter(vector<string> N);
		protected:
		};
	class Ammeter:public Resistor{
		public:
			Ammeter(string N1,string N2);
			Ammeter(vector<string> N);			
		protected:
		};
			
	class Inductor:public Component{
		public:
			Inductor(string N1,string N2,double L, double dT);
			bool Compute_Ih(bool e);
			bool Set_Value(double Val);
		};
			
	class Capacitor:public Component{
		public:
			Capacitor(string N1,string N2,double C, double dT);
			bool Compute_Ih(bool e);
			bool Set_Value(double Val);
		};
		
	class Current_Source:public Component{
		public:
			Current_Source(){};
			Current_Source(string N1,string N2,double Ic, double R);
			bool Compute_Ih(bool e);
		protected:
			double In;
		};	
			
	class DC_Source:public Current_Source{
		public:
			DC_Source(string N1,string N2,double Vol,double Res);
		};
		
	class AC_Source:public Current_Source{
		public:
			AC_Source(string N1,string N2,double Vol,double F, double A,double Res,double dT);
			bool Compute_Ih(bool e);
			void Set_Voltage(double);
			void Set_Frequency(double); 
			void Set_Angle(double);
			protected:
			double Frequency;
			double Phase;
			double T;
		};
		
	class Lossless_Line: public Component{
		public:
			Lossless_Line(string N1,string N2,double d, double l, double c, double dT);
			Lossless_Line(string N1,string N2,double d, double xl, double yc, double f, double dT);
			Lossless_Line(string N1,string N2, double Zc, double Tau, double dT);
			bool Compute_Ih(bool e);
			bool Set_Value(double p);
			bool Reset();
		protected:
			double Zc;
			double Tau;
			double k_eff;
			unsigned k_inf,k_sup;
			gsl_vector *Vm,*Vk,*Ikm,*Imk;
			double Vk_tau,Ikm_tau, Vm_tau,Imk_tau, alpha;
		};
		
	class Line: public Component{
		public:
			Line(string N1,string N2,double d, double l, double c, double r, double dT);
			Line(string N1,string N2, double Zc, double Tau, double r, double dT);
			Line(string N1,string N2,double d, double xl, double yc, double f, double r, double dT);
			bool Compute_Ih(bool e);
			bool Set_Value(double p);
			bool Reset();
		protected:
			double Zc;
			double R;
			double Tau;
			double k_eff;
			unsigned k_inf,k_sup;
			gsl_vector *Vm,*Vk,*Ikm,*Imk;
			double Vk_tau,Ikm_tau, Vm_tau,Imk_tau, alpha;
		};
		
//	class MLine: public Component{
//		public:
//			MLine(CableSet *C, vector<string> N, double l, double dT);
//			bool Compute_Ih(bool e);
//		protected:

//			double Zc;
//			double R;
//			double Tau;
//			unsigned k_inf,k_sup;
//			
//			
//			unsigned nModes;
//			vector<double> Zc_Mode;
//			vector<double> R_Mode;
//			vector<double> Tau_Mode;
//			vector<double> k_eff_Mode;
//			vector<unsigned> k_inf_Mode,k_sup_Mode;
//			gsl_matrix *Vm_Mode,*Vk_Mode,*Ikm_Mode,*Imk_Mode;
//			
//			gsl_vector *V_Mode;
//			gsl_vector *I_Mode;
//			gsl_vector *Ih_Mode;
//			gsl_matrix *Tv_inv,*Ti;			
//			vector<double> Vk_tau_Mode,Ikm_tau_Mode, Vm_tau_Mode,Imk_tau_Mode, alpha_Mode;

//			gsl_vector_view Vm,Vk,Ikm,Imk;
//			double Vk_tau,Ikm_tau, Vm_tau,Imk_tau, alpha;
//			
//			gsl_matrix *Zc_Mod,*R_Mod;
//			gsl_matrix *Zc_Phase,*R_Phase;
//			
//			gsl_matrix_view Reff_Send,Reff_Recv;
//			gsl_matrix_view Ti_Send,Ti_Recv;
//			gsl_matrix_view Tv_Send,Tv_Recv;
//			gsl_matrix *Gpr_Mode;

//		};
		
	class MLine2: public Component{
		public:
			MLine2(CableSet *C, vector<string> N, double l, double dT);
			bool Compute_Ih(bool e);
			bool Reset();
		protected:
			gsl_matrix *Tv_inv,*Ti;			
			gsl_matrix *Zc_Phase,*R_Phase;
			vector<Line*> Lines;
    	    gsl_matrix *Gpr_MODE;
			gsl_vector *V_Pri_MODE;
			gsl_vector *I_Hist_Pri_MODE;
			unsigned nModes;
		};

struct IM_ODE_Data{
	double u1a, u1b, u1c;
	double R1;
	double R2;
	double L1;
	double L2;
	double LH;
	double sig;
	double JJ;
	double KD;
	unsigned P;
	double dt;
	double MT;
	double w1, wl;
};

	class InductionMachine:public Component{
		public:
			InductionMachine(string N1, string N2, string N3, double r1, double r2, double l1, double l2, double lh, double jj, double kd, int p, double mt, double Dt);
			bool Compute_Ih(bool e);
			static int func (double t, const double Y[], double F[], void *params);
			double Get_Speed();
			void Set_Mec_Torque(double);
			double Get_Torque();
		protected:
			struct IM_ODE_Data data;
			double y[7];
			double T;
		};
	
	
}
#endif
