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

#include "support/support.h"
#include<math.h>
#include<vector>
#include <iostream>

#define mu_0 4E-7*M_PI //H/km
#define eps_0 8.854187817E-12 //F/m
#define gamma 0.5772156649 //constante de Euler
#define u_0 2E-4 // [H/km]

enum Cable_Models {BERGERON=1, JMARTI};
enum Installation_Mode {AIR=1, BURIED};

class Conductor{
		public:
			Conductor(double Rin, double Rout, double r, double mu);
			Conductor(){};
			double Get_Rout();
			double Get_Rin();
			double Get_mu();
			double Get_rho();
			double _Complex Inside_Impedance(double Freq);
			double _Complex Outside_Impedance(double Freq);
			double _Complex Mutual_Impedance(double Freq);
		protected:
			double rho,mu_r,Rin,Rout;
			bool Grounded;
			};
			
class Insulation{
		public:
			Insulation(double e, double mu);
			Insulation(){};
			double Get_Rout();
			double Get_Rin();
			double Get_mu();
			double Get_eps();
			bool Set_Rout(double r);
			bool Set_Rin(double r);
			double Potential_Coeff(double Freq);
			double _Complex Impedance(double Freq);
		protected:
			double eps,mu_r,Rin,Rout;
			};
			
class Cable{
		public:
			Cable(double R);
			bool Join(Conductor);
			bool Join(Insulation);
			bool Assembly();
			void Print();
			void Print(double);
			double Get_R();
			unsigned Get_N_Conductors();
			Conductor Get_Conductor(unsigned N);
			Insulation Get_Insulation(unsigned N);
		protected:
			bool Reorder_Conductors();
			std::vector<Conductor> vCond;
			std::vector<Insulation> vIso;
			double R;
			};
			
class CableSet{
		public:
			CableSet(double rho, double mu, unsigned Inst);
			~CableSet();
			bool Join(Cable,double x, double y);
			bool Compute_Parameters(double Freq);
			void Print();
			double Get_Speed(unsigned Mode);
			double Get_Zc_Mode(unsigned Mode);
			double Get_Zc_Phase(unsigned Mode);
			gsl_matrix* Get_Zc_Mode();
			gsl_matrix* Get_Zc_Phase();
			gsl_matrix* Get_R_Mode();
			gsl_matrix* Get_R_Phase();
			double Get_R_Mode(unsigned Mode);
			double Get_R_Phase(unsigned Mode);
			unsigned Get_N_Modes();
			gsl_matrix* Get_Ti();
			gsl_matrix* Get_Tv();
			gsl_matrix* Get_Ti_Inverse();
			gsl_matrix* Get_Tv_Inverse();
		protected:
			double _Complex Earth_Impedance(double Freq,unsigned C1);
			double _Complex Earth_Mutual_Impedance(double Freq,unsigned C1,unsigned C2);

			bool Compute_Z(double Freq);
			bool Compute_Y(double Freq);
			bool Modal_Trans(double Freq);
			double Distance(unsigned C1,unsigned C2);
			gsl_matrix_complex *Y,*Z,*Ym,*Zm;
			gsl_matrix *Zcf,*Zcm,*Rf,*Rm;
			gsl_vector_complex *Lambda;//,*Propm;
			gsl_matrix_complex *Tv_cplx,*Ti_cplx;
			gsl_matrix *Tv,*Tv_inv,*Ti,*Ti_inv;
			std::vector<double> Pos_x,Pos_y;
			std::vector<Cable> vCable;
			std::vector<double> Vel;
			double rho;
			double mu_r;
			unsigned Type;
			unsigned N;
			unsigned Installation;
	};

