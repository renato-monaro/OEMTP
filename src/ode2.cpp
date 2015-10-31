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

#include "ode2.h"

using namespace oemtp;

ODE::~ODE(){
  gsl_odeiv2_driver_free (d);
  }
  
bool ODE::Compute_Ih(bool e){
	for(unsigned k=0;k<Parameters.n_fonte;k++)
		Parameters.fonte[k]= gsl_vector_get(&View_V_Pri.vector,k);
      int s = gsl_odeiv2_driver_apply_fixed_step (d, &T, dT, 1, y);
   	for(unsigned k=0;k<sys.dimension;k++)
		gsl_vector_set(&View_I_Hist_Pri.vector,k,y[k]);
      if (s != GSL_SUCCESS){
          printf ("error: driver returned %d\n", s);
          return false;
          }
         
	return true;
	}
  
Inductor_ODE::Inductor_ODE(string N1, string N2, double ind,double Dt){
    Alias.push_back(N1);
	Alias.push_back(N2);
	Reff=gsl_matrix_alloc(1,1);
	Ic=gsl_vector_alloc(1);
	gsl_matrix_set(Reff,0,0,1E9);
	
	dT=Dt;
	y=(double*)malloc(1*sizeof(double));	
	Parameters.parametros=(double*)malloc(1*sizeof(double));
	Parameters.fonte=(double*)malloc(1*sizeof(double));
	Parameters.parametros[0]=ind;
	Parameters.n_fonte=1;
	Parameters.n_param=1;
	
	sys.function=func;
	sys.jacobian=NULL;
	sys.dimension=1;
	sys.params=&Parameters;
	
    d=gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4,1e-3, 1e-8, 1e-8);
    
	T=0;
	for(unsigned k=0;k<Parameters.n_fonte;k++)
		Parameters.fonte[k]=0;
   	for(unsigned k=0;k<sys.dimension;k++)
   		y[k]=0;
	}
	
int Inductor_ODE::func(double t, const double y[], double f[],void *params){
	struct ODE_Data *PARAM = (struct ODE_Data *)params;
	f[0] = (PARAM->fonte[0])/(PARAM->parametros[0]);
	return GSL_SUCCESS;
	}
	
DFIG_AB0::DFIG_AB0(string N1,string N2,string N3,string Nn,double Rs, double Rr, double Lm, double Ls, double Lr, unsigned p, double J, double D, double Dt){
    Alias.push_back(N1);Alias.push_back(Nn);
    Alias.push_back(N2);Alias.push_back(Nn);
    Alias.push_back(N3);Alias.push_back(Nn);
    
    Alias.push_back("TERRA");Alias.push_back("TERRA");
    Alias.push_back("TERRA");Alias.push_back("TERRA");
    Alias.push_back("TERRA");Alias.push_back("TERRA");
	Reff=gsl_matrix_alloc(6,6);
	Ic=gsl_vector_alloc(6);

//	gsl_matrix_set_identity(Reff);
//	gsl_matrix_scale(Reff,R_MAX);	

	gsl_matrix_set(Reff,0,0,2*(Ls-Lm)/Dt);
	gsl_matrix_set(Reff,1,1,2*(Ls-Lm)/Dt);
	gsl_matrix_set(Reff,2,2,2*(Ls-Lm)/Dt);		
	gsl_matrix_set(Reff,3,3,2*(Lr-Lm)/Dt);
	gsl_matrix_set(Reff,4,4,2*(Lr-Lm)/Dt);
	gsl_matrix_set(Reff,5,5,2*(Lr-Lm)/Dt);	

	pp=p;
	
	dT=Dt;
	T=0;
	
	double Sig=1-pow(Lm,2)/(Ls*Lr);
	
	Data.A=gsl_matrix_alloc(6,6);
	gsl_matrix_set_zero(Data.A);
	gsl_matrix_set(Data.A,0,0,-Rs/(Sig*Ls));
	gsl_matrix_set(Data.A,1,1,-Rs/(Sig*Ls));
	gsl_matrix_set(Data.A,2,2,-Rr/(Sig*Lr));
	gsl_matrix_set(Data.A,3,3,-Rr/(Sig*Lr));
	gsl_matrix_set(Data.A,0,2,Rs*Lm/(Sig*Ls*Lr));
	gsl_matrix_set(Data.A,1,3,Rs*Lm/(Sig*Ls*Lr));
	gsl_matrix_set(Data.A,2,0,Rr*Lm/(Sig*Ls*Lr));
	gsl_matrix_set(Data.A,3,1,Rr*Lm/(Sig*Ls*Lr));
	gsl_matrix_set(Data.A,4,4,-D/J);
	gsl_matrix_set(Data.A,5,4,1);
	
	Data.aux=(double*)malloc(3*sizeof(double));
	Data.aux[0]=(3.0/2.0)*(Lm/(Sig*Lr*Ls));
	Data.aux[1]=J;
	
	Data.V=gsl_vector_alloc(6);
	States=gsl_vector_alloc(6);
	
	Vab0_s=gsl_vector_subvector(Data.V,0,2);
	Vab0_r=gsl_vector_subvector(Data.V,2,2);
	Torque=gsl_vector_subvector(Data.V,4,1);
	
	
	Vabc_ant=gsl_vector_alloc(6);
	
	Phi=gsl_vector_subvector(States,0,4);
	W=gsl_vector_subvector(States,4,1);
	Theta=gsl_vector_subvector(States,5,1);
	
	Iab0=gsl_vector_alloc(4);
	Temp=gsl_vector_alloc(2);
	
	Iab0_s=gsl_vector_subvector(Iab0,0,2);
	Iab0_r=gsl_vector_subvector(Iab0,2,2);
	
	Clarke=gsl_matrix_alloc(2,3);
	gsl_matrix_set(Clarke,0,0,1.0);
	gsl_matrix_set(Clarke,0,1,-0.5);
	gsl_matrix_set(Clarke,0,2,-0.5);
	gsl_matrix_set(Clarke,1,0,0);
	gsl_matrix_set(Clarke,1,1,sqrt(3.0/4.0));
	gsl_matrix_set(Clarke,1,2,-sqrt(3.0/4.0));
	
	ClarkeI=gsl_matrix_alloc(3,2);
	gsl_matrix_set(ClarkeI,0,0,1.0);
	gsl_matrix_set(ClarkeI,0,1,0.0);
	gsl_matrix_set(ClarkeI,1,0,-0.5);
	gsl_matrix_set(ClarkeI,1,1,sqrt(3.0/4.0));
	gsl_matrix_set(ClarkeI,2,0,-0.5);
	gsl_matrix_set(ClarkeI,2,1,-sqrt(3.0/4.0));
	
	MI=gsl_matrix_alloc(2,2);
	M=gsl_matrix_alloc(2,2);
	
	Phi_I=gsl_matrix_alloc(4,4);
	gsl_matrix_set(Phi_I,0,0,-Lr/(Ls*Lr-pow(Lm,2)));
	gsl_matrix_set(Phi_I,0,2,Lm/(Ls*Lr-pow(Lm,2)));
	gsl_matrix_set(Phi_I,1,1,-Lr/(Ls*Lr-pow(Lm,2)));
	gsl_matrix_set(Phi_I,1,3,Lm/(Ls*Lr-pow(Lm,2)));
	gsl_matrix_set(Phi_I,2,2,Ls/(Ls*Lr-pow(Lm,2)));
	gsl_matrix_set(Phi_I,2,0,-Lm/(Ls*Lr-pow(Lm,2)));
	gsl_matrix_set(Phi_I,3,3,Ls/(Ls*Lr-pow(Lm,2)));
	gsl_matrix_set(Phi_I,3,1,-Lm/(Ls*Lr-pow(Lm,2)));
	
	sys.function=func;
	sys.jacobian=NULL;
	sys.dimension=6;
	sys.params=&Data;
	
	d=gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd,1e-3, 1e-6, 1e-9);
	}
	
DFIG_AB0::DFIG_AB0(string N1,string N2,string N3,string Nn,string Nr1,string Nr2,string Nr3,string Nrn, double Rs, double Rr, double Lm, double Ls, double Lr, unsigned p, double J, double D,double Dt){
    Alias.push_back(N1);Alias.push_back(Nn);
    Alias.push_back(N2);Alias.push_back(Nn);
    Alias.push_back(N3);Alias.push_back(Nn);
    
    Alias.push_back(Nr1);Alias.push_back(Nrn);
    Alias.push_back(Nr2);Alias.push_back(Nrn);
    Alias.push_back(Nr3);Alias.push_back(Nrn);
	Reff=gsl_matrix_alloc(6,6);
	Ic=gsl_vector_alloc(6);

//	gsl_matrix_set_identity(Reff);
//	gsl_matrix_scale(Reff,1E9);	

	gsl_matrix_set(Reff,0,0,2*(Ls-Lm)/Dt);
	gsl_matrix_set(Reff,1,1,2*(Ls-Lm)/Dt);
	gsl_matrix_set(Reff,2,2,2*(Ls-Lm)/Dt);		
	gsl_matrix_set(Reff,3,3,2*(Lr-Lm)/Dt);
	gsl_matrix_set(Reff,4,4,2*(Lr-Lm)/Dt);
	gsl_matrix_set(Reff,5,5,2*(Lr-Lm)/Dt);	

	pp=p;
	
	dT=Dt;
	T=0;
	
	double Sig=1-pow(Lm,2)/(Ls*Lr);
	
	Data.A=gsl_matrix_alloc(6,6);
	gsl_matrix_set_zero(Data.A);
	gsl_matrix_set(Data.A,0,0,-Rs/(Sig*Ls));
	gsl_matrix_set(Data.A,1,1,-Rs/(Sig*Ls));
	gsl_matrix_set(Data.A,2,2,-Rr/(Sig*Lr));
	gsl_matrix_set(Data.A,3,3,-Rr/(Sig*Lr));
	gsl_matrix_set(Data.A,0,2,Rs*Lm/(Sig*Ls*Lr));
	gsl_matrix_set(Data.A,1,3,Rs*Lm/(Sig*Ls*Lr));
	gsl_matrix_set(Data.A,2,0,Rr*Lm/(Sig*Ls*Lr));
	gsl_matrix_set(Data.A,3,1,Rr*Lm/(Sig*Ls*Lr));
	gsl_matrix_set(Data.A,4,4,-D/J);
	gsl_matrix_set(Data.A,5,4,1);
	
	Data.aux=(double*)malloc(3*sizeof(double));
	Data.aux[0]=(3.0/2.0)*(Lm/(Sig*Lr*Ls));
	Data.aux[1]=J;
	
	Data.V=gsl_vector_alloc(6);
	States=gsl_vector_alloc(6);
	
	Vab0_s=gsl_vector_subvector(Data.V,0,2);
	Vab0_r=gsl_vector_subvector(Data.V,2,2);
	Torque=gsl_vector_subvector(Data.V,4,1);
	
	
	Vabc_ant=gsl_vector_alloc(6);
	
	Phi=gsl_vector_subvector(States,0,4);
	W=gsl_vector_subvector(States,4,1);
	Theta=gsl_vector_subvector(States,5,1);
	
	Iab0=gsl_vector_alloc(4);
	Temp=gsl_vector_alloc(2);
	
	Iab0_s=gsl_vector_subvector(Iab0,0,2);
	Iab0_r=gsl_vector_subvector(Iab0,2,2);
	
	Clarke=gsl_matrix_alloc(2,3);
	gsl_matrix_set(Clarke,0,0,1.0);
	gsl_matrix_set(Clarke,0,1,-0.5);
	gsl_matrix_set(Clarke,0,2,-0.5);
	gsl_matrix_set(Clarke,1,0,0);
	gsl_matrix_set(Clarke,1,1,sqrt(3.0/4.0));
	gsl_matrix_set(Clarke,1,2,-sqrt(3.0/4.0));
	
	ClarkeI=gsl_matrix_alloc(3,2);
	gsl_matrix_set(ClarkeI,0,0,1.0);
	gsl_matrix_set(ClarkeI,0,1,0.0);
	gsl_matrix_set(ClarkeI,1,0,-0.5);
	gsl_matrix_set(ClarkeI,1,1,sqrt(3.0/4.0));
	gsl_matrix_set(ClarkeI,2,0,-0.5);
	gsl_matrix_set(ClarkeI,2,1,-sqrt(3.0/4.0));
	
	MI=gsl_matrix_alloc(2,2);
	M=gsl_matrix_alloc(2,2);
	
	Phi_I=gsl_matrix_alloc(4,4);
	gsl_matrix_set(Phi_I,0,0,-Lr/(Ls*Lr-pow(Lm,2)));
	gsl_matrix_set(Phi_I,0,2,Lm/(Ls*Lr-pow(Lm,2)));
	gsl_matrix_set(Phi_I,1,1,-Lr/(Ls*Lr-pow(Lm,2)));
	gsl_matrix_set(Phi_I,1,3,Lm/(Ls*Lr-pow(Lm,2)));
	gsl_matrix_set(Phi_I,2,2,Ls/(Ls*Lr-pow(Lm,2)));
	gsl_matrix_set(Phi_I,2,0,-Lm/(Ls*Lr-pow(Lm,2)));
	gsl_matrix_set(Phi_I,3,3,Ls/(Ls*Lr-pow(Lm,2)));
	gsl_matrix_set(Phi_I,3,1,-Lm/(Ls*Lr-pow(Lm,2)));
	
	sys.function=func;
	sys.jacobian=NULL;
	sys.dimension=6;
	sys.params=&Data;
	
	d=gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45,1e-3, 1e-9, 1e-6);
	}
	
bool DFIG_AB0::Set_Views(gsl_matrix_view GprV,gsl_vector_view Vv,gsl_vector_view IhV){
	View_Gpr=GprV;
	View_V_Pri=Vv;
	View_I_Hist_Pri=IhV;
	return true;
	}
	
bool DFIG_AB0::Compute_Ih(bool euler){
	double Angle;
	Angle=gsl_vector_get(&Theta.vector,0);
	gsl_matrix_set(M,0,0,cos(Angle));	
	gsl_matrix_set(M,0,1,-sin(Angle));
	gsl_matrix_set(M,1,0,sin(Angle));
	gsl_matrix_set(M,1,1,cos(Angle));
	
	Iabc_s=gsl_vector_subvector(&View_I_Hist_Pri.vector,0,3);
	Iabc_r=gsl_vector_subvector(&View_I_Hist_Pri.vector,3,3);
	Vabc_s=gsl_vector_subvector(&View_V_Pri.vector,0,3);
	Vabc_r=gsl_vector_subvector(&View_V_Pri.vector,3,3);
	
	gsl_blas_dgemv(CblasNoTrans, 1.0, Clarke, &Vabc_s.vector, 0.0,&Vab0_s.vector);
	gsl_blas_dgemv(CblasNoTrans, 1.0, Clarke, &Vabc_r.vector, 0.0,Temp);
	gsl_blas_dgemv(CblasNoTrans, 1.0, M, Temp, 0.0,&Vab0_r.vector);
	int status=gsl_odeiv2_driver_apply(d,&T,T+dT,States->data);
	//int status=gsl_odeiv2_driver_apply_fixed_step(d,&T,dT,1,States->data);
	if(status != GSL_SUCCESS){
		printf ("error, return value=%d\n", status);
		return false;
		}
	Angle=gsl_vector_get(&Theta.vector,0);
	gsl_matrix_set(MI,0,0,cos(Angle));	
	gsl_matrix_set(MI,0,1,sin(Angle));
	gsl_matrix_set(MI,1,0,-sin(Angle));
	gsl_matrix_set(MI,1,1,cos(Angle));
	
	gsl_blas_dgemv(CblasNoTrans, 1.0, Phi_I, &Phi.vector, 0.0,Iab0);
	gsl_blas_dgemv (CblasNoTrans, 1.0,&View_Gpr.matrix, Vabc_ant, 0.0, &View_I_Hist_Pri.vector); //Compesate
	gsl_blas_dgemv(CblasNoTrans, 1.0, ClarkeI, &Iab0_s.vector, 1.0,&Iabc_s.vector);
	gsl_blas_dgemv(CblasNoTrans, 1.0, MI, &Iab0_r.vector, 0.0,Temp);
	gsl_blas_dgemv(CblasNoTrans, 1.0, ClarkeI, Temp, 1.0,&Iabc_r.vector);
	gsl_vector_memcpy(Vabc_ant,&View_V_Pri.vector);
	return true;
	}


double DFIG_AB0::Get_Speed(){
	return gsl_vector_get(&W.vector,0)/pp;
	}
void DFIG_AB0::Set_Torque(double Tm){
	gsl_vector_set(&Torque.vector,0,Tm);
	}
	
double DFIG_AB0::Get_Torque(){
	return Data.aux[2];
	}
void DFIG_AB0::Set_Speed(double Wm){
	gsl_vector_set(&W.vector,0,Wm*pp);
	}
	
double DFIG_AB0::Get_Angle(){
	return gsl_vector_get(&Theta.vector,0);
	}
	
int DFIG_AB0::func(double t, const double y[], double f[],void *params){
	struct InductionMachineData Dados= *(struct InductionMachineData *)params;
	double Te,Tm;
	gsl_vector_view F=gsl_vector_view_array(f,6);
	gsl_vector_const_view Y=gsl_vector_const_view_array(y,6);
	
	gsl_matrix_set(Dados.A,2,3,-y[4]);
	gsl_matrix_set(Dados.A,3,2,y[4]);
	Te=Dados.aux[0]*cimag((y[2]-y[3]*I)*(y[0]+y[1]*I));
	Tm=gsl_vector_get(Dados.V,4);
	
	gsl_vector_memcpy(&F.vector,Dados.V);
	gsl_vector_set(&F.vector,4,(Te-Tm)/Dados.aux[1]);
	gsl_blas_dgemv(CblasNoTrans, 1.0, Dados.A, &Y.vector, 1.0,&F.vector); //F=A*Y
	Dados.aux[2]=Te;
	return GSL_SUCCESS;
	}
	
Synchronous_Machine::Synchronous_Machine(string Na, string Nb, string Nc, string Nn, string NF1,string NF2, double L, double Dt){
    Alias.push_back(Na); Alias.push_back(Nn);
    Alias.push_back(Nb); Alias.push_back(Nn); 
    Alias.push_back(Nc); Alias.push_back(Nn);
    Alias.push_back(NF1); Alias.push_back(NF2);
	
	Reff=gsl_matrix_alloc(4,4);
	Ic=gsl_vector_alloc(4);
	for(unsigned k=0;k<3;k++)
		gsl_matrix_set(Reff,k,k,1E9); //Usar valor para estabilizar Estator
	gsl_matrix_set(Reff,3,3,1E9); //Usar valor para estabilizar Campo
}

