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

#include "cable.h"

Conductor::Conductor(double ri, double ro, double r, double mu){
	rho=r;
	mu_r=mu;
	Rin=ri;
	Rout=ro;
	}
	
double Conductor::Get_Rout(){
	return Rout;
	}
	
double Conductor::Get_Rin(){
	return Rin;
	}
	
double Conductor::Get_mu(){
	return mu_r;
	}
double Conductor::Get_rho(){
	return rho;
	}
	
double _Complex Conductor::Inside_Impedance(double Freq){
	double _Complex m,Z,D, mr;
	m=csqrt(I*2*M_PI*Freq*mu_0*mu_r/rho);
	//mr=Rout*m;
	if(Rin>0){
		D=besseli(1,m*Rout)*besselk(1,m*Rin)-besseli(1,m*Rin)*besselk(1,m*Rout);
		Z=((m*rho)/(D*2*M_PI*Rin))*(besseli(0,m*Rin)*besselk(1,m*Rout)+besselk(0,m*Rin)*besseli(1,m*Rout));
		}
	else{
		Z=((m*rho)/(2*M_PI*Rout))*(besseli(0,m*Rout)/besseli(1,m*Rout));
		}
	//std::cout<<Freq<<" "<<Rout<<" "<<creal(m)<<","<<cimag(m)<<" "<<mr<<" "<<besseli(0,m*Rout)<<" "<<besseli(1,m*Rout)<<" "<<Z<<std::endl;
	return Z;
	}
	
double _Complex Conductor::Outside_Impedance(double Freq){
	double _Complex m,Z,D;
	m=csqrt(I*2*M_PI*Freq*mu_0*mu_r/rho);
	if(Rin>0){
		D=besseli(1,m*Rout)*besselk(1,m*Rin)-besseli(1,m*Rin)*besselk(1,m*Rout);
		Z=((m*rho)/(D*2*M_PI*Rout))*(besseli(0,m*Rout)*besselk(1,m*Rin)+besselk(0,m*Rout)*besseli(1,m*Rin));
	}
	else
		Z=0+0.0*I;
	return Z;
	}
	
double _Complex Conductor::Mutual_Impedance(double Freq){
	double _Complex m,Z,D;
	m=csqrt(I*2*M_PI*Freq*mu_0*mu_r/rho);
	D=besseli(1,m*Rout)*besselk(1,m*Rin)-besseli(1,m*Rin)*besselk(1,m*Rout);
	Z=rho/(D*2*M_PI*Rin*Rout);	
	return Z;
	}
			
Insulation::Insulation(double e, double mu){
	eps=e;
	mu_r=mu;
	}
	
double Insulation::Get_Rout(){
	return Rout;
	}
	
double Insulation::Get_Rin(){
	return Rin;
	}
	
bool Insulation::Set_Rout(double r){
	Rout=r;
	return true;
	}
	
bool Insulation::Set_Rin(double r){
	Rin=r;
	return true;
	}
	
double Insulation::Get_mu(){
	return mu_r;
	}
double Insulation::Get_eps(){
	return eps;
	}
		
double _Complex Insulation::Impedance(double Freq){
	return I*(2*M_PI*Freq*mu_0*mu_r)/(2*M_PI)*log(Rout/Rin);
	}

double Insulation::Potential_Coeff(double Freq){
	return log(Rout/Rin)/(2*M_PI*eps_0*eps);
	}

Cable::Cable(double r){
	R=r;
	}
	
double Cable::Get_R(){
	return R;
	}
	
bool Bundle::Join(Conductor c,double x, double y){
	vCond.push_back(c);
	Pos_x.push_back(x);
	Pos_y.push_back(y);
	return true;
	}

double Bundle::Get_RMG(){
	double rmg=1;
	for(unsigned k=0;k<vCond.size();k++){
		rmg*=vCond[k].Get_Rout();
	for(unsigned j=0;j<vCond.size();j++){
		if(k!=j)
			rmg*=sqrt(pow(Pos_x[k]-Pos_x[j],2)+pow(Pos_y[k]-Pos_y[j],2));
		}
	}
	return pow(rmg,1.0/pow(vCond.size(),2));
	}
	
bool Cable::Join(Insulation c){
	vIso.push_back(c);
	return true;
	}

bool Cable::Join(Conductor c){
	vCond.push_back(c);
	return true;
	}

	
bool Cable::Reorder_Conductors(){
	for(unsigned k=0;k<vCond.size();k++){
		for(unsigned j=0;j<vCond.size();j++){
 			if(vCond[k].Get_Rout()<vCond[j].Get_Rout()){        
				Conductor c=vCond[j];   
				vCond[j]=vCond[k];
				vCond[k]=c;
				}
			}
		}
	return true;
	}
	
bool Cable::Assembly(){
	if(!Reorder_Conductors())
		return false;
	for(unsigned k=0;k<vIso.size();k++){	
		vIso[k].Set_Rin(vCond[k].Get_Rout());
		if(k<(vCond.size()-1))
			vIso[k].Set_Rout(vCond[k+1].Get_Rin());
		else
			vIso[k].Set_Rout(R);
		}
	return true;
	}

unsigned  Bundle::Get_N_Conductors(){
	return vCond.size();
	}
	
Conductor Bundle::Get_Conductor(unsigned N){
		return vCond[N];
	}
	
Insulation Cable::Get_Insulation(unsigned N){
	return vIso[N];
	}
	
void Cable::Print(){
	for(unsigned k=0;k<vCond.size();k++){
		std::cout<<"Conductor #"<<k<<std::endl;
		std::cout<<"\t Internal Radius[m]:"<<vCond[k].Get_Rin()<<" External Radius[m]:"<<vCond[k].Get_Rout()<<std::endl;
		std::cout<<"\t Resisitivity[ohm.m]:"<<vCond[k].Get_rho()<<" Relative Permissivity[H/m]:"<<vCond[k].Get_mu()<<std::endl;
		}
	for(unsigned k=0;k<vIso.size();k++){
		std::cout<<"Insulation #"<<k<<std::endl;
		std::cout<<"\t Internal Radius[m]:"<<vIso[k].Get_Rin()<<" External Radius[m]:"<<vIso[k].Get_Rout()<<std::endl;
		std::cout<<"\t Permeability[F/m]:"<<vIso[k].Get_eps()<<" Relative Permissivity[H/m]:"<<vIso[k].Get_mu()<<std::endl;
		}
	}
	
void Cable::Print(double Freq){
	double _Complex zi,zo,zm;
	double p;
	for(unsigned k=0;k<vCond.size();k++){
		std::cout<<"Conductor #"<<k<<std::endl;
		std::cout<<"\t Internal Radius[m]:"<<vCond[k].Get_Rin()<<" External Radius[m]:"<<vCond[k].Get_Rout()<<std::endl;
		std::cout<<"\t Resisitivity[ohm.m]:"<<vCond[k].Get_rho()<<" Relative Permissivity[H/m]:"<<vCond[k].Get_mu()<<std::endl;
		zi=vCond[k].Inside_Impedance(Freq);
		zo=vCond[k].Outside_Impedance(Freq);
		zm=vCond[k].Mutual_Impedance(Freq);
		std::cout<<"\t Inside Impedance[ohm/m]:"<<creal(zi)<<" "<<cimag(zi)<<"*i Outside Impedance[ohm/m]:"<<creal(zo)<<" "<<cimag(zo)<<"*i Mutual Impedance[ohm/m]:"<<creal(zm)<<" "<<cimag(zm)<<"*i"<<std::endl;
		}
	for(unsigned k=0;k<vIso.size();k++){
		std::cout<<"Insulation #"<<k<<std::endl;
		std::cout<<"\t Internal Radius[m]:"<<vIso[k].Get_Rin()<<" External Radius[m]:"<<vIso[k].Get_Rout()<<std::endl;
		std::cout<<"\t Permeability[F/m]:"<<vIso[k].Get_eps()<<" Relative Permissivity[H/m]:"<<vIso[k].Get_mu()<<std::endl;
		zi=vIso[k].Impedance(Freq);
		p=vIso[k].Potential_Coeff(Freq);
		std::cout<<"\t Impedance[ohm/m]:"<<creal(zi)<<" "<<cimag(zi)<<"*i Potencial Coeeficient[1/F]:"<<p<<std::endl;
		}
	}
CableSet::CableSet(double r,double u, unsigned Inst){
	rho=r;
	mu_r=u;
	Type=BERGERON;
	Installation=Inst;
	}

bool CableSet::Join(Cable c,double x, double y){
	c.Assembly();
	vCable.push_back(c);
	Pos_x.push_back(x);
	Pos_y.push_back(y);
	}


double CableSet::Distance(unsigned C1,unsigned C2){
	if((C1<vCable.size())&&(C2<vCable.size()))
		return sqrt(pow(Pos_x[C1]-Pos_x[C2],2)+pow(Pos_y[C1]-Pos_y[C2],2));
	return 0;
	}

double _Complex CableSet::Earth_Impedance(double Freq,unsigned C1){
	if(C1<vCable.size()){
		double _Complex z,m;
		m=csqrt(I*2*M_PI*Freq*mu_0*mu_r/rho);
		z=(rho*cpow(m,2)/(2*M_PI))*(-clog(gamma*m*vCable[C1].Get_R()/2.0)+0.5-m*Pos_y[C1]*4.0/3.0);
		//z=((rho*m)/(2*M_PI*r))*(besselk(0,m*r)/besselk(1,m*r));
		return z;
		}
	return 0+0*I;
	}

double _Complex CableSet::Earth_Mutual_Impedance(double Freq,unsigned C1, unsigned C2){
	if((C1<vCable.size())&&(C2<vCable.size())){
		double _Complex z,m;
		m=csqrt(I*2*M_PI*Freq*mu_0*mu_r/rho);
		z=(rho*cpow(m,2)/(2*M_PI))*(-clog(gamma*m*Distance(C1,C2)/2.0)+0.5-m*(Pos_y[C1]+Pos_y[C2])*2.0/3.0);
		//z=(rho*besselk(0,m*d12))/(2*M_PI*h1*h2*besselk(1,m*h1)*besselk(1,m*h2));
		return z;
		}
	return 0+0*I;
	}
bool CableSet::Compute_Parameters(double Freq){
	N=0;
	for(unsigned k=0;k<vCable.size();k++)
		N+=vCable[k].Get_N_Conductors();
		
	Z=gsl_matrix_complex_alloc(N,N);
	Y=gsl_matrix_complex_alloc(N,N);
	Gamma=gsl_matrix_complex_alloc(N,N);

	Zm=gsl_matrix_complex_alloc(N,N);
	Ym=gsl_matrix_complex_alloc(N,N);
	
	Zcf=gsl_matrix_alloc(N,N);
	Zcm=gsl_matrix_complex_alloc(N,N);

	Rf=gsl_matrix_alloc(N,N);
	Rm=gsl_matrix_alloc(N,N);

	
	Lambda=gsl_vector_complex_alloc(N);
	Tv_cplx=gsl_matrix_complex_alloc(N,N);
	Ti_cplx=gsl_matrix_complex_alloc(N,N);
	Tv=gsl_matrix_alloc(N,N);
	Ti=gsl_matrix_alloc(N,N);
	Ti_inv=gsl_matrix_alloc(N,N);
	Tv_inv=gsl_matrix_alloc(N,N);

	Compute_Z(Freq);
	Compute_Y(Freq);
	Modal_Trans(Freq);
	return true;
	}
bool CableSet::Compute_Y(double Freq){
	double Pp,Pm;
	int Nc=0,Noff=0;
	Insulation iso;
	for(unsigned k=0;k<vCable.size();k++){
		Pp=0;
		Pm=0;
		Nc=vCable[k].Get_N_Conductors()-1;
		for(int n=Nc;n>=0;n--){
			iso=vCable[k].Get_Insulation(n);
			Pp+=iso.Potential_Coeff(Freq);
			matrix_set(Y,n+Noff,n+Noff,Pp);
			//Y[N*(n+Noff)+n+Noff]=Pp;
			for(int l=0;l<n;l++){
				matrix_set(Y,l+Noff,n+Noff,Pp);
				matrix_set(Y,n+Noff,l+Noff,Pp);
				//Y[N*(l+Noff)+n+Noff]=Pp;
				//Y[N*(n+Noff)+l+Noff]=Pp;
				}
			//if(Installation==AIR){
			//	}
			}
		Noff+=Nc+1;
		}
	Matrix_Inverse(Y);
	Matrix_Scale(Y,2*M_PI*Freq*I);

	return true;
	}
bool CableSet::Compute_Z(double Freq){
	double _Complex Zp,Zmt;
	int Nc=0,Noff=0;
	Conductor c_in,c_out;
	Insulation i_out;
	for(unsigned k=0;k<vCable.size();k++){
		Zp=0;
		Zmt=0;
		Nc=vCable[k].Get_N_Conductors()-1;
		for(int n=Nc;n>=0;n--){
			c_in=vCable[k].Get_Conductor(n);
			i_out=vCable[k].Get_Insulation(n);
			Zp+=c_in.Outside_Impedance(Freq)+i_out.Impedance(Freq);
			if(n<Nc){ //not last conductor
				c_out=vCable[k].Get_Conductor(n+1);
				Zp+=c_out.Inside_Impedance(Freq)-2*c_out.Mutual_Impedance(Freq);
				}
			matrix_set(Z,n+Noff,n+Noff,Zp);
			//Z[N*(n+Noff)+n+Noff]=Zp;
			}
		for(int n=Nc;n>=1;n--){
			c_in=vCable[k].Get_Conductor(n);
			i_out=vCable[k].Get_Insulation(n);
			Zmt+=c_in.Outside_Impedance(Freq)+i_out.Impedance(Freq)-c_in.Mutual_Impedance(Freq);
			if(n<(Nc)){ //not last conductor
				c_out=vCable[k].Get_Conductor(n+1);
				Zmt+=c_out.Inside_Impedance(Freq)-c_out.Mutual_Impedance(Freq);
				}
			for(int l=0;l<n;l++){
				matrix_set(Z,l+Noff,n+Noff,Zmt);
				matrix_set(Z,n+Noff,l+Noff,Zmt);
//					Z[N*(l+Noff)+n+Noff]=Zmt;
//					Z[N*(n+Noff)+l+Noff]=Zmt;
				}
			}
		Noff+=Nc+1;
		}
unsigned NoffCol=0,NoffLin=0;
	for(unsigned k=0;k<vCable.size();k++){
			for(int m=0;m<vCable.size();m++){
				if(m==k){
					Zp=Earth_Impedance(Freq,k);
					}
				else{
					Zp=Earth_Mutual_Impedance(Freq,k,m);
					}
			for(int c1=0;c1<vCable[k].Get_N_Conductors();c1++){
				for(int c2=0;c2<vCable[m].Get_N_Conductors();c2++){
					matrix_acc(Z,c1+NoffLin,c2+NoffCol,Zp);
					//Z[N*(c1+NoffLin)+c2+NoffCol]+=Zp;
					}
				}
			NoffCol+=vCable[m].Get_N_Conductors();
			}
		NoffCol=0;
		NoffLin+=vCable[k].Get_N_Conductors();
		}
	return true;
	}
	
	
bool CableSet::Modal_Trans(double Freq){

//matrix_set(Z,0,0,1.2569087E-04+1.4086245E-03*I);
//matrix_set(Z,0,1,1.0721643E-04+1.3292159E-03*I);
//matrix_set(Z,0,2,1.0144081E-04+1.2475746E-03*I);
//matrix_set(Z,1,0,1.0721643E-04+1.3292159E-03*I);
//matrix_set(Z,1,1,4.0780345E-04+1.3278125E-03*I);
//matrix_set(Z,1,2,1.0144081E-04+1.2475746E-03*I);
//matrix_set(Z,2,0,1.0144081E-04+1.2475746E-03*I);
//matrix_set(Z,2,1,1.0144081E-04+1.2475746E-03*I);
//matrix_set(Z,2,2,2.1228107E-04+1.2271508E-03*I);

	gsl_matrix_complex *P;
	P=gsl_matrix_complex_alloc(N,N);
	Matrix_Product(Z,Y,P);
    Eigen(P, Tv_cplx, Lambda);
//    printf("\nTv\n");
//    Matrix_Print(Tv_cplx);
//   	Vector_Print(Lambda);
//    gsl_eigen_nonsymmv_sort(Lambda,Tv_cplx,GSL_EIGEN_SORT_ABS_ASC);
//    printf("\nTv\n");
//    Matrix_Print(Tv_cplx);
//	Vector_Print(Lambda);
	
   	Matrix_Product(Y,Z,P);
    Eigen(P, Ti_cplx, Lambda);
   	printf("\nTi\n");
	Matrix_Print(Ti_cplx);
//	Vector_Print(Lambda);
//    gsl_eigen_nonsymmv_sort(Lambda,Ti_cplx,GSL_EIGEN_SORT_ABS_ASC);
	printf("\nTi\n");
	Matrix_Print(Ti_cplx);
//	Vector_Print(Lambda);
                   

//    printf("\nTv\n");
//    Matrix_Print(Tv);

	Matrix_Copy_Real(Tv_cplx, Tv);
	Matrix_Copy_Real(Ti_cplx, Ti);

	Matrix_Inverse(Tv_inv,Tv);
//	printf("\nTv_i\n");
//	Matrix_Print(Tv_inv);
	
	Matrix_Inverse(Ti_inv,Ti);
//	printf("\nTi_i\n");
//	Matrix_Print(Ti_inv);
	
	

//            
//	matrix_set(Tv,0,0,1.0120724+5.3212771*I);
//    matrix_set(Tv,0,1,.2425654+8.8083966*I);
//    matrix_set(Tv,0,2,1.0000000);
//    matrix_set(Tv,1,0,1.0285814+5.1042885*I);
//    matrix_set(Tv,1,1,1.0000000);
//    matrix_set(Tv,1,2,.1509504-8.9793113*I);
//    matrix_set(Tv,2,0,1.0000000);
//    matrix_set(Tv,2,1,.0699121-102.9568563*I);
//    matrix_set(Tv,2,2,.0100537-96.3878629*I);
//    
//    matrix_set(Ti,0,0,.0029093+171.8526780*I);
//    matrix_set(Ti,0,1,.0694309+72.4406461*I);
//    matrix_set(Ti,0,2,.9900077-4.0475933*I);
//    matrix_set(Ti,1,0,.1540592+170.8862574*I);
//    matrix_set(Ti,1,1,1.0278826-3.4713315*I);
//    matrix_set(Ti,1,2,.9021608-177.4307344*I);
//    matrix_set(Ti,2,0,1.0403098-.0189639*I);
//    matrix_set(Ti,2,1,.2787170-160.7558572*I);
//    matrix_set(Ti,2,2,.7874568+178.4792275*I);
//	Matrix_Inverse(Tv_inv,Tv);
//	Matrix_Inverse(Ti_inv,Ti);
//	
//	Matrix_Transform(Ti_inv,Y,Tv,Ym);
//	
//	
//Matrix_Product(gsl_complex_rect(1.0,0.0),CblasTrans,Tv,CblasNoTrans,Y,gsl_complex_rect(0.0,0.0),P);
//Matrix_Product(gsl_complex_rect(1.0,0.0),CblasNoTrans,P,CblasNoTrans,Tv,gsl_complex_rect(0.0,0.0),Ym);

	Matrix_Transform(Tv_inv,Z,Ti,Zm);
	Matrix_Transform(Ti_inv,Y,Tv,Ym);
		
//			printf("Zm\n");
//	Matrix_Print(Zm);
//	printf("Ym\n");
//	Matrix_Print(Ym);

	for(unsigned k=0;k<N;k++){
		Vel.push_back((2*M_PI*Freq)/cimag(csqrt(vector_get(Lambda,k))));
		matrix_set(Gamma,k, k, (csqrt((matrix_get(Zm,k,k))*matrix_get(Ym,k,k))));
		matrix_set(Zcm,k, k, (csqrt((matrix_get(Zm,k,k))/matrix_get(Ym,k,k))));
		matrix_set(Rm,k,k,creal(matrix_get(Zm,k,k)));
		}
	//Matrix_Transform(Tv,Zcm,Ti_inv,Zcf);
	Matrix_Transform(Tv,Rm,Ti_inv,Rf);
	gsl_matrix_complex_free(P);
	return true;
	}

void CableSet::Print(){
	printf("Z\n");
	Matrix_Print(Z);
	printf("Y\n");
	Matrix_Print(Y);
	printf("Tv\n");
	Matrix_Print(Tv);
	printf("Ti\n");
	Matrix_Print(Ti);
	printf("Zm\n");
	Matrix_Print(Zm);
	printf("Ym\n");
	Matrix_Print(Ym);
	printf("Zcm\n");
	Matrix_Print(Zcm);
	printf("Gamma\n");
	Matrix_Print(Gamma);
	printf("\n");
//	printf("Zcf\n");
//	Matrix_Print(Zcf);
//	printf("\n");
//	printf("Velocity\n");
//	for(unsigned k=0;k<N;k++){
//		printf("%E ",Vel[k]);
//		}	
//	printf("\n");
	}
	
gsl_matrix* CableSet::Get_Ti(){
	return Ti;
	}
gsl_matrix* CableSet::Get_Tv(){
	return Tv;
	}
gsl_matrix* CableSet::Get_Ti_Inverse(){
	return Ti_inv;
	}
gsl_matrix* CableSet::Get_Tv_Inverse(){
	return Tv_inv;
	}
			

gsl_matrix_complex* CableSet::Get_Zc_Mode(){
	return Zcm;
	}
gsl_matrix* CableSet::Get_Zc_Phase(){
	return Zcf;
	}
	
gsl_matrix* CableSet::Get_R_Mode(){
	return Rm;
	}
gsl_matrix* CableSet::Get_R_Phase(){
	return Rf;
	}
	
	double CableSet::Get_Speed(unsigned Mode){
	if(Mode<Vel.size())
		return Vel[Mode];
	return 0;
	}
	
double _Complex CableSet::Get_Zc_Mode(unsigned Mode){
	if((Mode<Zcm->size1)&&(Mode<Zcm->size2))
		return g2c(gsl_matrix_complex_get(Zcm,Mode,Mode));
	return 0;
	}

double CableSet::Get_R_Mode(unsigned Mode){
	if((Mode<Rm->size1)&&(Mode<Rm->size2))
		return gsl_matrix_get(Rm,Mode,Mode);
	return 0;
	}
		
double CableSet::Get_Zc_Phase(unsigned Mode){
	if((Mode<Zcf->size1)&&(Mode<Zcf->size2))
		return gsl_matrix_get(Zcf,Mode,Mode);
	return 0;
	}

double CableSet::Get_R_Phase(unsigned Mode){
	if((Mode<Rf->size1)&&(Mode<Rf->size2))
		return gsl_matrix_get(Rf,Mode,Mode);
	return 0;
	}
	
unsigned CableSet::Get_N_Modes(){
	return N;
	}
	
CableSet::~CableSet(){
	if(N>0){
	//free(Pi);free(P0);free(Zi);free(Z0);
	free(Y);free(Z);free(Ym);free(Zm);
	free(Lambda);free(Tv);free(Tv_inv);free(Ti);free(Ti_inv);
	}
 }
	

OverHeadLineSet::OverHeadLineSet(double r,double u){
	rho=r;
	mu_r=u;
	Type=BERGERON;
	Compute_abcd(carson_terms);
	}

bool OverHeadLineSet::Join(Conductor c,double x, double y, double sag){
	Bundle B1;
	B1.Join(c,x,y);
	vBundle.push_back(B1);
	Pos_x.push_back(x);
	Pos_y.push_back(y+sag/3.0);
	}

bool OverHeadLineSet::Join(Bundle b,double x, double y, double sag){
	vBundle.push_back(b);
	Pos_x.push_back(x);
	Pos_y.push_back(y+sag/3.0);
	}

double OverHeadLineSet::Angle(unsigned C1,unsigned C2){
	if((C1<vConductor.size())&&(C2<vConductor.size()))
		return asin((Pos_x[C2]-Pos_x[C1])/Distance_Image(C1,C2));
	return 0;
	}

double OverHeadLineSet::Distance_Image(unsigned C1,unsigned Ci2){
	if((C1<vConductor.size())&&(Ci2<vConductor.size()))
		return sqrt(pow(Pos_x[C1]-Pos_x[Ci2],2)+pow(Pos_y[C1]+Pos_y[Ci2],2));
	return 0;
	}

double OverHeadLineSet::Distance(unsigned C1,unsigned Ci2){
	if((C1<vConductor.size())&&(Ci2<vConductor.size()))
		return sqrt(pow(Pos_x[C1]-Pos_x[Ci2],2)+pow(Pos_y[C1]-Pos_y[Ci2],2));
	return 0;
	}

double _Complex OverHeadLineSet::External_Impedance(double Freq,unsigned C1,unsigned C2){
	if((C1<vConductor.size())&&(C2<vConductor.size())){
		double _Complex Z;
		if(C1==C2)
			Z=I*Freq*mu_0*log(Distance_Image(C1,C2)/vConductor[C1].Get_Rout());
		else
			Z=I*Freq*mu_0*log(Distance_Image(C1,C2)/Distance(C1,C2));
		return Z;
		}
	return 0+0*I;
	}
double _Complex OverHeadLineSet::Earth_Impedance(double Freq,unsigned C1,unsigned C2){
	if((C1<vConductor.size())&&(C2<vConductor.size())){
		double _Complex z,m;
		double a, dR,dX, theta,D;
		a=4*M_PI*sqrt(5)*1E-4*Distance_Image(C1,C2)*sqrt(Freq/rho);
		theta=Angle(C1,C2);
		dR=M_PI/8;
		dX=0.5*(0.6159315-log(a));
		for(unsigned k=1;k<carson_terms;k+=4){
		dR+=-b[k]*pow(a,k)*cos(theta*k);
		dR+=b[k+1]*((c[k+1]-log(a))*pow(a,k+1)*cos((k+1)*theta)+theta*pow(a,k+1)*sin((k+1)*theta));
		dR+=b[k+2]*pow(a,k+2)*cos(theta*(k+2));
		dR-=d[k+3]*pow(a,k+3)*cos(theta*(k+3));
		
		dX+=b[k]*pow(a,k)*cos(theta*k);
		dX-=d[k+1]*pow(a,k+1)*cos(theta*(k+1));
		dX-=b[k+2]*pow(a,k+2)*cos(theta*(k+2));
		dX-=b[k+3]*((c[k+3]-log(a))*pow(a,k+3)*cos((k+3)*theta)+theta*pow(a,k+3)*sin((k+3)*theta));

		}
		dR*=4*2*M_PI*Freq*1E-7;
		dX*=4*2*M_PI*Freq*1E-7;
		//z=((rho*m)/(2*M_PI*r))*(besselk(0,m*r)/besselk(1,m*r));
		return dR+dX*I;;
		}
	return 0+0*I;
	}


bool OverHeadLineSet::Compute_abcd(int terms){
	b.resize(terms+1);
	d.resize(terms+1);
	c.resize(terms+1);
	b[1]=(sqrt(2.0)/6);
	b[2]=(1.0/16);
	c[2]=1.3659315;
	d[1]=b[1]*M_PI/4;
	d[2]=b[2]*M_PI/4;
	double sig=1;
	for(unsigned k=3;k<=terms;k++){
		if((k%4)==1)
			sig=-sig;
		b[k]=sig*fabs(b[k-2])*1.0/(k*(k+2));
		}
	for(unsigned k=0;k<=terms;k++){
		d[k]=b[k]*M_PI/4;
		}
	for(unsigned k=4;k<=terms;k++)
		c[k]=c[k-2]+1.0/k+1.0/(k+2);
	for(unsigned k=0;k<=terms;k++){
		std::cout<<"B["<<k<<"]="<<b[k]<<" C["<<k<<"]="<<c[k]<<" D["<<k<<"]="<<d[k]<<std::endl;
		}
	return true;
	}

OverHeadLineSet::~OverHeadLineSet(){
	if(N>0){
	//free(Pi);free(P0);free(Zi);free(Z0);
	free(Y);free(Z);free(Ym);free(Zm);
	free(Lambda);free(Tv);free(Tv_inv);free(Ti);free(Ti_inv);
	}
}


bool OverHeadLineSet::Compute_Parameters(double Freq){
	N=vConductor.size();
		
	Z=gsl_matrix_complex_alloc(N,N);
	Gamma=gsl_matrix_complex_alloc(N,N);
	Y=gsl_matrix_complex_alloc(N,N);
;

	Zm=gsl_matrix_complex_alloc(N,N);
	Ym=gsl_matrix_complex_alloc(N,N);
	
	Zcf=gsl_matrix_alloc(N,N);
	Zcm=gsl_matrix_complex_alloc(N,N);

	Rf=gsl_matrix_alloc(N,N);
	Rm=gsl_matrix_alloc(N,N);

	
	Lambda=gsl_vector_complex_alloc(N);
	Tv_cplx=gsl_matrix_complex_alloc(N,N);
	Ti_cplx=gsl_matrix_complex_alloc(N,N);
	Tv=gsl_matrix_alloc(N,N);
	Ti=gsl_matrix_alloc(N,N);
	Ti_inv=gsl_matrix_alloc(N,N);
	Tv_inv=gsl_matrix_alloc(N,N);

	Compute_Z(Freq);
	Compute_Y(Freq);
	Modal_Trans(Freq);
	return true;
	}


bool OverHeadLineSet::Compute_Z(double Freq){
	double _Complex Zint,Zext,Zearth;
	for(unsigned k=0;k<N;k++)
	for(unsigned j=k;j<N;j++){
		if(k==j)
			Zint=vConductor[k].Inside_Impedance(Freq);
		else
			Zint=0+0*I;
		Zearth=Earth_Impedance(Freq,k,j);
		Zext=External_Impedance(Freq,k,j);
		//std::cout<<k<<" "<<j<<" "<<creal(Zext)<<","<<cimag(Zext)<<"\t"<<creal(Zint)<<","<<cimag(Zint)<<"\t"<<creal(Zearth)<<","<<cimag(Zearth)<<std::endl;
		matrix_set(Z,k,j,Zext+Zint+Zearth);
		matrix_set(Z,j,k,Zext+Zint+Zearth);
		}
	return true;
	}


bool OverHeadLineSet::Compute_Y(double Freq){
	double _Complex p;
	for(unsigned k=0;k<N;k++)
		for(unsigned j=k;j<N;j++){
			if(k==j)
				p=log(Distance_Image(k,j)/vConductor[k].Get_Rout());	
			else
				p=log(Distance_Image(k,j)/Distance(k,j));
			p/=2*M_PI*eps_0;
			matrix_set(Y,k,j,p);
			matrix_set(Y,j,k,p);
			}
	Matrix_Inverse(Y);
	Matrix_Scale(Y,2*M_PI*Freq*I);
	return true;
	}

bool OverHeadLineSet::Compute_Zequ(double D){
	gsl_matrix_complex *Temp1,*Temp2,*Temp3;
	N=vConductor.size();
	Theta1=gsl_matrix_complex_alloc(N,N);
	Theta2=gsl_matrix_complex_alloc(N,N);
	Theta3=gsl_matrix_complex_alloc(N,N);

	Temp1=gsl_matrix_complex_alloc(N,N);
	Temp2=gsl_matrix_complex_alloc(N,N);
	Temp3=gsl_matrix_complex_alloc(N,N);

	A=gsl_matrix_complex_alloc(N,N);
	B=gsl_matrix_complex_alloc(N,N);
	for(unsigned k=0;k<N;k++){
		matrix_set(Theta1,k,k,csinh(matrix_get(Gamma,k,k)*D/matrix_get(Zcm,k,k)));
		matrix_set(Theta2,k,k,ccosh(matrix_get(Gamma,k,k)*D));
		matrix_set(Theta3,k,k,matrix_get(Zcm,k,k)*csinh(matrix_get(Gamma,k,k)*D));
		}
	/*std::cout<<"-----Theta1-----"<<std::endl;
	for(unsigned k=0;k<N;k++){
		for(unsigned j=0;j<N;j++){
			double _Complex Imp=matrix_get(Theta1,k,j);
			std::cout<<creal(Imp)<<","<<cimag(Imp)<<"\t";
			}
		std::cout<<std::endl;
		}
	std::cout<<"-----Theta2-----"<<std::endl;
	for(unsigned k=0;k<N;k++){
		for(unsigned j=0;j<N;j++){
			double _Complex Imp=matrix_get(Theta2,k,j);
			std::cout<<creal(Imp)<<","<<cimag(Imp)<<"\t";
			}
		std::cout<<std::endl;
		}
	std::cout<<"-----Theta3-----"<<std::endl;
	for(unsigned k=0;k<N;k++){
		for(unsigned j=0;j<N;j++){
			double _Complex Imp=matrix_get(Theta3,k,j);
			std::cout<<creal(Imp)<<","<<cimag(Imp)<<"\t";
			}
		std::cout<<std::endl;
		}*/
	Matrix_Inverse(Theta3);
	Matrix_Transform(Ti_cplx,Theta2,Theta3,Temp1);
	Matrix_Product(c2g(1.0+0*I), CblasNoTrans, Temp1, CblasTrans, Ti_cplx,c2g(0.0+0*I) ,A);

	std::cout<<"-----A-----"<<std::endl;
	for(unsigned k=0;k<N;k++){
		for(unsigned j=0;j<N;j++){
			double _Complex Imp=matrix_get(A,k,j);
			std::cout<<creal(Imp)<<","<<cimag(Imp)<<"\t";
			}
		std::cout<<std::endl;
		}
	double _Complex Sa=0+0*I;
	for(unsigned k=0;k<N;k++)
		for(unsigned j=0;j<N;j++)
			Sa+=matrix_get(A,k,j);
	std::cout<<"-----Sa="<<creal(Sa)<<","<<cimag(Sa)<<std::endl;

	Matrix_Transform(Ti_cplx,Theta1,Temp2,B); 
	Matrix_Product(Temp1,Theta2,Temp3);
	Matrix_Product(c2g(-1.0+0*I), CblasNoTrans,Temp3, CblasConjTrans, Ti_cplx,c2g(-1.0+0*I) ,B);
	

	std::cout<<"-----B-----"<<std::endl;
	for(unsigned k=0;k<N;k++){
		for(unsigned j=0;j<N;j++){
			double _Complex Imp=matrix_get(B,k,j);
			std::cout<<creal(Imp)<<","<<cimag(Imp)<<"\t";
			}
		std::cout<<std::endl;
		}
	double _Complex Sb=0+0*I;
	for(unsigned k=0;k<N;k++)
		for(unsigned j=0;j<N;j++)
			Sb+=matrix_get(B,k,j);
	std::cout<<"-----Sb="<<creal(Sb)<<","<<cimag(Sb)<<std::endl;
	
	double _Complex ZcEq,GammaEq;
	D/=1000;
	GammaEq=(1.0/D)*cacosh(-Sa/Sb);
	ZcEq=-1.0/(Sb*csinh(GammaEq*D));

	std::cout<<"-----GammaEq="<<creal(GammaEq)<<","<<cimag(GammaEq)<<std::endl;
	std::cout<<"-----ZcEq="<<creal(ZcEq)<<","<<cimag(ZcEq)<<std::endl;
	return true;
	}
/*
bool OverHeadLineSet::Print(){
	std::cout<<"-----Z-----"<<std::endl;
	for(unsigned k=0;k<N;k++){
		for(unsigned j=0;j<N;j++){
			double _Complex Imp=matrix_get(Z,k,j);
			std::cout<<creal(Imp)<<","<<cimag(Imp)<<"\t";
			}
		std::cout<<std::endl;
		}
	std::cout<<"-----Y-----"<<std::endl;
	for(unsigned k=0;k<N;k++){
		for(unsigned j=0;j<N;j++){
			double _Complex Imp=matrix_get(Y,k,j);
			std::cout<<creal(Imp)<<","<<cimag(Imp)<<"\t";
			}
		std::cout<<std::endl;
	}


	std::cout<<"-----Zm-----"<<std::endl;
	for(unsigned k=0;k<N;k++){
		for(unsigned j=0;j<N;j++){
			double _Complex Imp=matrix_get(Zm,k,j);
			std::cout<<creal(Imp)<<","<<cimag(Imp)<<"\t";
			}
		std::cout<<std::endl;
	}
	std::cout<<"-----Ym-----"<<std::endl;
	for(unsigned k=0;k<N;k++){
		for(unsigned j=0;j<N;j++){
			double _Complex Imp=matrix_get(Ym,k,j);
			std::cout<<creal(Imp)<<","<<cimag(Imp)<<"\t";
			}
		std::cout<<std::endl;
	}
	return true;
	}*/


