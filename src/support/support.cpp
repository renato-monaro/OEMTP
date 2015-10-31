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

#include "support.h"

double _Complex g2c(gsl_complex T){
	return GSL_REAL(T)+GSL_IMAG(T)*I; 
	}

gsl_complex c2g(double _Complex a){
	gsl_complex T;
	GSL_SET_COMPLEX(&T, creal(a),cimag(a));
	return T;
	}

void vector_set(gsl_vector_complex *V,unsigned n, double _Complex a){
	gsl_vector_complex_set(V,n,c2g(a));
	}

void vector_set(gsl_vector *V,unsigned n, double  a){
	gsl_vector_set(V,n,a);
	}
double _Complex vector_get(gsl_vector_complex *V,unsigned n){
	return g2c(gsl_vector_complex_get(V,n)); 
	}
	
double vector_get(gsl_vector *V,unsigned n){
	return gsl_vector_get(V,n); 
	}

void vector_acc(gsl_vector_complex *V,unsigned n, double _Complex a){
	gsl_complex T1;
	T1=gsl_vector_complex_get(V,n);
	gsl_vector_complex_set(V,n,gsl_complex_add(T1,c2g(a)));
	}

void matrix_set(gsl_matrix_complex *V,unsigned l, unsigned c, double _Complex a){
	gsl_matrix_complex_set(V,l,c,c2g(a));
	}
	
void matrix_set(gsl_matrix *V,unsigned l, unsigned c, double a){
	gsl_matrix_set(V,l,c,a);
	}

double _Complex matrix_get(gsl_matrix_complex *V,unsigned l, unsigned c){
	return g2c(gsl_matrix_complex_get(V,l,c));
	}
	
double matrix_get(gsl_matrix *V,unsigned l, unsigned c){
	return gsl_matrix_get(V,l,c);
	}

void matrix_acc(gsl_matrix_complex *V,unsigned l, unsigned c, double _Complex a){
	gsl_complex T1;
	T1=gsl_matrix_complex_get(V,l,c);
	gsl_matrix_complex_set(V,l,c,gsl_complex_add(T1,c2g(a)));
	}

void matrix_div(gsl_matrix_complex *V,unsigned l, unsigned c, double _Complex a){
	double _Complex T1;
	T1=matrix_get(V,l,c);
	gsl_matrix_complex_set(V,l,c,c2g(T1/a));
	}

double _Complex besselk(double nu, double _Complex z){
  double _Complex res;
  int kode=1,n=1,nz,ierr;
  double zr=creal(z);
  double zi=cimag(z);
  double cyr[1],cyi[1];

zbesk_(&zr,&zi,&nu,&kode,&n,cyr,cyi,&nz,&ierr);
  if(ierr!=0){
    printf("error!\n");
    }
  res=cyr[0]+I*cyi[0];
  return res;
}

double _Complex besseli(double nu, double _Complex z){
	double _Complex res;
  int kode=1,n=1,nz,ierr;
  double zr=creal(z);
  double zi=cimag(z);
  double cyr[1],cyi[1];
  
  zbesi_(&zr,&zi,&nu,&kode,&n,cyr,cyi,&nz,&ierr);
  if(ierr!=0){
    printf("error!\n");
  }
  res=cyr[0]+I*cyi[0];
  return res;
}

gsl_complex besseli(double nu, gsl_complex z){
	return c2g(besseli(nu,g2c(z)));
	}
gsl_complex besselk(double nu, gsl_complex z){
	return c2g(besselk(nu,g2c(z)));
	}
	
bool Matrix_Inverse(gsl_matrix *A){
	int s;
	gsl_matrix *T;
	gsl_permutation *p;
	
	T=gsl_matrix_alloc(A->size1,A->size2);
	gsl_matrix_memcpy(T, A);
	
  	p=gsl_permutation_alloc (T->size1);
  	
	gsl_linalg_LU_decomp (T, p, &s);
    gsl_linalg_LU_invert (T, p,A);
    
	gsl_permutation_free(p);
    gsl_matrix_free(T);
	return true;
	}
	
bool Matrix_Inverse(gsl_matrix *A, gsl_matrix *B){
	int s;
	gsl_matrix *T;
	gsl_permutation *p;
	
	T=gsl_matrix_alloc(B->size1,B->size2);
	gsl_matrix_memcpy(T, B);
	
  	p=gsl_permutation_alloc (T->size1);
  	
	gsl_linalg_LU_decomp (T, p, &s);
    gsl_linalg_LU_invert (T, p,A);
    
	gsl_permutation_free(p);
    gsl_matrix_free(T);
	return true;
	}
	
bool Matrix_Inverse(gsl_matrix_complex *A){
	int s;
	gsl_matrix_complex *T;
	gsl_permutation *p;
	
	T=gsl_matrix_complex_alloc(A->size1,A->size2);
	gsl_matrix_complex_memcpy(T, A);
	
  	p=gsl_permutation_alloc (T->size1);
  	
	gsl_linalg_complex_LU_decomp (T, p, &s);
    gsl_linalg_complex_LU_invert (T, p,A);
    
	gsl_permutation_free(p);
    gsl_matrix_complex_free(T);
	return true;
	}
	
bool Matrix_Inverse(gsl_matrix_complex *A, gsl_matrix_complex *B){
	int s;
	gsl_matrix_complex *T;
	gsl_permutation *p;
	
	T=gsl_matrix_complex_alloc(B->size1,B->size2);
	gsl_matrix_complex_memcpy(T, B);
	
  	p=gsl_permutation_alloc (T->size1);
  	
	gsl_linalg_complex_LU_decomp (T, p, &s);
    gsl_linalg_complex_LU_invert (T, p,A);
    
	gsl_permutation_free(p);
    gsl_matrix_complex_free(T);
	return true;
	}
	
	
	
bool Matrix_Scale(gsl_matrix_complex *A, gsl_complex s){
	gsl_matrix_complex_scale(A, s);
	return true;
	}
	
bool Matrix_Scale(gsl_matrix_complex *A, double _Complex s){
	gsl_matrix_complex_scale(A, c2g(s));
	return true;
	}
	
bool Matrix_Scale(gsl_matrix *A, double s){
	gsl_matrix_scale(A, s);
	return true;
	}
	

bool Matrix_Product(gsl_complex alpha, CBLAS_TRANSPOSE_t tA, gsl_matrix_complex *A, CBLAS_TRANSPOSE_t tB, gsl_matrix_complex *B, gsl_complex beta ,gsl_matrix_complex *C ){
	gsl_blas_zgemm(tA, tB,alpha,A,B,beta,C);
	return true;
	}
bool Matrix_Product(double alpha,  CBLAS_TRANSPOSE_t tA, gsl_matrix *A,  CBLAS_TRANSPOSE_t tB, gsl_matrix *B, double beta ,gsl_matrix *C ){
	gsl_blas_dgemm(tA, tB,alpha,A,B,beta,C);
	return true;
	}

bool Matrix_Product(gsl_matrix_complex *A, gsl_matrix_complex *B, gsl_matrix_complex *C ){
	return Matrix_Product(gsl_complex_rect(1.0,0.0),CblasNoTrans,A,CblasNoTrans,B,gsl_complex_rect(0.0,0.0),C);
	}
bool Matrix_Product(gsl_matrix *A, gsl_matrix *B,gsl_matrix *C){
	return Matrix_Product(1.0,CblasNoTrans,A,CblasNoTrans,B,0.0,C);
	}
	
bool Matrix_Transform(gsl_matrix_complex *T1,gsl_matrix_complex *A,gsl_matrix_complex *T2,gsl_matrix_complex *M){
	gsl_matrix_complex *T;
	T=gsl_matrix_complex_alloc(T1->size1,A->size2);
	Matrix_Product(T1,A,T);
	Matrix_Product(T,T2,M);
	return true;
	}
bool Matrix_Transform(gsl_matrix *T1,gsl_matrix *A,gsl_matrix *T2, gsl_matrix *M){
	gsl_matrix *T;
	T=gsl_matrix_alloc(T1->size1,A->size2);
	Matrix_Product(T1,A,T);
	Matrix_Product(T,T2,M);
	return true;
	}
	
bool Matrix_Transform(gsl_matrix *M1,gsl_matrix_complex *A,gsl_matrix *M2,gsl_matrix_complex *M){
	gsl_matrix_complex *T,*T1,*T2;
	T1=gsl_matrix_complex_alloc(M1->size1,M1->size2);
	T2=gsl_matrix_complex_alloc(M2->size1,M2->size2);
	T=gsl_matrix_complex_alloc(M1->size1,A->size2);
	Matrix_Copy_Real(M1,T1);
	Matrix_Copy_Real(M2,T2);
	Matrix_Product(T1,A,T);
	Matrix_Product(T,T2,M);
	return true;
	}
	
bool Matrix_Transpose(gsl_matrix_complex *A, gsl_matrix_complex *B){
 	gsl_matrix_complex_transpose_memcpy(A,B);
	return true;
	}
bool Matrix_Transpose(gsl_matrix *A, gsl_matrix *B){
 	gsl_matrix_transpose_memcpy(A,B);
	return true;
	}
	
bool Matrix_Copy(gsl_matrix_complex *A, gsl_matrix_complex *B){
 	gsl_matrix_complex_memcpy(A,B);
	return true;
	}
bool Matrix_Copy(gsl_matrix *A, gsl_matrix *B){
 	gsl_matrix_memcpy(A,B);
	return true;
	}
	
bool Matrix_Copy_Real(gsl_matrix *A, gsl_matrix_complex *B){
	if((A->size1==B->size1)&&(A->size2==B->size2)){
 		for(unsigned k=0;k<A->size1;k++)
 		 		for(unsigned l=0;l<A->size2;l++)
 		 			matrix_set(B,k,l,matrix_get(A,k,l)+0.0*I);
		return true;
		}
	return false;
	}
	
bool Matrix_Copy_Real(gsl_matrix_complex *A, gsl_matrix *B){
	if((A->size1==B->size1)&&(A->size2==B->size2)){
 		for(unsigned k=0;k<A->size1;k++)
 		 		for(unsigned l=0;l<A->size2;l++)
 		 			matrix_set(B,k,l,creal(matrix_get(A,k,l)));
		return true;
		}
	return false;
	}
	
	
bool Matrix_Transpose(gsl_matrix_complex *A){
 	gsl_matrix_complex_transpose(A);
	return true;
	}
	
bool Matrix_Transpose(gsl_matrix *A){
 	gsl_matrix_transpose(A);
	return true;
	}
	
	
int Eigen(gsl_matrix_complex *A, gsl_matrix_complex *vec, gsl_vector_complex *val) {

	gsl_matrix_complex *T;
	T=gsl_matrix_complex_alloc(A->size1,A->size2);
	
//	Matrix_Copy(T,A);
	Matrix_Transpose(T,A);
	char jobvl='N';
	char jobvr='V';
	int n=T->size1;
	int lda=T->tda;
	int ldvl=T->tda;
	int ldvr=T->tda;
	int info;
	int lwork=2*n;
	double *rwork;
	double _Complex *work;
	
	rwork=(double*)malloc(2*n*sizeof(double));
	work=(double _Complex*)malloc(sizeof(double _Complex)*lwork);
	zgeev_(&jobvl,&jobvr,&n, (double _Complex*)T->data, &lda, (double _Complex*)val->data,NULL, &ldvl,(double _Complex*)vec->data, &ldvr, work, &lwork, rwork, &info);
	Matrix_Transpose(vec);
//	for(int k=0;k<A->size2;k++){
//		double _Complex Scale;
//		for(int l=0;l<n;l++)
//			Scale+=cpow(cabs(matrix_get(vec,l,k)),2);
//		printf("NORM %E %E\n",creal(Scale),cimag(Scale));
//		for(int l=0;l<n;l++){
//			matrix_div(vec,l,k,csqrt(Scale));
//			}
//		}
	gsl_matrix_complex_free(T);
	return info;
	}


void Matrix_Print(gsl_matrix_complex *M){
	double _Complex T;
	for(unsigned l=0;l<M->size1;l++){
		for(unsigned k=0;k<M->size2;k++){
			T=matrix_get(M,l,k);
			printf("%E %E*i\t",creal(T),cimag(T));
			}
		printf("\n");
		}
	}
	
void Matrix_Print(gsl_matrix *M){
	double T;
	for(unsigned l=0;l<M->size1;l++){
		for(unsigned k=0;k<M->size2;k++){
			T=matrix_get(M,l,k);
			printf("%E\t",T);
			}
		printf("\n");
		}
	}
	
void Vector_Print(gsl_vector_complex *M){
	double _Complex T;
	for(unsigned l=0;l<M->size;l++){
		T=vector_get(M,l);
		printf("%E %E*i\t",creal(T),cimag(T));
		}
	printf("\n");
	}
	
void Vector_Print(gsl_vector *M){
	double T;
	for(unsigned l=0;l<M->size;l++){
		T=vector_get(M,l);
		printf("%E\t",T);
		}
	printf("\n");
	}


