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

#ifndef SUPPORT_H
#define SUPPORT_H

#include <complex.h>
#include <math.h>
#include <cstdio>
#include<cstdlib>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include<gsl/gsl_vector.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>

double _Complex g2c(gsl_complex T);
gsl_complex c2g(double _Complex a);

void vector_set(gsl_vector_complex *V,unsigned n, double _Complex a);
void vector_set(gsl_vector *V,unsigned n, double a);

void vector_acc(gsl_vector_complex *V,unsigned n, double _Complex a);

double _Complex vector_get(gsl_vector_complex *V,unsigned n);
double vector_get(gsl_vector *V,unsigned n);


void matrix_set(gsl_matrix_complex *V,unsigned l, unsigned c, double _Complex a);
void matrix_set(gsl_matrix *V,unsigned l, unsigned c, double a);

void matrix_acc(gsl_matrix_complex *V,unsigned l, unsigned c, double _Complex a);

void matrix_div(gsl_matrix_complex *V,unsigned l, unsigned c, double _Complex a);

double _Complex matrix_get(gsl_matrix_complex *V,unsigned l, unsigned c);
double matrix_get(gsl_matrix *V,unsigned l, unsigned c);


extern "C" void zbesk_(double*, double*, double*, int*, int*, double*, double*, int*, int*);
extern "C" void zbesi_(double*, double*, double*, int*, int*, double*, double*, int*, int*);
extern "C" void zgeev_(char *jobvl, char *jobvr, int *n, 
	double _Complex *a, int *lda, double _Complex *w, double _Complex *vl, 
	int *ldvl, double _Complex *vr, int *ldvr, double _Complex *work, 
	int *lwork, double *rwork, int *info);

double _Complex besselk(double nu, double _Complex z);
gsl_complex besselk(double nu, gsl_complex z);
double _Complex besseli(double nu, double _Complex z);
gsl_complex besseli(double nu, gsl_complex z);
//int Eigen(int N, double _Complex *A, double _Complex *vec, double _Complex *val);
int Eigen(gsl_matrix_complex *A, gsl_matrix_complex *vec, gsl_vector_complex *val);

bool Matrix_Inverse(gsl_matrix_complex *A);
bool Matrix_Inverse(gsl_matrix *A);

bool Matrix_Inverse(gsl_matrix_complex *A, gsl_matrix_complex *B);
bool Matrix_Inverse(gsl_matrix *A,gsl_matrix *B);

bool Matrix_Product(gsl_complex alpha, CBLAS_TRANSPOSE_t tA, gsl_matrix_complex *A, CBLAS_TRANSPOSE_t tB, gsl_matrix_complex *B, gsl_complex beta ,gsl_matrix_complex *C );
bool Matrix_Product(double alpha,  CBLAS_TRANSPOSE_t tA, gsl_matrix *A, CBLAS_TRANSPOSE_t tB, gsl_matrix *B, double beta ,gsl_matrix *C );

bool Matrix_Product(gsl_matrix_complex *A, gsl_matrix_complex *B, gsl_matrix_complex *C );
bool Matrix_Product(gsl_matrix *A, gsl_matrix *B,gsl_matrix *C);

bool Matrix_Copy(gsl_matrix_complex *A, gsl_matrix_complex *B);
bool Matrix_Copy(gsl_matrix *A, gsl_matrix *B);
bool Matrix_Copy_Real(gsl_matrix *A, gsl_matrix_complex *B);
bool Matrix_Copy_Real(gsl_matrix_complex *A, gsl_matrix *B);

bool Matrix_Transpose(gsl_matrix_complex *A, gsl_matrix_complex *B);
bool Matrix_Transpose(gsl_matrix *A, gsl_matrix *B);
bool Matrix_Transpose(gsl_matrix_complex *A);
bool Matrix_Transpose(gsl_matrix *A);

bool Matrix_Transform(gsl_matrix_complex *T1,gsl_matrix_complex *A,gsl_matrix_complex *T2, gsl_matrix_complex *M);
bool Matrix_Transform(gsl_matrix *T1,gsl_matrix *A,gsl_matrix *T2, gsl_matrix *M);
bool Matrix_Transform(gsl_matrix *M1,gsl_matrix_complex *A,gsl_matrix *M2,gsl_matrix_complex *M);

bool Matrix_Scale(gsl_matrix_complex *A, gsl_complex s);
bool Matrix_Scale(gsl_matrix_complex *A, double _Complex s);
bool Matrix_Scale(gsl_matrix *A, double s);

bool Matrix_Copy_Real(gsl_matrix_complex *A, gsl_matrix *B);

void Matrix_Print(gsl_matrix *M);
void Matrix_Print(gsl_matrix_complex *M);

void Vector_Print(gsl_vector *M);
void Vector_Print(gsl_vector_complex *M);
#endif
