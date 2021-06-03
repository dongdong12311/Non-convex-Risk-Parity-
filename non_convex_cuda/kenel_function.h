#ifndef KENEL_FUNCTION_H
#define KENEL_FUNCTION_H
#include<cuda_runtime_api.h>
#include<handler_cuda_error.h>
#include<Common.h>
#include <cublas_v2.h>
#include<stdio.h>

static const double alpha_1   =  -1.0 ;
static const double alpha1  =  1.0 ;
static const double beta   =  0.0 ;


__global__ void update_f(double*f,double*R,double * z,double *u,double *rho,int n);



__global__ void update_x_hat(double * xhat, double*x, double *zold, double *gamma, int n);


__global__ void update_tt_cu(double * tt, double *x_hat, double *u,double *rho,int n);



__global__ void update_z_cu(double * z, double *tt, double * lb, double *ub,
                            int n);

__global__ void update_u_new(double *u_new, double * u_cu, double *rho, double *x_hat_cu, double *z_cu, int n);

__global__ void generate_invs(double *f,double *d,double *rho,int n);

__global__ void update_l0_hat0_and_so_on(double * l0_cu, double *lhat0_cu,
                         double *Bv0_cu,double *Au0_cu,
                                        double *u_new,double *u_cu,
                         double *x_hat_cu,double *z_old,
                                        double * z_cu,double *rho,int n);

__global__ void update_hat_and_so_on(double *lhat0_cu,
                         double *Bv0_cu,double *Au0_cu,
                                         double *u_cu,
                         double *x_hat_cu,double *z_old,
                                         double * z_cu,double *rho,int n);

__global__ void prepare_for_norm(double *norm_cu, double *norm_cu2, double *x_cu,
                                 double *z_cu, double *rho,
                                 double *z_old, int n);

void aradmm_estimate(cublasHandle_t handle,
                     double* tau,
                               double *gamma,
                               double *Au,
                               double *Au0,
                               double *l_hat,
                               double *l_hat0,
                               double *Bv,
                               double *Bv0,
                               double *l,
                               double * l0,int n);
#endif // KENEL_FUNCTION_H
