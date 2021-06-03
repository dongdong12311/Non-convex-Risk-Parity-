#ifndef DFP_CU_H
#define DFP_CU_H
#include<cuda_runtime.h>
#include"account.h"
#include"costfunction.h"
__global__   void dfp_choose(double *res,int indss_s,int indss_e,
                                      int indsb_s,int indsb_e,
                                      int indsm_s, int indsm_e,
                                      double *tt,
                                      double *rho,
                                      double *lb,double *ub,
                                      Account * account,int n);
// choose three types of functions to solve the minumum value

__global__ void my_newton(double* res,
    double* tt,
    double* rho, int n);

__global__ void calcost(double *res,double *x_cu,int indss_s,int indss_e,
                           int indsb_s,int indsb_e,
                           int indsm_s, int indsm_e,
                           Account * account,int n);
__global__ void callogcost(double *res,double *x_cu,int n);
#endif // DFP_CU_H
