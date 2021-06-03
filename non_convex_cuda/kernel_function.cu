#include<cuda_runtime_api.h>
#include<handler_cuda_error.h>
#include<Common.h>
#include <cublas_v2.h>
#include<stdio.h>
#include"kenel_function.h"
#include<math.h>
const double gamma0 = 1.5;
const double gmh = 1.9 ;
const double gmg = 1.1 ;
const double minval = 1e-16;
const double orthval = 0.2;

double curv_adaptive_BB(double al_h,double de_h){

    double tmph = de_h/al_h;
    double tau_h;
    if (tmph > .8)
        tau_h = de_h;
    else
        tau_h = al_h - 0.8*de_h;



    return tau_h;

}
void aradmm_estimate(cublasHandle_t handle,
                     double* rho,
                               double *gamma_,
                               double *Au,
                               double *Au0,
                               double *l_hat,
                               double *l_hat0,
                               double *Bv,
                               double *Bv0,
                               double *l,
                               double * l0,int n){

    //Au0 = -Au + Au0
    cublasDaxpy(handle,n,&alpha_1,Au,1,Au0,1);
    cublasDaxpy(handle,n,&alpha_1,Bv,1,Bv0,1);
    //lhat0 = - lhat + lhat0
    cublasDaxpy(handle,n,&alpha_1,l_hat,1,l_hat0,1);
    cublasDaxpy(handle,n,&alpha_1,l ,1,l0,1);

    double ul_hat = 0.0;
    cublasDdot (  handle,   n,Au0, 1,
                               l_hat0, 1,&ul_hat);

    double vl = 0.0;
    cublasDdot (  handle,   n,Bv0, 1,
                               l0, 1,&vl);
    //cublasDasum(handle, n, Au0 ,1, &ul_hat);
    double dl = 0.0;
    double dl_hat = 0.0;
    cublasDnrm2(  handle,   n, l_hat0, 1, &dl_hat);
    cublasDnrm2(  handle,   n, l0, 1, &dl);

    double du = 0.0;
    double dv = 0.0;
    cublasDnrm2(  handle,   n, Au0, 1, &du);
    cublasDnrm2(  handle,   n, Bv0, 1, &dv);


//    printf("ul_hat=%.4f\n",ul_hat);
//    printf("vl=%.4f\n",vl);
//    printf("dl=%.4f\n",dl);
//    printf("dlhat=%.4f\n",dl_hat);
//    printf("du=%.4f\n",du);
//    printf("dv=%.4f\n",dv);
    bool hflag = false;
    bool gflag = false;
    double tau = 0.0;
    double bb_h,bb_g;
    if (ul_hat > orthval*du*dl_hat + minval){
        hflag = true;
        double al_h = dl_hat*dl_hat/ul_hat;
        double de_h = ul_hat/(du*du);
        bb_h = curv_adaptive_BB(al_h, de_h);
    }
    if (vl > orthval*dv*dl + minval){
        gflag = true;
        double al_g = dl*dl/vl;
        double de_g = vl/(dv*dv);
        bb_g = curv_adaptive_BB(al_g, de_g);
    }
    double gamma;
    if (hflag && gflag){
        double ss_h = sqrt(bb_h);
        double ss_g = sqrt(bb_g);
        gamma = std::min(1 + 2.0/(ss_h/ss_g + ss_g/ss_h), gamma0);
        tau = ss_h*ss_g;
    }else if (hflag){
        gamma = gmh;
        tau = bb_h;
    }
    else if (gflag){
        gamma = gmg;
        tau = bb_g;
    }
    else
        gamma = gamma0;

    if(hflag || gflag){
        cudaMemcpy(rho,&tau,sizeof(double),cudaMemcpyHostToDevice);
        cudaMemcpy(rho,&tau,sizeof(double),cudaMemcpyHostToDevice);
    }


}


__global__ void update_f(double*f,double*R,double * z,double *u,double *rho,int n){
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    if (x>=n)
        return;
    f[x] = -R[x] - *rho * z[x] + u[x];
}

__global__ void update_u_new(double *u_new, double * u_cu,double *rho,double *x_hat_cu,
                             double *z_cu,int n){
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    if (x>=n)
        return;
    u_new[x] =u_cu[x] + (*rho) *(x_hat_cu[x] - z_cu[x]);


}
__global__ void update_x_hat(double * xhat, double*x,double *zold,double *gamma,int n){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i>=n)
        return;
    xhat[i] = *gamma * x[i] + (1 - *gamma) * zold[i];
}
__global__ void update_tt_cu(double * tt, double *x_hat, double *u,double *rho,int n){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i>=n)
        return;
    tt[i] =  x_hat[i] +  u[i]/ (*rho);
}
__global__ void update_z_cu(double * z, double *tt, double * lb, double *ub,
                            int n){
    // z =   min(max(tt ,lb),ub) ;
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i>=n)
        return;
    z[i] =  fmax(tt[i],lb[i]);
    z[i] =  fmin(z[i],ub[i]);

}
__global__ void generate_invs(double *f,double *d,double *rho,int n){
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    if (x>=n)
        return;
    f[x*n+x] = 1.0/(d[x] + *rho - 1.0);
}
__global__ void update_l0_hat0_and_so_on(double * l0_cu, double *lhat0_cu,
                         double *Bv0_cu,double *Au0_cu,
                                        double *u_new,double *u_cu,
                         double *x_hat_cu,double *z_old,
                                         double * z_cu,double *rho,int n){

    int x = blockIdx.x * blockDim.x + threadIdx.x;
    if (x>=n)
        return;
    l0_cu[x] = u_new[x];
    lhat0_cu[x] = u_cu[x] + (*rho) * (x_hat_cu[x] - z_old[x]);
    Bv0_cu[x] =  z_cu[x];
    Au0_cu[x] =  -x_hat_cu[x];
}

__global__ void update_hat_and_so_on(double *lhat0_cu,
                         double *Bv0_cu,double *Au0_cu,
                                         double *u_cu,
                         double *x_hat_cu,double *z_old,
                                         double * z_cu,double *rho,int n){

    int x = blockIdx.x * blockDim.x + threadIdx.x;
    if (x>=n)
        return;

    lhat0_cu[x] = u_cu[x] + (*rho) * (x_hat_cu[x] - z_old[x]);
    Bv0_cu[x] =  z_cu[x];
    Au0_cu[x] =  -x_hat_cu[x];
}
__global__ void prepare_for_norm(double *norm_cu,double *norm_cu2,double *x_cu,
                                 double *z_cu,double *rho,
                                 double *z_old,int n){


    int x = blockIdx.x * blockDim.x + threadIdx.x;
    if (x>=n)
        return;
    norm_cu2[x]  = x_cu[x] - z_cu[x];
    norm_cu[x]  = *rho *(z_old[x] - z_cu[x]);

}
