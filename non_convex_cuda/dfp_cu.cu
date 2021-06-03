#include<cuda_runtime_api.h>
#include<handler_cuda_error.h>
#include<Common.h>
#include <cublas_v2.h>
#include<stdio.h>
#include"dfp_cu.h"

__device__ double newton_solver(LogCost * cost, double *tt,double *rho,
                      double *x){
                      
            double er = 1;
            double tol=1e-12;
            int iter=0;
            double a;
            double xn;
            double grad, hess;
            while (er>tol){
                
                grad = cost->grad(x,tt,rho);
                hess = cost->hess(x,tt,rho);
                a = grad/hess;
                xn= *x - a ;
                er = abs( a / *x);
                *x=xn;
                iter++;
            }
            return *x; 
                      
 }
__global__ void my_newton(double *res,
                           double *tt,
                           double *rho,int n){

    int x = blockIdx.x * blockDim.x + threadIdx.x;
    double x0 = 1e-5;
    if (x>=n)
        return;
    LogCost cost;
    res[x] =  newton_solver( &cost, &tt[x], rho, &x0);

    return ;

}
__device__ double dfp_solver(Cost * cost, double *tt,double *rho,
                      double *x0,double *w_last,double *C){
    int maxk=100;

    double rho0=0.55;
    const double sigma=0.4;
    double epsilon=1e-10;
    int k=0;
    double Hk=1.0;
    double x= *x0;
    double gk = 0.0;
    double dk = 0.0;

    while(k < maxk){

        gk = cost->grad(x0,w_last,tt,rho,C);
//        printf("gk = %.4f\n",gk);
//        printf(" k = %d\n",k);

        if(  (gk*gk) < epsilon )
            break;
        dk = -Hk * gk;
        int m=0;
        int mk=0;
        while( m < 15){
            double p = pow(rho0,m);
            double temp = (*x0) + p * dk;
            double temp1 = cost->cost(&temp,w_last,tt,rho,C);
            double temp2 = cost->cost(x0,w_last,tt,rho,C) + sigma * p * gk * dk;
            if( temp1 < temp2){
                mk=m;break;
        }
            m++;
        }
        double p = pow(rho0,mk);
        x = *x0 + p * dk;

        double sk= x - *x0;
        double yk=cost->grad(&x,w_last,tt,rho,C) - gk;
        Hk = sk / (yk);
        k++;
        *x0 = x;

    }

    return x;
}
__global__ void dfp_choose(double *res,int indss_s,int indss_e,
                           int indsb_s,int indsb_e,
                           int indsm_s, int indsm_e,
                           double *tt,
                           double *rho,
                           double *lb,double *ub,
                           Account * account,int n){

    int x = blockIdx.x * blockDim.x + threadIdx.x;
    double x0 = 0.0;

    if (x>=n)
        return;
    if ((x < indss_e) && ( x>= indss_s)  ){
        Ss_cost cost_ss;
        res[x] =  dfp_solver( &cost_ss, &tt[x], rho,
                           &x0, &account->w_last[x], &account->C);
        if(res[x] < lb[x])
            res[x] = lb[x];
        if(res[x] > ub[x])
            res[x] = ub[x];
        return ;
    }
    if ( (x < indsb_e) && (x >= indsb_s)){
        Sb_cost cost_sb;
        res[x] =  dfp_solver( &cost_sb, &tt[x], rho,
                           &x0, &account->w_last[x], &account->C);
        if(res[x] < lb[x])
            res[x] = lb[x];
        if(res[x] > ub[x])
            res[x] = ub[x];
        return ;
    }
    if ( (x < indsm_e) && (x >= indsm_s)){

        Sm_cost cost_sm;
        res[x] =  dfp_solver( &cost_sm, &tt[x], rho,
                           &x0, &account->w_last[x], &account->C);
        if(res[x] < lb[x])
            res[x] = lb[x];
        if(res[x] > ub[x])
            res[x] = ub[x];
        return;
    }
    return ;





}

__global__ void callogcost(double *res,double *x_cu,int n){
    int x = blockIdx.x * blockDim.x + threadIdx.x;



    if (x>=n)
        return;
    res[x] =  log(x_cu[x]);

}
__global__ void calcost(double *res,double *x_cu,int indss_s,int indss_e,
                           int indsb_s,int indsb_e,
                           int indsm_s, int indsm_e,
                           Account * account,int n){

    int x = blockIdx.x * blockDim.x + threadIdx.x;



    if (x>=n)
        return;
    res[x] = 0.0;
    if ((x < indss_e) && ( x>= indss_s)  ){

        double temp = pow(( account->C * ( x_cu[x] -  account->w_last[x] )+142.1)/1004.0,2.0);
        res[x] = 0.01477 * exp(-temp) * max(( x_cu[x] -  account->w_last[x] ),0.0);

        return ;
    }
    if ( (x < indsb_e) && (x >= indsb_s)){
        double temp = pow(( account->C * ( x_cu[x] -  account->w_last[x] )+1621.0)/1627.0,2.0);
        res[x] = 0.02079 * exp(-temp) * max(( x_cu[x] -  account->w_last[x] ),0.0);

        return ;
    }
    if ( (x < indsm_e) && (x >= indsm_s)){

        double temp = pow(( account->C * ( x_cu[x] -  account->w_last[x] )+198.9)/648.0, 2.0);
        res[x] = 0.02079 * exp(-temp) * max(( x_cu[x] -  account->w_last[x] ),0.0);

        return;
    }
    return ;





}
