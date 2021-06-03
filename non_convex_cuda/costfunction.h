#ifndef COSTFUNCTION_H
#define COSTFUNCTION_H
#include<cuda_runtime.h>
#include<cuda_runtime_api.h>
#include <cublas_v2.h>
class LogCost
{
public:
    __device__ virtual  double grad(double* w, double* t, double* rho) {
        double g = -1.0 / (*w) + (*rho) * (*w - *t);
        return g;
    }
    __device__ virtual  double hess(double* w, double* t, double* rho) {
        double h = 1.0 / ((*w) * (*w)) + (*rho);
        return h;
    }
};
class Cost
{
public:
    __device__ virtual  double cost(double *w,double *w_last,double *t,double *rho,double *C);
    __device__ virtual  double grad(double *w,double *w_last,double *t,double *rho,double *C);
};
class Ss_cost: public Cost
{
public:

    __device__ virtual double cost(double *w,double *w_last,double *t,double *rho,double *C) override{
        double temp = pow((*C * ( *w - *w_last )+142.1)/1004.0,2.0);
        double f = 0.01477 * exp(-temp) * max(*w-*w_last,0.0) +
                *rho/2.0 * (*w-*t)*(*w-*t);
        return f;
    }
    __device__ virtual double grad(double *w,double *w_last,double *t,double *rho,double *C) override{
        double temp = pow((*C * ( *w - *w_last )+142.1)/1004.0,2.0);
        double g = 0.01477  * exp(-temp) * max(*w-*w_last,0.0)
                * (-2.0) * ((*C * ( *w - *w_last )+142.1)/1004.0) * (*C/1004.0) +
                *rho  * (*w-*t) ;
        return g;
    }
};
class Sb_cost: public Cost
{
public:
    __device__ virtual double cost(double *w,double *w_last,double *t,double *rho,double *C) override{
        double temp = pow((*C * ( *w - *w_last )+1621.0)/1627.0, 2.0);
        double f = 0.02079 * exp(-temp) * max(*w-*w_last,0.0) +
                *rho/2.0 * (*w-*t)*(*w-*t);
        return f;
    }
    __device__ virtual double grad(double *w,double *w_last,double *t,double *rho,double *C) override{
        double temp = pow((*C * ( *w - *w_last )+1621.0)/1627.0, 2.0);
        double g = 0.02079 * exp(-temp) * max(*w-*w_last,0.0)
                * (-2.0) * ((*C * ( *w - *w_last )+1621.0)/1627.0) * (*C/1627.0)+
                *rho  * (*w-*t) ;
        return g;
    }
};
class Sm_cost: public Cost
{
public:
    __device__ virtual double cost(double *w,double *w_last,double *t,double *rho,double *C) override{
        double temp = pow((*C * ( *w - *w_last )+198.9)/648.0, 2.0);
        double f = 0.02079 * exp(-temp) * max(*w-*w_last,0.0) +
                *rho/2.0 * (*w-*t)*(*w-*t);
        return f;
    }
    __device__ virtual double grad(double *w,double *w_last,double *t,double *rho,double *C) override{
        double temp = pow((*C * ( *w - *w_last )+198.9)/648.0, 2.0);
        double g = 0.02079 * exp(-temp) * max(*w-*w_last,0.0)
                * (-2.0) * ((*C * ( *w - *w_last )+198.9)/648.0) * (*C/648.0)+
                *rho  * (*w-*t) ;
        return g;
    }
};
#endif // COSTFUNCTION_H
