#ifndef PORTFOLIOSOLVER_H
#define PORTFOLIOSOLVER_H
#include"portfoliomodel.h"
#define _RECORD_
using namespace Eigen;
using namespace Eigen::internal;
using namespace Eigen::Architecture;
#include <cuda_runtime.h>
#include<cublas.h>
#include <cusolverDn.h>
#include<parameters.h>
void matrix_times(cublasHandle_t handle,double *f,double *g,int mf,int nf, int kg,double *res);




class SolverRecord
{
public:
    std::vector<double> errors;
    std::vector<double> values;
    int maxrecord;
    SolverRecord( );
};
class portfoliosolver
{
public:
    portfoliosolver();
    int RegisterModel(PortfolioModel *m){ this->model = m; return 0;}

    int SolveModel(const Option &);
    int init();
    SolverRecord record;
    double result[10000];
private:
    MatrixXd P;
    double* P_cu_;
    double* Pt_cu_;
    double* A_cu_;
    double* APt_cu_;
    double* PAt_cu_;
    MatrixXd Pt;
    MatrixXd APt;
    MatrixXd PAt;
    double* S;
    PortfolioModel *model;

};

#endif // PORTFOLIOSOLVER_H
