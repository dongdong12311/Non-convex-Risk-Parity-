#include "portfoliomodel.h"
#include<stdlib.h>
#include<iostream>
PortfolioModel::PortfolioModel(int m, int n)
{
    this->m = m;
    this->n = n;
    Sigma = MatrixXd(n,n);
    A  = MatrixXd(m,n);
    b  = MatrixXd(m,1); // m * 1
    Er = MatrixXd(n,1); // n * 1
    lb = MatrixXd(n,1); // n * 1
    ub = MatrixXd(n,1); // n * 1
}
int PortfolioModel::RegisterIndex(std::vector<int> &sb,std::vector<int> &sm,std::vector<int> &ss){
    this->n_Sb = sb;
    this->n_Sm = sm;
    this->n_Ss = ss;
    return 0;

}
int PortfolioModel::RegisterParameters(double C,double *w_last){
    this->C = C;
    this->w_last = w_last;
}

int PortfolioModel::RegisterData(std::vector<double> & Sigma,std::vector<double> & A,
                  std::vector<double> &b,
                  std::vector<double> &Er,std::vector<double> &lb,
                                  std::vector<double> &ub){

    if(this->n * this->n != Sigma.size()){
        printf("error Sigma size");

        return -1;
    }
    if(this->n * this->m != A.size()){
        printf("error A size");
        return -1;
    }

    for (int i = 0 ; i< this->n;++i){
        this->lb(i,0) = lb[i];
        this->ub(i,0) = ub[i];
        this->Er(i,0) = Er[i];
        for(int j = 0;j<n;++j){
            this->Sigma(i,j) = Sigma[i*n + j];
        }
    }
    for (int i = 0 ; i< m;++i){
        this->b(i,0) = b[i];
        for(int j = 0;j<n;++j){
            this->A(i,j) = A[i*n + j];
        }
    }

    return 0;


}
