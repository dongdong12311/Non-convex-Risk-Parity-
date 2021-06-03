#ifndef PORTFOLIOMODEL_H
#define PORTFOLIOMODEL_H
#include<iostream>
#include <Eigen/Dense>
#include<vector>

using Eigen::MatrixXd;
using Eigen::MatrixXi;


class PortfolioModel
{
public:
    PortfolioModel(int m,int n);
    MatrixXd Sigma; // sigma n * n
    MatrixXd A; // m * n
    MatrixXd b; // m * 1
    MatrixXd Er; // n * 1
    MatrixXd lb; // n * 1
    MatrixXd ub; // n * 1
    int m;
    int n;
    std::vector<int> n_Ss;
    std::vector<int> n_Sb;
    std::vector<int> n_Sm;
    double C;
    double *w_last;
public:

    int RegisterData(std::vector<double> & Sigma,std::vector<double> & A,
                      std::vector<double> &b,
                      std::vector<double> &Er,std::vector<double> &lb,
                      std::vector<double> &ub);
    int RegisterIndex(std::vector<int> & ,std::vector<int> & ,std::vector<int> & );
    int RegisterParameters(double C,double *w_last);

};


#endif // PORTFOLIOMODEL_H
