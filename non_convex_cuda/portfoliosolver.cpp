#include "portfoliosolver.h"
#include<Common.h>
#include<Eigen/SVD>
#include <cuda_runtime.h>
#include<cublas.h>
#include <cusolverDn.h>
 #include<handler_cuda_error.h>
using namespace Eigen;
using namespace Eigen::internal;
using namespace Eigen::Architecture;

int portfoliosolver::init(){
    int n = this->model->n;
    MatrixXd I = MatrixXd::Identity(n, n);
    MatrixXd H = this->model->Sigma * 2 + I;



    int  lwork = 0;


    int *devInfo = NULL;double *d_work = NULL;
    cusolverDnHandle_t cusolverH = NULL;
    CHECK(cudaMalloc ((void**)&this->Pt_cu_, sizeof(double) * n * n));
    CHECK(cudaMalloc ((void**)&this->P_cu_, sizeof(double) * n * n));
    CHECK(cudaMalloc ((void**)&this->S, sizeof(double) * n));
    CHECK(cudaMalloc ((void**)&devInfo, sizeof(int)));
    CHECK(cudaMemcpy(this->Pt_cu_, H.data(), sizeof(double) * n * n, cudaMemcpyHostToDevice));
    cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR; // compute eigenvalues and eigenvectors.
    cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;
    cusolverStatus_t cusolver_status = CUSOLVER_STATUS_SUCCESS;
    cusolver_status = cusolverDnCreate(&cusolverH);
    assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
    cudaError_t cudaStat1 = cudaSuccess;
    cudaError_t cudaStat2 = cudaSuccess;

    cusolver_status = cusolverDnDsyevd_bufferSize(
        cusolverH,
        jobz,
        uplo,
        n,
        this->Pt_cu_,
        n,
        this->S,
        &lwork);

    assert (cusolver_status == CUSOLVER_STATUS_SUCCESS);
    cudaStat1 = cudaMalloc((void**)&d_work, sizeof(double)*lwork);
    assert(cudaSuccess == cudaStat1);
    cusolver_status = cusolverDnDsyevd(
        cusolverH,
        jobz,
        uplo,
        n,
        this->Pt_cu_,
        n,
        this->S,
        d_work,
        lwork,
        devInfo);
    cudaStat1 = cudaDeviceSynchronize();
    assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
    assert(cudaSuccess == cudaStat1);

    //cudaStat2 = cudaMemcpy(U.data(), this->Pt_cu_, sizeof(double)*n*n, cudaMemcpyDeviceToHost);

    // this->P_cu_
    double const alpha(1.0);
    double const beta(0.0);
    cublasHandle_t handle;
    cublasCreate(&handle);
    cublasDgeam(handle, CUBLAS_OP_T, CUBLAS_OP_N, n, n, &alpha, this->Pt_cu_, n, &beta,
                this->Pt_cu_, n, this->P_cu_, n);

    assert(cudaSuccess == cudaStat2);

    int m = this->model->m;
    CHECK(cudaMalloc((void**)&this->A_cu_, sizeof(double) * m * n));
    CHECK(cudaMemcpy(this->A_cu_, this->model->A.data(), sizeof(double) * m * n,
                     cudaMemcpyHostToDevice));
    // Acu * Pt_cu
    CHECK(cudaMalloc((void**)&this->APt_cu_, sizeof(double) * m * n));
    CHECK(cudaMalloc((void**)&this->PAt_cu_, sizeof(double) * m * n));
    matrix_times(handle, this->A_cu_, this->Pt_cu_, m, n, n, this->APt_cu_);


    cublasDgeam(handle, CUBLAS_OP_T, CUBLAS_OP_T, n, m, &alpha, this->APt_cu_, m, &beta,
                this->APt_cu_, m, this->PAt_cu_, n);

}
SolverRecord::SolverRecord( ):maxrecord(50000){
    this->errors.reserve(maxrecord);
    this->values.reserve(maxrecord);
}
