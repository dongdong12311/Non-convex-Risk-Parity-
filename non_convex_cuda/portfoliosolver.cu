#include "portfoliosolver.h"
#include<cuda_runtime_api.h>
#include<handler_cuda_error.h>
#include<Common.h>
#include <cublas_v2.h>
#include<stdio.h>
#include"kenel_function.h"
#include"dfp_cu.h"
double objective_log(dim3 &block,cublasHandle_t handle,
                     double * sigma_cu, double *norm_cu,double *x_cu,int n);

portfoliosolver::portfoliosolver(){

    this->record = SolverRecord();

}


void show_res(double *s,int n){
    double *res = (double *)malloc(n*sizeof(double));
    CHECK(cudaMemcpy(res,s,n*sizeof(double),cudaMemcpyDeviceToHost));
    for (int i = 0 ; i< n;++i){
        printf("%.6f\t",res[i]);
    }

    free(res);
}
void matrix_times(cublasHandle_t handle,double *f,double *g,int mf,int nf, int kg,double *res){

    cublasDgemm(handle,CUBLAS_OP_N,CUBLAS_OP_N,mf,kg,nf,
                 &alpha1 ,f,mf,g,nf,&beta,res,mf);
}
void matrix_timesT(cublasHandle_t handle,double *f,double *g,int mf,int nf, int kg,double *res){

    cublasDgemm(handle,CUBLAS_OP_N,CUBLAS_OP_T,mf,kg,nf,
                 &alpha1 ,f,mf,g,nf,&beta,res,mf);
}
int portfoliosolver::SolveModel(const Option & opts){



    cublasHandle_t handle;
    cublasCreate(&handle);

    //prepare for choose portfolio
    Ss_cost *ss_cost; CHECK(cudaMalloc(&ss_cost, sizeof(Ss_cost)));
    Sb_cost *sb_cost; CHECK(cudaMalloc(&sb_cost, sizeof(Sb_cost)));
    Sm_cost *sm_cost; CHECK(cudaMalloc(&sm_cost, sizeof(Sm_cost)));
//    int ss_s = model->n_Ss[0];
//    int ss_e = model->n_Ss[model->n_Ss.size()-1];
//    int sb_s = model->n_Sb[0];
//    int sb_e = model->n_Sb[model->n_Sb.size()-1];
//    int sm_s = model->n_Sm[0];
//    int sm_e = model->n_Sm[model->n_Sm.size()-1];



    //
    int m = model->m;
    int n = model->n;int n_size = n * sizeof(double);

    // norm
    double eps_pri = 0.0;
    double eps_dual = 0.0;
    double r_norm = 0.0;
    double s_norm = 0.0;
    double sqrtn = sqrt(n)*opts.ABSTOL;
    double *norm_cu;CHECK(cudaMalloc(&norm_cu,n_size ));
    double *norm_cu2;CHECK(cudaMalloc(&norm_cu2,n_size ));

    double *z_cu,*z_old;
    double *u_cu;
    double *u_new;

    double *f_cu;
    double *R_cu;
    double *P_cu, *Pt_cu;
    double *Pb_cu, *x2_cu;

    double *APt_cu;
    double *PAt_cu;
    double *invs_cu;
    double *S_cu;
    double * Ainv;
    int * info ;


    double *tt_cu;
    CHECK(cudaMalloc(&tt_cu,n*sizeof(double)));

    double *b_cu;
    CHECK(cudaMalloc(&b_cu,n*sizeof(double)));

    CHECK(cudaMalloc(&info, sizeof(int)));

    CHECK(cudaMalloc(&z_cu,n*sizeof(double)));
    CHECK(cudaMalloc(&z_old,n*sizeof(double)));
    double *x_cu, *x_hat_cu;
    CHECK(cudaMalloc(&x_cu,n*sizeof(double)));
    CHECK(cudaMalloc(&x_hat_cu,n*sizeof(double)));

    CHECK(cudaMalloc(&Ainv,m*m*sizeof(double)));



    CHECK(cudaMalloc(&u_cu,n*sizeof(double)));
    CHECK(cudaMalloc(&u_new,n*sizeof(double)));
    CHECK(cudaMalloc(&f_cu,n*sizeof(double)));
    CHECK(cudaMalloc(&R_cu,n*sizeof(double)));

    CHECK(cudaMalloc(&invs_cu,n*n*sizeof(double)));
    CHECK(cudaMalloc(&Pb_cu,n*sizeof(double)));
    CHECK(cudaMalloc(&x2_cu,n*sizeof(double)));

    CHECK(cudaMalloc(&APt_cu,n*m*sizeof(double)));

    CHECK(cudaMalloc(&PAt_cu,n*m*sizeof(double)));

    CHECK(cudaMemcpy(R_cu,this->model->Er.data(),n*sizeof(double),cudaMemcpyHostToDevice));


    Pt_cu = this->Pt_cu_;
    P_cu = this->P_cu_;
    S_cu = this->S;



    APt_cu = this->APt_cu_;
    PAt_cu = this->PAt_cu_;

    CHECK(cudaMemcpy(b_cu,this->model->b.data(),n*sizeof(double),cudaMemcpyHostToDevice));


    double *sigma_cu; CHECK(cudaMalloc(&sigma_cu,n*n*sizeof(double)));
    CHECK(cudaMemcpy(sigma_cu,this->model->Sigma.data(),n*n*sizeof(double),cudaMemcpyHostToDevice));


    double *lb_cu;CHECK(cudaMalloc(&lb_cu,n*sizeof(double)));
    double *ub_cu;CHECK(cudaMalloc(&ub_cu,n*sizeof(double)));

    double *l0_cu;CHECK(cudaMalloc(&l0_cu,n_size ));
    double *lhat0_cu;CHECK(cudaMalloc(&lhat0_cu,n_size ));
    double *Bv0_cu;CHECK(cudaMalloc(&Bv0_cu,n_size ));
    double *Au0_cu;CHECK(cudaMalloc(&Au0_cu,n_size ));
    double *lhat_cu;CHECK(cudaMalloc(&lhat_cu,n_size ));
    double *Bv_cu;CHECK(cudaMalloc(&Bv_cu,n_size ));
    double *Au_cu;CHECK(cudaMalloc(&Au_cu,n_size ));



    CHECK(cudaMemcpy(lb_cu,this->model->lb.data(),n*sizeof(double),cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(ub_cu,this->model->ub.data(),n*sizeof(double),cudaMemcpyHostToDevice));

    double *temp1_cu; // m * n
    double *temp2_cu; // m * m
    CHECK(cudaMalloc(&temp1_cu,m * n*sizeof(double)));
    CHECK(cudaMalloc(&temp2_cu, m * m*sizeof(double)));

    int freq = opts.adp_freq;
    int siter = std::max(opts.adp_start_iter-1, 1);
    int eiter = std::min(opts.adp_end_iter, opts.max_iter)+1;


    double *gamma_;CHECK(cudaMalloc(&gamma_, sizeof(double)));
    cudaMemcpy(gamma_, &opts.gamma,  sizeof(double), cudaMemcpyHostToDevice);


    double *rho;
    double rhotemp = 1.0;
    CHECK(cudaMalloc(&rho, sizeof(double)));
    cudaMemcpy(rho, &rhotemp,  sizeof(double ), cudaMemcpyHostToDevice);

    cublasStatus_t stat;

    double **Ap, **Aip;
    double **A1 ,**Ai;
    CHECK(cudaMalloc(&Ap, sizeof(double *)));
    A1 = (double **)malloc(sizeof(double *));
    CHECK(cudaMalloc(&Aip, sizeof(double *)));
    Ai = (double **)malloc(sizeof(double *));
    A1[0] = temp2_cu;
    Ai[0] = Ainv;
    cudaMemcpy(Ap, A1,  sizeof(double *), cudaMemcpyHostToDevice);
    cudaMemcpy(Aip, Ai,  sizeof(double *), cudaMemcpyHostToDevice);


    Account *account_cu;CHECK(cudaMalloc(&account_cu, sizeof(Account)));
    //CHECK(cudaMemcpy(&account_cu->C, &this->model->C, sizeof(double),cudaMemcpyHostToDevice));
    //CHECK(cudaMemcpy(&account_cu->w_last[0], this->model->w_last, n*sizeof(double),cudaMemcpyHostToDevice));

    dim3 block(n,n/128+ 1);

    printf("start solving");
    int iter = 0;
    for (iter = 1; iter < opts.max_iter;++iter){


        update_f<<<block,128>>>(f_cu,R_cu,z_cu,u_cu,rho,n);


        // Pb_b = - P_cu * f
        cublasDgemv(handle, CUBLAS_OP_N, n, n, &alpha_1,
                    P_cu, n, f_cu, 1, &beta, Pb_cu,1);

        // create invs n * n
        generate_invs<<<block,128>>>(invs_cu,S_cu, rho, n);

        // temp1 = APt*invs  m*n
        matrix_times(handle,APt_cu,invs_cu,m,n,n,temp1_cu);

        // temp2 = temp1*PAt
        matrix_times(handle,temp1_cu,PAt_cu,m,n,m,temp2_cu);

        // temp2 = inv(temp)
        stat = cublasDmatinvBatched(handle, m,  Ap, m,  Aip, m, info, 1);

        // temp_b = APt * invs * Pb_b = temp1_cu * Pb_b n*1
        cudaMemcpy(x2_cu, b_cu,  n * sizeof(double), cudaMemcpyDeviceToDevice);
        cublasDgemv(handle, CUBLAS_OP_N, m, n, &alpha1,
                    temp1_cu, m, Pb_cu, 1, &alpha_1, x2_cu,1);

        // x2 = Ainv * temp_b  m * 1
        cublasDgemv(handle, CUBLAS_OP_N, m, m, &alpha1,
                    Ainv, m, x2_cu, 1, &beta, x2_cu,1);

        // Pb_cu = Pb - PAt*x2
        cublasDgemv(handle, CUBLAS_OP_N, n, m, &alpha_1,
                    PAt_cu, n, x2_cu, 1, &alpha1, Pb_cu,1);
        cublasDgemv(handle, CUBLAS_OP_N, n, n, &alpha_1,
                    invs_cu, n, Pb_cu, 1, &beta, Pb_cu,1);
        cublasDgemv(handle, CUBLAS_OP_N, n, n, &alpha_1,
                    Pt_cu, n, Pb_cu, 1, &beta, x_cu,1);


        // zold = z
        cudaMemcpy(z_old, z_cu,  n_size, cudaMemcpyDeviceToDevice);

        // x_hat = gamma_*x + (1 - gamma_)*zold;
        update_x_hat<<<block,128>>>(x_hat_cu, x_cu,z_old,gamma_,n);

        // tt_cu =   x_hat_cu + u/rho;
        update_tt_cu<<<block,128>>>(tt_cu,x_hat_cu,u_cu,rho,n);

        // update z = min(max(x_hat + u/rho ,lb),ub) ;
        // for log problem there is no need
        
        //update_z_cu<<<block,128>>>(z_cu,tt_cu,lb_cu,ub_cu,n);


        //dfp_choose<<<block,128>>>(z_cu,ss_s,ss_e,sb_s,sb_e,sm_s,sm_e,tt_cu,rho,
        //                          lb_cu,ub_cu,account_cu,n);

        cublasDnrm2(handle,n,z_cu,1,&r_norm); //r_norm
//        if( isnan(r_norm))  {
//            printf("z_cu  ");
//            show_res(z_cu,n);
//            printf("nan ! error");
//            exit(-1);
//        }

        my_newton<<<block,128>>>(z_cu,
                            tt_cu,
                            rho,n);
        // unew = u + rho* (x_hat - z);
        update_u_new<<<block,128>>>(u_new, u_cu,rho,x_hat_cu,z_cu,n);

        // objective function do not need

        // norms
        // x_cu - z_cu; z_cu - zold_cu; norm(x) norm(z) norm(u_new)

        if(iter == 1){
             update_l0_hat0_and_so_on<<<block,128>>>(l0_cu,lhat0_cu,Bv0_cu,Au0_cu,
                                                     u_new,u_cu,x_hat_cu,z_old,
                                                     z_cu,rho,n);
        }else{
            if( ((iter%freq) == 0) && (iter > siter)  && (iter < eiter) ){
                update_hat_and_so_on<<<block,128>>>(lhat_cu,Bv_cu,Au_cu,u_cu,x_hat_cu,z_old,z_cu,rho,n);
                aradmm_estimate(handle,rho,gamma_,Au_cu,Au0_cu,
                                lhat_cu,lhat0_cu,Bv_cu,Bv0_cu,u_new,l0_cu,n);

                update_l0_hat0_and_so_on<<<block,128>>>(l0_cu,lhat0_cu,Bv0_cu,Au0_cu,
                                                        u_new,u_cu,x_hat_cu,z_old,
                                                        z_cu,rho,n);
            }
        }
        // u = unew
        cudaMemcpy(u_cu, u_new,  n_size ,cudaMemcpyDeviceToDevice);

        // calculate norms


        cublasDnrm2(handle,n,x_cu,1,&eps_dual); //xnorm
        cublasDnrm2(handle,n,z_cu,1,&eps_pri);  // znorm
        eps_pri = sqrtn + opts.RELTOL * max(eps_dual,eps_pri);
        cublasDnrm2(handle,n,u_new,1,&eps_dual);
        eps_dual  = sqrtn + opts.RELTOL *eps_dual;

        // norm(x-z)
        prepare_for_norm<<<block,128>>>(norm_cu,norm_cu2,x_cu,z_cu,rho,z_old,n);

        cublasDnrm2(handle,n,norm_cu2,1,&r_norm); //r_norm
        cublasDnrm2(handle,n,norm_cu,1,&s_norm); //s_ norm

        if( isnan(r_norm))  {
            printf("x_cu  ");
            show_res(x_cu,n);

            printf("z_cu  ");
            show_res(z_cu,n);

            printf("norm_cu2  ");
            show_res(norm_cu2,n);

            printf("nan ! error");
            exit(-1);
        }
        printf("%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t  \n", iter,
                     r_norm ,  eps_pri ,
                     s_norm ,  eps_dual );



#ifdef _RECORD_

        this->record.errors[iter] = r_norm;
        double value = 0.0;

//        calculate the objective value
        //cublasDgemv(handle, CUBLAS_OP_N, n, n, &alpha1,
//                    sigma_cu, n, x_cu, 1, &beta, norm_cu,1);
        //cublasDdot (handle,  n,x_cu, 1,norm_cu, 1,&this->record.values[iter]);
        //cublasDdot (handle,  n,x_cu, 1,R_cu, 1,&value);

        //this->record.values[iter] -= value;
        this->record.values[iter] = objective_log(block,handle, sigma_cu,norm_cu,z_cu,n);

//        calcost<<<block,128>>>(norm_cu, x_cu,ss_s,ss_e,sb_s,sb_e,sm_s,sm_e,
//                                    account_cu,n);
        cublasDasum(handle, n, norm_cu ,1, &value);
        //this->record.values[iter] += value;

         printf("value = %.4f\n",this->record.values[iter]);

#endif
        if ( (r_norm  <  eps_pri)  && (s_norm  <  eps_dual ))
            break;





    }

    //calculate objective
    CHECK(cudaMemcpy(&this->result[0], x_cu, n_size,cudaMemcpyDeviceToHost));





    return iter;
}
double objective_log(dim3 &block,cublasHandle_t handle, double * sigma_cu, double *norm_cu,double *x_cu,int n){
    double result = 0.0;
    double value = 0.0;
    cublasDgemv(handle, CUBLAS_OP_N, n, n, &alpha1,
                sigma_cu, n, x_cu, 1, &beta, norm_cu,1);
    cublasDdot (handle,  n,x_cu, 1,norm_cu, 1,&result);
    callogcost<<<block,128>>>(norm_cu, x_cu,n);
    cublasDasum(handle, n, norm_cu ,1, &value);
    result += value;
    return result;
}
