#include<iostream>
#include<vector>
#include <stdio.h>
#include<Common.h>
#include<string>
#include<portfoliomodel.h>
#include<portfoliosolver.h>
template <typename T>
void show_vector( T &data){
    for (int i = 0;i<data.size();++i){
        std::cout << data[i] << " ";
    }std::cout << std::endl;
}
#include <ctime>
#include <chrono>
#include<fstream>
#include<string>
void output_result(double *data,int n, const char *filename){
    FILE *file = NULL;

    // output the w
    file = fopen(filename, "w+");
    if(file){
        for (int i = 0 ; i < n; ++i)
         {

            fprintf(file, "%.16f, ", data[i]) ;
        }
        fclose(file);
    }
    else{
       printf("error for output\n") ;
    }


}

int main(int argc, char **argv){
    int size;
    if (argc <= 1)
        size = 10;
    else
        size = atoi(argv[1]);
    std::string file_path = "../data/portfoliodata/" + std::to_string(size) +"/";

    std::string Apath = file_path + "A_admm.mat";
    std::string Erpath = file_path + "Er_admm.mat";
    std::string Sigmapath = file_path + "Sigma_admm.mat";
    std::string bpath = file_path + "b_admm.mat";
    std::string lbpath = file_path + "lb_admm.mat";
    std::string ubpath = file_path + "ub_admm.mat";
    std::string n_Sbpath = file_path + "n_Sb.mat";

    std::string n_Smpath = file_path + "n_Sm.mat";
    std::string n_Sspath = file_path + "n_Ss.mat";
    std::string Cpath = file_path + "C.mat";
    std::string w_lastpath = file_path + "w_last.mat";

    printf("Reading Data\n");

    std::vector<double> A  = readFile<double>(Apath.c_str());
    std::vector<double> Er = readFile<double>(Erpath.c_str());
    std::vector<double> Sigma = readFile<double>(Sigmapath.c_str());
    std::vector<double> b = readFile<double>(bpath.c_str());
    std::vector<double> lb = readFile<double>(lbpath.c_str());
    std::vector<double> ub = readFile<double>(ubpath.c_str());
    std::vector<int> n_Sb = readFile(n_Sbpath.c_str());
    for (int i = 0; i < n_Sb.size();++i) n_Sb[i] = n_Sb[i] - 1;
    std::vector<int> n_Sm = readFile(n_Smpath.c_str());
    for (int i = 0; i < n_Sm.size();++i) n_Sm[i] = n_Sm[i] - 1;
    std::vector<int> n_Ss = readFile(n_Sspath.c_str());
    for (int i = 0; i < n_Ss.size();++i) n_Ss[i] = n_Ss[i] - 1;

    std::vector<double> C = readFile<double>(Cpath.c_str());
    std::vector<double> w_last = readFile<double>(w_lastpath.c_str());


    int m = b.size();
    int n = lb.size();

    printf("Create model\n");

    PortfolioModel model = PortfolioModel(m,n);
    int flag;

    printf("Register data\n");

    flag = model.RegisterData(Sigma,A,b,Er,lb,ub);
    if (flag)
        return -1;

    printf("RegisterIndex\n");

    model.RegisterIndex(n_Sb,n_Sm,n_Ss);
    model.RegisterParameters(C[0],w_last.data());

    printf("Create Solver\n");

    portfoliosolver Solver = portfoliosolver();

    Solver.RegisterModel(&model);

    printf("Initialization\n");

    Solver.init();

    printf("Start solving\n");

    Option opts = Option();
    auto start = std::chrono::steady_clock::now();
    int iter = Solver.SolveModel(opts);
    auto end = std::chrono::steady_clock::now();




    auto diff = end - start;
    double t = (double)std::chrono::duration <double, std::milli> (diff).count()/1000.0;
    std::cout << t << " s" << std::endl;

    // output the objective value and the error
    printf("Outputs the result\n");
    std::string output_path = "./output/";




#ifdef _RECORD_

    std::string output_solution = output_path + "solution_of_" + std::to_string(size) + ".csv";
    output_result(&Solver.result[0],n,output_solution.c_str());
    std::string output_value = output_path + "true_value_on_" + std::to_string(size) + ".csv";
    output_result(&Solver.record.values[0],iter,output_value.c_str());
    std::string output_error = output_path + "error_on_" + std::to_string(size) + ".csv";
    output_result(&Solver.record.errors[0],iter,output_error.c_str());

#else
    // if we do not record the value we record the time
    std::string output_time = output_path + "time_spent_on_" + std::to_string(size) + ".csv";
    output_result(&t,1,output_time.c_str());
#endif




    return 0;
}
