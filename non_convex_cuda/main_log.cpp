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
    int rng;
    if (argc <= 2) {
        size = 10;
        rng = 1;
    }    
    else
    {
        size = atoi(argv[1]);
        rng = atoi(argv[2]);
    }
    std::string file_path = "../data/logportfoliodata/" + std::to_string(rng) +"/" + std::to_string(size) +"/";

    std::string Apath = file_path + "Aeq.mat";

    std::string Sigmapath = file_path + "Sigma.mat";
    std::string bpath = file_path + "beq.mat";  
    std::string lbpath = file_path + "lb.mat";
 
 

    printf("Reading Data\n");

    std::vector<double> A  = readFile<double>(Apath.c_str());
    
    std::vector<double> Sigma = readFile<double>(Sigmapath.c_str());
    std::vector<double> b = readFile<double>(bpath.c_str());
    std::vector<double> lb = readFile<double>(lbpath.c_str());
    //std::vector<double> ub = readFile<double>(ubpath.c_str());


 
//    for (int i = 0 ; i < b.size();++i)
//    std::cout << b[i] << std::endl;

    int m = b.size();
    int n = lb.size();
    std::vector<double> Er(n, 0);
    std::vector<double> ub(n, 1);
    printf("Create model\n");

    PortfolioModel model = PortfolioModel(m,n);
    int flag;

    printf("Register data\n");

    flag = model.RegisterData(Sigma,A,b,Er,lb,ub);
    if (flag)
        return -1;

 

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
    std::string output_path = "./output_log/";




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
