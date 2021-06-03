#ifndef PARAMETERS_H
#define PARAMETERS_H
class Option
{
public:
        int adp_freq;
        int adp_start_iter;
        int adp_end_iter;
        double orthval;
        double beta_scale;
        double res_scale;
        double gamma;
        int max_iter;
        double ABSTOL;
        double RELTOL;
        Option(){
            adp_freq = 2;
            adp_start_iter = 10;
            adp_end_iter = 1000;
            max_iter = 2000;
            orthval = 0.2;
            beta_scale = 2.0;
            res_scale = 0.1;
            gamma = 1.0;
            ABSTOL = 1e-10;
            RELTOL = 1e-10;
        }
};
#endif // PARAMETERS_H
