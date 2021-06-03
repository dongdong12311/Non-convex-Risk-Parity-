#include "admmoption.h"

AdmmOption::AdmmOption()
{
    adp_freq = 2;
    adp_start_iter = 10;
    adp_end_iter = 1000;
    max_iter = 2000;
    orthval = 0.2;
    beta_scale = 2;
    res_scale = 0.1;
    gamma = 1;
    ABSTOL = 1e-5;
    RELTOL = 1e-5;
}
