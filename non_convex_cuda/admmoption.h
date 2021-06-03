#ifndef ADMMOPTION_H
#define ADMMOPTION_H


class AdmmOption
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
    AdmmOption();
};

#endif // ADMMOPTION_H
