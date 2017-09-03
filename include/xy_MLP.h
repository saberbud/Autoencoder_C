#ifndef xy_MLP_H
#define xy_MLP_H

#include <vector>
#include "xy_W_bag.h"

class xy_MLP
{
public:
    xy_MLP();

    xy_W_bag calc_prob(xy_W_bag &A, xy_W_bag &B);
    xy_W_bag calc_targ(xy_W_bag &A, xy_W_bag &B);
    double calc_err(xy_W_bag &A, xy_W_bag &B);

    std::vector<xy_W_bag> calc_d_class(std::vector<xy_W_bag> &A);
    std::vector<xy_W_bag> calc_d_base(std::vector<xy_W_bag> &A);

    void init_Ps(const std::vector<double> &A);

    void input_weights(std::vector<xy_W_bag> &A);
    std::vector<xy_W_bag> &Get_weights() {return ws;};

    void input_data(xy_W_bag &A);
    std::vector<xy_W_bag> &Get_probs() {return wps;};

    void input_targ(xy_W_bag &A);
    xy_W_bag &Get_targ() {return targorg;};

    std::vector<xy_W_bag> &Get_cbag() {return cbag;};

    void eval_probs_base();
    void eval_probs_targ();

    void eval_d_class();
    void eval_d_base();

    void eval_err();

    int check_run();

    std::vector<xy_W_bag> add_d_W(std::vector<xy_W_bag> &A, double &x);
    std::vector<xy_W_bag> dW_to_W(std::vector<xy_W_bag> &A);
    std::vector<xy_W_bag> W_to_dW(std::vector<xy_W_bag> &A, std::vector<xy_W_bag> &B);
    std::vector<xy_W_bag> mat_sum(double &fa, std::vector<xy_W_bag> &A, double &fb, std::vector<xy_W_bag> &B);
    double s_prod(std::vector<xy_W_bag> &A, std::vector<xy_W_bag> &B);

    void xy_minimize(double &f_type, double &f_leng);

private:
    std::vector<xy_W_bag> ws;  //w1,2,3,...class
    std::vector<xy_W_bag> wps;  //w_probs: data,1,2,3...
    xy_W_bag targorg;  //original target

    std::vector<xy_W_bag> cbag;  //for checking

    xy_W_bag IX_xy;

    double error;

    int P_init;
    double P_INT;
    double P_EXT;
    double P_MAX;
    double P_RATIO;
    double P_SIG;
    double P_RHO;

};
#endif
