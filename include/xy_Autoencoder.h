#ifndef xy_Autoencoder_H
#define xy_Autoencoder_H

#include <vector>
#include <string>
#include "xy_W_bag.h"
#include "xydata_mini_batch.h"
#include "xyrbmunit.h"
#include "xy_MLP.h"

class xy_Autoencoder
{
public:
    xy_Autoencoder();

    void train_data_in(const std::vector<std::vector<double>> &data, const std::vector<std::vector<double>> &targ, unsigned int &nb);

    void test_data_in(const std::vector<std::vector<double>> &data, const std::vector<std::vector<double>> &targ, unsigned int &nb);

    void Set_net_strucuture(const std::vector<unsigned int> &st);

    void init_RBMs();

    void Set_RBM_para(std::vector<double> &para, unsigned int &id);

    int train_RBMs(unsigned int &nep, std::vector<unsigned int> &sw, std::vector<std::vector<double>> &swp);

    void Set_mlp_para(std::vector<double> &para);

    void train_mlp(unsigned int &nep, unsigned int &n_o, double &leng);

    void set_mlp_w_W();

    std::vector<unsigned int> n_test_miss();

    xy_W_bag calc_target(xy_W_bag &A);

    unsigned int n_miss_match(const std::vector<std::vector<double>> &A, const std::vector<std::vector<double>> &B);

    void show_ws(std::vector<xy_W_bag> &A, unsigned int &id);

    void Set_one_W(xy_W_bag &A, unsigned int &id);

    unsigned int max_idx(const std::vector<double> &A);

    std::vector<double> rand_vec(unsigned int &ny);
    std::vector<std::vector<double>> rand_mat(unsigned int &nx, unsigned int &ny);

    void read_W(std::string &pref);

    void print_W(double &type, std::string &pref);

    void print_to_file(xy_W_bag &A, std::string &name);

private:
    xydata_mini_batch train_data;
    xydata_mini_batch test_data;

    std::vector<xy_W_bag> Weights;
    std::vector<xy_W_bag> train_wps;
    std::vector<xy_W_bag> test_wps;

    std::vector<xyRBMUNIT> rbmv;
    xy_MLP mlp;

};
#endif
