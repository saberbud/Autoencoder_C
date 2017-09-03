#ifndef xyRBMUNIT_H
#define xyRBMUNIT_H

#include <vector>
#include "xyNeuron.h"

class xyRBMUNIT
{
public:
    xyRBMUNIT(unsigned int v_size=1, unsigned int h_size=1);

    void Set_v_size(unsigned int size);
    void Set_h_size(unsigned int size);
    unsigned int GetvSize() const {return xyvneurons.size();};
    unsigned int GethSize() const {return xyhneurons.size();};

    void init_h_w();
    void init_vh_w();
    void init_v_bias();
    void init_h_bias();
    void init_prods();

    void init_vh_inc();

    void show_v_bias();
    void show_h_bias();
    void show_vishid_w();
    void show_poshidstates();

    void showdata(const std::vector<std::vector<double>> &data);
    void test_in_call();
    void Set_train_para(const std::vector<double> &para);
    int rbm_train_one_batch(const std::vector<std::vector<double>> &data);

    double Get_err() const {return err;};
    const std::vector<std::vector<double>> Get_poshidprobs() const {return poshidprobs;};
    const std::vector<std::vector<double>> Get_vishid() const {return vishid_w;};
    std::vector<double> Get_v_bias();
    std::vector<double> Get_h_bias();

private:
    std::vector<xyNeuron> xyvneurons;
    std::vector<xyNeuron> xyhneurons;
    std::vector<double> hid_w;
    std::vector<std::vector<double>> vishid_w;
    std::vector<std::vector<double>> posprods;
    std::vector<std::vector<double>> negprods;
    std::vector<std::vector<double>> poshidprobs;
    std::vector<std::vector<double>> poshidstates;
    std::vector<std::vector<double>> vishid_inc;

    double err;

    double epsilonw;
    double epsilonvb;
    double epsilonhb;
    double weightcost;
    double momentum;

};
#endif
