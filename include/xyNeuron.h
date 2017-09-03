#ifndef xy_NEURON_H
#define xy_NEURON_H

#include <vector>

class xyNeuron
{
public:
    xyNeuron(unsigned int size = 1);

    void SetWeights(const std::vector<double> &w);
    std::vector<double> GetWeights() const;

    void Set_bias(const double &b_in);
    double Get_bias() const {return bias;};

    void Set_biasinc(const double &b_in);
    double Get_biasinc() const {return biasinc;};

    void Set_size(unsigned int size);
    unsigned int GetSize() const {return poshidprobs.size();};

    void Set_p_prob(unsigned int id, const double &val);
    double Get_p_prob(unsigned int id) const {return poshidprobs.at(id);};

    void Set_n_prob(unsigned int id, const double &val);
    double Get_n_prob(unsigned int id) const {return neghidprobs.at(id);};

private:
    double sigmoid(const double &a, const double &p) const;
    double dsigmoid(const double &a, const double &p) const;

private:
    std::vector<double> poshidprobs;
    std::vector<double> neghidprobs;
    double bias;
    double biasinc;
};

#endif
