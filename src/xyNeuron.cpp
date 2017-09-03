
#include "xyNeuron.h"
#include <iostream>
#include <cmath>



// constructor
xyNeuron::xyNeuron(unsigned int size)
: bias(0.)
{
    poshidprobs.resize(size);
    neghidprobs.resize(size);
}


void xyNeuron::Set_bias(const double &b_in)
{
  bias=b_in;
}

void xyNeuron::Set_biasinc(const double &b_in)
{
  biasinc=b_in;
}

void xyNeuron::Set_size(unsigned int size)
{
    poshidprobs.resize(size);
    neghidprobs.resize(size);
}

void xyNeuron::Set_p_prob(unsigned int id, const double &val)
{
    if(id>=poshidprobs.size()||id<0){
      std::cout << "Set_p_prob out of idex range" << std::endl;
      return;
    }
    poshidprobs.at(id)=val;
}

void xyNeuron::Set_n_prob(unsigned int id, const double &val)
{
    if(id>=neghidprobs.size()||id<0){
      std::cout << "Set_n_prob out of idex range" << std::endl;
      return;
    }
    neghidprobs.at(id)=val;
}

// sigmoid function for output
inline double xyNeuron::sigmoid(const double &a, const double &p)
const
{
	return 1./(1. + std::exp(-a/p));
}

inline double xyNeuron::dsigmoid(const double &sig, const double &p)
const
{
    // sigmoid f(x, a) derivative f'(x, a) = 1/a * f(x, a) * (1 - f(x, a))
    return 1./p * sig*(1. - sig);
}
