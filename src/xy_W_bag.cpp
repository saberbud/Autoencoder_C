#include "xy_W_bag.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>


// constructor
xy_W_bag::xy_W_bag()
{
    // place holder
}

void xy_W_bag::W_in(const std::vector<std::vector<double>> &data)
{
    bag_W.clear();
    bag_W=data;
}

void xy_W_bag::b_in(const std::vector<double> &data)
{
    bag_b.clear();
    bag_b=data;
}

void xy_W_bag::d_W_in(const std::vector<std::vector<double>> &data)
{
    bag_d_W.clear();
    bag_d_W=data;
}

void xy_W_bag::d_b_in(const std::vector<double> &data)
{
    bag_d_b.clear();
    bag_d_b=data;
}

int xy_W_bag::Wb_comb()
{
    bag_Wb.clear();
    unsigned int nx=bag_W.size();
    unsigned int ny=bag_b.size();

    for(unsigned int i=0;i<nx;i++){
      if(bag_W.at(i).size()!=ny){
        std::cout << "W and b dimension mismatch." << std::endl;
        return 1;
      }
    }

    bag_Wb=bag_W;
    bag_Wb.push_back(bag_b);

    return 0;
}
