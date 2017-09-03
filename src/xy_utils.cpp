#include "xy_utils.h"
#include <iostream> 
#include <cmath>

double max(const std::vector<double> &delta){
  double temp=delta.at(0);
  for(unsigned int i = 1; i < delta.size(); ++i){
    if(delta.at(i)>temp)temp=delta.at(i);
  }
  return temp;
}

double absmax(const std::vector<double> &delta){
  double temp=std::fabs(delta.at(0));
  double tempa;
  for(unsigned int i = 1; i < delta.size(); ++i){
    tempa=std::fabs(delta.at(i));
    if(tempa>temp)temp=tempa;
  }
  return temp;
}

double abs_av_delta(const std::vector<double> &delta){
  double temp=0.;
  double tempa;
  for(unsigned int i = 0; i < delta.size(); ++i){
    tempa=std::fabs(delta.at(i));
    temp=temp+tempa*tempa;
  }
  temp=temp/delta.size();
  temp=std::sqrt(temp);

  return temp;
}

void show_matrix(const std::vector<std::vector<double>> &data)
{
    unsigned int nx=data.size();
    unsigned int ny=0;
    if(nx>0)ny=data.at(0).size();
    std::cout << "nx= " << nx << " ; ny= " << ny << std::endl;

    for(unsigned int j=0;j<nx;j++){
      for(unsigned int i=0;i<ny;i++){
        std::cout << data.at(j).at(i) << " ";
      }
      std::cout << std::endl;
    }
}

std::vector<std::vector<double>> matrix_product(const std::vector<std::vector<double>> &A,const std::vector<std::vector<double>> &B)
{
    std::vector<std::vector<double>> P;
    P.clear();

    unsigned int nx=A.size();
    unsigned int ny=0;
    if(nx>0)ny=A.at(0).size();

    unsigned int nx2=B.size();
    unsigned int ny2=0;
    if(nx2>0)ny2=B.at(0).size();

    if(ny!=nx2){
      std::cout << "Matrix product input dimensions do not match." << std::endl;
      return P;
    }

    std::vector<double> tempv;
    double temp;
    for(unsigned int i=0;i<nx;i++){
      tempv.clear();
      for(unsigned int k=0;k<ny2;k++){
        temp=0.;
        for(unsigned int j=0;j<ny;j++){
          temp=temp+A.at(i).at(j)*B.at(j).at(k);
        }
        tempv.push_back(temp);
      }
      P.push_back(tempv);
    }


    return P;
}
