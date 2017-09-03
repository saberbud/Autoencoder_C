#include "xydata_mini_batch.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <cmath>


// constructor
xydata_mini_batch::xydata_mini_batch()
{
    data_batch.clear();
    mini_batch.clear();
    targ_batch.clear();
    mini_batch_t.clear();
}

void xydata_mini_batch::data_in(const std::vector<std::vector<double>> &data)
{
    data_batch.clear();
    data_batch=data;
}

void xydata_mini_batch::targ_in(const std::vector<std::vector<double>> &data)
{   
    targ_batch.clear();
    targ_batch=data;
}

void xydata_mini_batch::data_shuffle()
{
    unsigned int size=data_batch.size();
    unsigned int sizep=size-1;

    if(size<1){
      std::cout << "Data not been input." << std::endl;
      return;
    }

    if(size!=targ_batch.size()){
      std::cout << "Target not been input." << std::endl;
      return;
    }

    std::mt19937 rng;
    rng.seed(std::random_device()());
    std::uniform_int_distribution<unsigned int> uni_dist(0., sizep);

    unsigned int tempi;
    unsigned int ny;
    std::vector<double> tempv;

    if(size>0)ny=data_batch.at(0).size();

    for(unsigned int i=0;i<size;i++){
      if(ny!=data_batch.at(i).size()){
        std::cout << "data ny not constant. Not do shuffle" << std::endl;
        return;
      }

      tempi=uni_dist(rng);
      // std::cout << tempi << std::endl;

      tempv=data_batch.at(i);
      data_batch.at(i)=data_batch.at(tempi);
      data_batch.at(tempi)=tempv;

      tempv=targ_batch.at(i);
      targ_batch.at(i)=targ_batch.at(tempi);
      targ_batch.at(tempi)=tempv;
    }

}

xy_W_bag xydata_mini_batch::Get_mini_batch(unsigned int &id)
{
    xy_W_bag P;

    unsigned int size=mini_batch.size();
    unsigned int sizep=size-1;

    // std::cout << "size= " << size << " ; sizep= " << sizep << " , id= " << id << std::endl;

    if(id>sizep || id>size){
      std::cout << "Get_mini_batch: index out of range." << std::endl;
      return P;
    }

    P=mini_batch.at(id);

    return P;
}

xy_W_bag xydata_mini_batch::Get_mini_batch_t(unsigned int &id)
{   
    xy_W_bag P;
    
    unsigned int size=mini_batch_t.size();
    unsigned int sizep=size-1;
    
    // std::cout << "size= " << size << " ; sizep= " << sizep << " , id= " << id << std::endl;
    
    if(id>sizep || id>size){
      std::cout << "Get_mini_batch: index out of range." << std::endl;
      return P;
    }
    
    P=mini_batch_t.at(id);
    
    return P;
}

void xydata_mini_batch::make_mini_batch(unsigned int &nb)
{
    unsigned int size=data_batch.size();

    if(size<1){
      std::cout << "No data input, cannot make mini batch." << std::endl;
      return;
    }

    if(size!=targ_batch.size()){
      std::cout << "Target not been input, cannot make mini batch" << std::endl;
      return;
    }

    unsigned int nres;
    nres=size%nb;
    // std::cout << "nres= " << nres << std::endl;

    if(nres!=0){
      std::cout << "N batches must be a fraction int of data dimension" << std::endl;
      return;
    }

    unsigned int ncase;
    ncase=size/nb;
    // std::cout << "ncase= " << ncase << std::endl;

    mini_batch.resize(nb);
    mini_batch_t.resize(nb);
    // std::cout << "# of mini batches= " << mini_batch.size() << std::endl;

    std::vector<std::vector<double>> tempm;
    std::vector<std::vector<double>> tempm_t;
    unsigned int tempi;
    for(unsigned int i=0;i<nb;i++){
      tempm.clear();
      tempm_t.clear();
      for(unsigned int j=0;j<ncase;j++){
        tempi=i*ncase+j;
        // std::cout << tempi << std::endl;
        tempm.push_back(data_batch.at(tempi));
        tempm_t.push_back(targ_batch.at(tempi));
      }
      mini_batch.at(i).W_in(tempm);
      mini_batch_t.at(i).W_in(tempm_t);
    }


}
