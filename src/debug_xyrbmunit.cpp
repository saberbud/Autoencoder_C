#include "xyrbmunit.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <cmath>

// constructor
xyRBMUNIT::xyRBMUNIT(unsigned int v_size, unsigned int h_size)
: epsilonw(0.1),epsilonvb(0.1),epsilonhb(0.1),weightcost(0.0002),momentum(0.5)
{
    xyvneurons.resize(v_size);
    xyhneurons.resize(h_size);
}


void xyRBMUNIT::Set_v_size(unsigned int size)
{
    xyvneurons.resize(size);
}

void xyRBMUNIT::Set_h_size(unsigned int size)
{   
    xyhneurons.resize(size);
}

void xyRBMUNIT::init_h_w()
{
    std::mt19937 rng;
    rng.seed(std::random_device()());
    // std::uniform_real_distribution<double> uni_dist(-1., 1.);
    std::normal_distribution<double> norm_dist(0., 1.);

    unsigned int hsize=xyhneurons.size();
    double temp;

    hid_w.clear();

    for(unsigned int i=0;i<hsize;i++){
      // hid_w.push_back(uni_dist(rng));
      temp=norm_dist(rng);
      temp=temp*0.1;
      hid_w.push_back(temp);
    }


    for(unsigned int i=0;i<hsize;i++){
      std::cout << hid_w.at(i) << std::endl;
    }
    
}

void xyRBMUNIT::init_vh_w(){
    std::mt19937 rng;
    rng.seed(std::random_device()());
    std::normal_distribution<double> norm_dist(0., 1.);

    unsigned int vsize=xyvneurons.size();
    unsigned int hsize=xyhneurons.size();
    double temp;

    vishid_w.clear();
    for(unsigned int j=0;j<vsize;j++){
      hid_w.clear();
      for(unsigned int i=0;i<hsize;i++){
        temp=norm_dist(rng);
        temp=temp*0.1;
        hid_w.push_back(temp);
      }
      vishid_w.push_back(hid_w);
    }

    vishid_inc.clear();
    for(unsigned int j=0;j<vsize;j++){
      hid_w.clear();
      for(unsigned int i=0;i<hsize;i++){
        temp=0.;
        hid_w.push_back(temp);
      }
      vishid_inc.push_back(hid_w);
    }

    std::cout << "vishid_w" << std::endl;
    for(unsigned int j=0;j<vsize;j++){
      for(unsigned int i=0;i<hsize;i++){
        std::cout << vishid_w.at(j).at(i) << " ";
      }
      std::cout << std::endl;
    }

    std::cout << "vishid_inc" << std::endl;
    for(unsigned int j=0;j<vsize;j++){
      for(unsigned int i=0;i<hsize;i++){
        std::cout << vishid_inc.at(j).at(i) << " ";
      }
      std::cout << std::endl;
    }

}

void xyRBMUNIT::init_v_bias()
{
    double a=0.;
    unsigned int vsize=xyvneurons.size();
    for(unsigned int j=0;j<vsize;j++){
      xyvneurons.at(j).Set_bias(a);
      xyvneurons.at(j).Set_biasinc(a);
    }

    for(unsigned int j=0;j<vsize;j++){
      std::cout << j << " v bias " << xyvneurons.at(j).Get_bias() << " | " << xyvneurons.at(j).Get_biasinc() << std::endl;
    }

}

void xyRBMUNIT::init_h_bias()
{
    double a=0.;
    unsigned int hsize=xyhneurons.size();
    for(unsigned int j=0;j<hsize;j++){
      xyhneurons.at(j).Set_bias(a);
      xyhneurons.at(j).Set_biasinc(a);
    }

    for(unsigned int j=0;j<hsize;j++){
      std::cout << j << " h bias " << xyhneurons.at(j).Get_bias() << " | " << xyhneurons.at(j).Get_biasinc() << std::endl;
    }

}

void xyRBMUNIT::init_prods()
{
    unsigned int vsize=xyvneurons.size();
    unsigned int hsize=xyhneurons.size();
    double temp=0.;
    std::vector<double> data;

    posprods.clear();
    negprods.clear();
    for(unsigned int j=0;j<vsize;j++){
      data.clear();
      for(unsigned int i=0;i<hsize;i++){
        data.push_back(temp);
      }
      posprods.push_back(data);
      negprods.push_back(data);
    }

    for(unsigned int j=0;j<vsize;j++){
      for(unsigned int i=0;i<hsize;i++){
        std::cout << " || " << posprods.at(j).at(i) << " | " << negprods.at(j).at(i);
      }
      std::cout << std::endl;
    }
}

void xyRBMUNIT::init_vh_inc()
{
    unsigned int vsize=xyvneurons.size();
    unsigned int hsize=xyhneurons.size();
    double temp=0.;
    std::vector<double> data;

    vishid_inc.clear();
    for(unsigned int j=0;j<vsize;j++){
      data.clear();
      for(unsigned int i=0;i<hsize;i++){
        data.push_back(temp);
      }
      vishid_inc.push_back(data);
    }

    for(unsigned int j=0;j<vsize;j++){
      for(unsigned int i=0;i<hsize;i++){
        std::cout << vishid_w.at(j).at(i) << " | " << vishid_inc.at(j).at(i) << " || ";
      }
      std::cout << std::endl;
    }
}

void xyRBMUNIT::show_v_bias()
{
    unsigned int vsize=xyvneurons.size();
    std::cout<< "Show visbias, # of neurons: " << vsize << std::endl;

    for(unsigned int k=0;k<vsize;k++){
      std::cout << xyvneurons.at(k).Get_bias() << " ";
    }
    std::cout << std::endl;
}

void xyRBMUNIT::show_h_bias()
{
    unsigned int hsize=xyhneurons.size();
    std::cout<< "Show hidbias, # of neurons: " << hsize << std::endl;

    for(unsigned int j=0;j<hsize;j++){
      std::cout << xyhneurons.at(j).Get_bias() << " ";
    }
    std::cout << std::endl;
}

void xyRBMUNIT::show_vishid_w()
{
    unsigned int vsize=xyvneurons.size();
    unsigned int hsize=xyhneurons.size();
    std::cout<< "Show vishid_w, size: " << vsize << " X " << hsize << std::endl;

    for(unsigned int k=0;k<vsize;k++){
      for(unsigned int j=0;j<hsize;j++){
        std::cout << vishid_w.at(k).at(j) << " ";
      }
      std::cout << std::endl;
    }
}

void xyRBMUNIT::show_poshidstates()
{
    unsigned int numcases=poshidstates.size();
    unsigned int hsize=xyhneurons.size();
    std::cout<< "Show poshidstates: " << numcases << " X " << hsize << std::endl;

    for(unsigned int i=0;i<numcases;i++){
      for(unsigned int j=0;j<hsize;j++){
        std::cout << poshidstates.at(i).at(j) << " ";
      }
      std::cout << std::endl;
    }
}

void xyRBMUNIT::showdata(const std::vector<std::vector<double>> &data)
{
    unsigned int numcases=data.size();
    unsigned int numdims;
    if(numcases>0)numdims=data.at(0).size();
    std::cout << "numcases= " << numcases << " ; numdims= " << numdims << std::endl;

    for(unsigned int j=0;j<numcases;j++){
      for(unsigned int i=0;i<numdims;i++){
        std::cout << data.at(j).at(i) << " ";
      }
      std::cout << std::endl;
    }
    test_in_call();
}

void xyRBMUNIT::test_in_call()
{
    double a=1.;
    double b=std::exp(-a);
    std:: cout << "test_in_call " << b << std::endl;
}

void xyRBMUNIT::Set_train_para(const std::vector<double> &para)
{
    unsigned int size=para.size();
    if(size!=5){
      std::cout << "para input size incorrect, will use default." << std::endl;
      return;
    }
    epsilonw=para.at(0);
    epsilonvb=para.at(1);
    epsilonhb=para.at(2);
    weightcost=para.at(3);
    momentum=para.at(4);
}

int xyRBMUNIT::rbm_train_one_batch(const std::vector<std::vector<double>> &data)
{
    std::mt19937 rng;
    rng.seed(std::random_device()());
    std::uniform_real_distribution<double> uni_dist(0., 1.);

    unsigned int numcases=data.size();
    unsigned int numdims=0;
    if(numcases>0)numdims=data.at(0).size();

    unsigned int vsize=xyvneurons.size();
    if(vsize!=numdims){
      std::cout << "Data dimension does not match neuron layer size." << std::endl;
      return 1;
    }

    unsigned int hsize=xyhneurons.size();
    if(hsize<1){
      std::cout << "Hid neuron layer size < 1." << std::endl;
      return 1;
    }

    double d,vh,b,A,P;
    for(unsigned int j=0;j<hsize;j++){
      xyhneurons.at(j).Set_size(numcases);
      for(unsigned int i=0;i<numcases;i++){
        A=0.;
        for(unsigned int k=0;k<vsize;k++){
          d=data.at(i).at(k);
          vh=vishid_w.at(k).at(j);
          A=A+d*vh;
          // std::cout << "d,vh,A " << d << " " << vh << " " << A << std::endl;
        }
        b=xyhneurons.at(j).Get_bias();
        A=A+b;
        P=1./(1.+std::exp(-A));
        // std::cout << "A,P: " << A << " | " << P << std::endl;
        xyhneurons.at(j).Set_p_prob(i,P);
      }
    }

    std::cout << "poshidprobs" << std::endl;
    for(unsigned int j=0;j<hsize;j++){
      for(unsigned int i=0;i<numcases;i++){
        std::cout << xyhneurons.at(j).Get_p_prob(i) << " ";
      }
      std::cout << std::endl;
    }
    

    for(unsigned int k=0;k<vsize;k++){
      for(unsigned int j=0;j<hsize;j++){
        A=0.;
        for(unsigned int i=0;i<numcases;i++){
          d=data.at(i).at(k);
          vh=xyhneurons.at(j).Get_p_prob(i);
          A=A+d*vh;
        }
        posprods.at(k).at(j)=A;
      }
    }

    std::cout << "posprods" << std::endl;
    for(unsigned int k=0;k<vsize;k++){
      for(unsigned int j=0;j<hsize;j++){
        std::cout << posprods.at(k).at(j) << " ";
      }
      std::cout << std::endl;
    }

    std::vector<double> poshidact;
    poshidact.clear();
    poshidact.resize(hsize);
    for(unsigned int j=0;j<hsize;j++){
      A=0.;
      for(unsigned int i=0;i<numcases;i++){
        vh=xyhneurons.at(j).Get_p_prob(i);
        A=A+vh;
      }
      poshidact.at(j)=A;
    }

    std::cout << "poshidact" << std::endl;
    for(unsigned int j=0;j<hsize;j++){
      std::cout << poshidact.at(j) << " ";
    }
    std::cout << std::endl;


    std::vector<double> posvisact;
    posvisact.clear();
    posvisact.resize(vsize);
    for(unsigned int k=0;k<vsize;k++){
      A=0.;
      for(unsigned int i=0;i<numcases;i++){
        vh=data.at(i).at(k);
        A=A+vh;
      }
      posvisact.at(k)=A;
    }

    std::cout << "posvisact" << std::endl;
    for(unsigned int k=0;k<vsize;k++){
      std::cout << posvisact.at(k) << " ";
    }
    std::cout << std::endl;

    // std::vector<std::vector<double>> poshidstates;
    poshidstates.clear();

    std::vector<double> tempv;
    double tempr,tempp,tempbi;
    for(unsigned int i=0;i<numcases;i++){
      tempv.clear();
      for(unsigned int j=0;j<hsize;j++){
        tempr=uni_dist(rng);
        tempp=xyhneurons.at(j).Get_p_prob(i);
        tempbi=0.;
        if(tempp>tempr)tempbi=1.;
        tempv.push_back(tempbi);
        std::cout << "tempr,p,bi: " << tempr << " " << tempp << " " << tempbi << std::endl;
      }
      poshidstates.push_back(tempv);
    }

    std::cout << "poshidstates" << std::endl;
    for(unsigned int i=0;i<numcases;i++){
      for(unsigned int j=0;j<hsize;j++){
        std::cout << poshidstates.at(i).at(j) << " ";
      }
      std::cout << std::endl;
    }


    std::vector<std::vector<double>> negdata;
    negdata.clear();
    for(unsigned int i=0;i<numcases;i++){
      tempv.clear();
      for(unsigned int k=0;k<vsize;k++){
        A=0.;
        for(unsigned int j=0;j<hsize;j++){
          vh=vishid_w.at(k).at(j);
          d=poshidstates.at(i).at(j);
          A=A+vh*d;
        }
        b=xyvneurons.at(k).Get_bias();
        A=A+b;
        P=1./(1.+std::exp(-A));
        tempv.push_back(P);
      }
      negdata.push_back(tempv);
    }

    std::cout << "negdata" << std::endl;
    for(unsigned int i=0;i<numcases;i++){
      for(unsigned int k=0;k<vsize;k++){
        std::cout << negdata.at(i).at(k) << " ";
      }
      std::cout << std::endl;
    }


    for(unsigned int j=0;j<hsize;j++){
      for(unsigned int i=0;i<numcases;i++){
        A=0.;
        for(unsigned int k=0;k<vsize;k++){
          d=negdata.at(i).at(k);
          vh=vishid_w.at(k).at(j);
          A=A+d*vh;
        }
        b=xyhneurons.at(j).Get_bias();
        A=A+b;
        P=1./(1.+std::exp(-A));
        xyhneurons.at(j).Set_n_prob(i,P);
      }
    }

    std::cout << "neghidprobs" << std::endl;
    for(unsigned int j=0;j<hsize;j++){
      for(unsigned int i=0;i<numcases;i++){
        std::cout << xyhneurons.at(j).Get_n_prob(i) << " ";
      }
      std::cout << std::endl;
    }


    for(unsigned int k=0;k<vsize;k++){
      for(unsigned int j=0;j<hsize;j++){
        A=0.;
        for(unsigned int i=0;i<numcases;i++){
          d=negdata.at(i).at(k);
          vh=xyhneurons.at(j).Get_n_prob(i);
          A=A+d*vh;
        }
        negprods.at(k).at(j)=A;
      }
    }

    std::cout << "negprods" << std::endl;
    for(unsigned int k=0;k<vsize;k++){
      for(unsigned int j=0;j<hsize;j++){
        std::cout << negprods.at(k).at(j) << " ";
      }
      std::cout << std::endl;
    }


    std::vector<double> neghidact;
    neghidact.clear();
    neghidact.resize(hsize);
    for(unsigned int j=0;j<hsize;j++){
      A=0.;
      for(unsigned int i=0;i<numcases;i++){
        vh=xyhneurons.at(j).Get_n_prob(i);
        A=A+vh;
      }
      neghidact.at(j)=A;
    }


    std::cout << "neghidact" << std::endl;
    for(unsigned int j=0;j<hsize;j++){
      std::cout << neghidact.at(j) << " ";
    }
    std::cout << std::endl;


    std::vector<double> negvisact;
    negvisact.clear();
    negvisact.resize(vsize);
    for(unsigned int k=0;k<vsize;k++){
      A=0.;
      for(unsigned int i=0;i<numcases;i++){
        vh=negdata.at(i).at(k);
        A=A+vh;
      }
      negvisact.at(k)=A;
    }


    std::cout << "negvisact" << std::endl;
    for(unsigned int k=0;k<vsize;k++){
      std::cout << negvisact.at(k) << " ";
    }
    std::cout << std::endl;

    std::cout << std::endl;
    std::cout << "epsilonw,epsilonvb,epsilonhb,weightcost,momentum: " << epsilonw << " " << epsilonvb << " " << epsilonhb << " " << weightcost << " " << momentum << std::endl;

    err=0.;
    for(unsigned int i=0;i<numcases;i++){
      for(unsigned int k=0;k<vsize;k++){
        err=err+std::pow(negdata.at(i).at(k)-data.at(i).at(k),2);
      }
    }
    std::cout << "err= " << err << std::endl;


    for(unsigned int k=0;k<vsize;k++){
      d=xyvneurons.at(k).Get_biasinc();
      vh=posvisact.at(k)-negvisact.at(k);
      A=momentum*d+epsilonvb*vh/numcases;

      xyvneurons.at(k).Set_biasinc(A);
    }

    std::cout << "visbiasinc" << std::endl;
    for(unsigned int k=0;k<vsize;k++){
      std::cout << xyvneurons.at(k).Get_biasinc() << std::endl;
    }

    for(unsigned int j=0;j<hsize;j++){
      d=xyhneurons.at(j).Get_biasinc();
      vh=poshidact.at(j)-neghidact.at(j);
      A=momentum*d+epsilonhb*vh/numcases;

      xyhneurons.at(j).Set_biasinc(A);
    }

    std::cout << "hidbiasinc" << std::endl;
    for(unsigned int j=0;j<hsize;j++){
      std::cout << xyhneurons.at(j).Get_biasinc() << std::endl;
    }


    for(unsigned int k=0;k<vsize;k++){
      for(unsigned int j=0;j<hsize;j++){
        d=vishid_inc.at(k).at(j);
        vh=posprods.at(k).at(j)-negprods.at(k).at(j);
        b=vishid_w.at(k).at(j);
        A=momentum*d+epsilonw*vh/numcases-weightcost*b;

        vishid_inc.at(k).at(j)=A;
      }
    }

    std::cout << "vishidinc" << std::endl;
    for(unsigned int k=0;k<vsize;k++){
      for(unsigned int j=0;j<hsize;j++){
        std::cout << vishid_inc.at(k).at(j) << " ";
      }
      std::cout << std::endl;
    }


    //update
    for(unsigned int k=0;k<vsize;k++){
      d=xyvneurons.at(k).Get_biasinc();
      A=xyvneurons.at(k).Get_bias();

      A=A+d;
      xyvneurons.at(k).Set_bias(A);
    }

    for(unsigned int j=0;j<hsize;j++){
      d=xyhneurons.at(j).Get_biasinc();
      A=xyhneurons.at(j).Get_bias();

      A=A+d;
      xyhneurons.at(j).Set_bias(A);
    }

    for(unsigned int k=0;k<vsize;k++){
      for(unsigned int j=0;j<hsize;j++){
        d=vishid_inc.at(k).at(j);
        A=vishid_w.at(k).at(j);

        A=A+d;
        vishid_w.at(k).at(j)=A;
      }
    }


    std::cout<< "updated visbias" << std::endl;
    for(unsigned int k=0;k<vsize;k++){
      std::cout << xyvneurons.at(k).Get_bias() << " ";
    }
    std::cout << std::endl;

    std::cout<< "updated hidbias" << std::endl;
    for(unsigned int j=0;j<hsize;j++){
      std::cout << xyhneurons.at(j).Get_bias() << " ";
    }
    std::cout << std::endl;

    std::cout<< "updated vishid_w" << std::endl;
    for(unsigned int k=0;k<vsize;k++){
      for(unsigned int j=0;j<hsize;j++){
        std::cout << vishid_w.at(k).at(j) << " ";
      }
      std::cout << std::endl;
    }


    return 0;
}
