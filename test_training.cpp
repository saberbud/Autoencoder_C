
#include "include/xy_utils.h"
#include "include/xyNeuron.h"
#include "include/xyrbmunit.h"
#include "include/xydata_mini_batch.h"
#include "include/xy_W_bag.h"
#include "include/xy_MLP.h"
#include "include/xy_Autoencoder.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(int /*argc*/, char * /*argv*/ [])
{

  vector<double> tempv;
  vector<double> tempv2;
  vector<vector<double>> tempm;
  vector<vector<double>> temp_data;
  vector<vector<double>> temp_targ;

  vector<vector<double>> tr_data;
  vector<vector<double>> tr_targ;

  unsigned int nbat=100;
  unsigned int idx=1;
  xy_W_bag b_check;

  xy_Autoencoder xyA;

  cout << endl;
  cout << "Set net structure" << endl;
  std::vector<unsigned int> strc={2157,2200,1000,500,2};
  xyA.Set_net_strucuture(strc);

  ifstream infile;
  // string infname="xy_PRad_2200MeV_evnets_good_20000_50000.txt";
  // string infname="xy_PRad_2200MeV_evnets_good_angmax_gt_5_eclmax_1800_2500_1k_0_3016961.txt";
  string infname="xy_PRad_2200MeV_evnets_sim_angmax_gt_5_eclmax_1800_2500_1k_0_373069.txt";
  unsigned int ndis=5000;  //display reading process

  double temp_d;
  // unsigned int Dx=3000;
  unsigned int Dx=1000;
  unsigned int Dy=2157;
  tempv2={0.,1.};
  infile.open(infname);
  if(!infile.is_open()){
    cout << "Input file not exist: " << infname << endl;
    return 0;
  }

  tempm.clear();
  temp_targ.clear();
  for(unsigned int i=0;i<Dx;i++){
    tempv.clear();
    if(i%ndis==0){
      cout << "Reading row " << i << endl;
    }
    for(unsigned int j=0;j<Dy;j++){
      infile >> temp_d;
      tempv.push_back(temp_d);
    }
    tempm.push_back(tempv);
    temp_targ.push_back(tempv2);
  }
  infile.close();

  // Dx=3000;
  Dx=1000;
  Dy=2157;
  tempv2={1.,0.};
  // infname="xy_PRad_2200MeV_evnets_cosmic_20000_50000.txt";
  infname="xy_PRad_2200MeV_evnets_emptytar_angmax_gt_5_eclmax_1800_2500_1k_0_12121883.txt";
  infile.open(infname);
  if(!infile.is_open()){
    cout << "Input file not exist: " << infname << endl;
    return 0;
  }
  for(unsigned int i=0;i<Dx;i++){
    tempv.clear();
    if(i%ndis==0){
      cout << "Reading row " << i << endl;
    }
    for(unsigned int j=0;j<Dy;j++){
      infile >> temp_d;
      tempv.push_back(temp_d);
    }
    tempm.push_back(tempv);
    temp_targ.push_back(tempv2);
  }
  infile.close();

  unsigned int sdx=tempm.size();
  unsigned int sdy=tempm.at(0).size();
  unsigned int sdy2=temp_targ.at(0).size();
  cout << "Train data size= " << sdx << " X " << sdy << endl;
  for(unsigned int i=0;i<sdx;i++){
    if(sdy!=tempm.at(i).size()){
      cout << "y dimension " << sdy << " != " << tempm.at(i).size() << endl;
    }

    if(sdy2!=temp_targ.at(i).size()){
      cout << "y dimension " << sdy2 << " != " << temp_targ.at(i).size() << endl;
    }
  }

  unsigned int ix=999;
  unsigned int iy=718;
  cout << tempm.at(ix).at(iy) << endl;
  cout << temp_targ.at(ix).at(0) << " " << temp_targ.at(ix).at(1) << endl;

  ix=1999;iy=26;
  cout << tempm.at(ix).at(iy) << endl;
  cout << temp_targ.at(ix).at(0) << " " << temp_targ.at(ix).at(1) << endl;
  

  // show_matrix(temp_targ);

  tr_data=tempm;
  tr_targ=temp_targ;

  // nbat=60;  //3k+3k in
  nbat=20;  //1k+1k
  xyA.train_data_in(tempm,temp_targ,nbat);


  // Dx=1000;
  Dx=500;
  Dy=2157;
  tempv2={0.,1.};
  tempm.clear();
  temp_targ.clear();
  // infname="xy_PRad_2200MeV_evnets_good_0_10000.txt";
  // infname="xy_PRad_2200MeV_evnets_good_angmax_gt_5_eclmax_1800_2500_500_3016970_4535147.txt";
  infname="xy_PRad_2200MeV_evnets_sim_angmax_gt_5_eclmax_1800_2500_500_374000_563649.txt";
  infile.open(infname);
  if(!infile.is_open()){
    cout << "Input file not exist: " << infname << endl;
    return 0;
  }
  for(unsigned int i=0;i<Dx;i++){
    tempv.clear();
    if(i%ndis==0){
      cout << "Reading row " << i << endl;
    }
    for(unsigned int j=0;j<Dy;j++){
      infile >> temp_d;
      tempv.push_back(temp_d);
    }
    tempm.push_back(tempv);
    temp_targ.push_back(tempv2);
  }
  infile.close();


  // Dx=1000;
  Dx=500;
  Dy=2157;
  tempv2={1.,0.};
  // infname="xy_PRad_2200MeV_evnets_cosmic_0_10000.txt";
  infname="xy_PRad_2200MeV_evnets_emptytar_angmax_gt_5_eclmax_1800_2500_500_12121885_17485610.txt";
  infile.open(infname);
  if(!infile.is_open()){
    cout << "Input file not exist: " << infname << endl;
    return 0;
  }
  for(unsigned int i=0;i<Dx;i++){
    tempv.clear();
    if(i%ndis==0){
      cout << "Reading row " << i << endl;
    }
    for(unsigned int j=0;j<Dy;j++){
      infile >> temp_d;
      tempv.push_back(temp_d);
    }
    tempm.push_back(tempv);
    temp_targ.push_back(tempv2);
  }
  infile.close();

  sdx=tempm.size();
  sdy=tempm.at(0).size();
  sdy2=temp_targ.at(0).size();
  cout << "Test data size= " << sdx << " X " << sdy << endl;
  for(unsigned int i=0;i<sdx;i++){
    if(sdy!=tempm.at(i).size()){
      cout << "y dimension " << sdy << " != " << tempm.at(i).size() << endl;
    }

    if(sdy2!=temp_targ.at(i).size()){
      cout << "y dimension " << sdy2 << " != " << temp_targ.at(i).size() << endl;
    }
  }


  ix=499;
  iy=2154;
  cout << tempm.at(ix).at(iy) << endl;
  cout << temp_targ.at(ix).at(0) << " " << temp_targ.at(ix).at(1) << endl;

  ix=999;
  iy=2142;
  cout << tempm.at(ix).at(iy) << endl;
  cout << temp_targ.at(ix).at(0) << " " << temp_targ.at(ix).at(1) << endl;
  

  nbat=10;
  xyA.test_data_in(tempm,temp_targ,nbat);

  //cosmic test
  Dx=107;
  Dy=2157;
  tempv2={1.,0.};
  tempm.clear();
  temp_targ.clear();
  infname="xy_PRad_2200MeV_evnets_cosmic_angmax_gt_5_eclmax_1800_2500_107_0_61559.txt";
  infile.open(infname);
  if(!infile.is_open()){
    cout << "Input file not exist: " << infname << endl;
    return 0;
  }
  for(unsigned int i=0;i<Dx;i++){
    tempv.clear();
    if(i%ndis==0){
      cout << "Reading row " << i << endl;
    }
    for(unsigned int j=0;j<Dy;j++){
      infile >> temp_d;
      tempv.push_back(temp_d);
    }
    tempm.push_back(tempv);
    temp_targ.push_back(tempv2);
  }
  infile.close();

  vector<vector<double>> temp_data_1;
  vector<vector<double>> temp_targ_1;
  temp_data_1=tempm;
  temp_targ_1=temp_targ;

  ix=106;
  iy=875;
  cout << temp_data_1.at(ix).at(iy) << endl;
  cout << temp_targ_1.at(ix).at(0) << " " << temp_targ_1.at(ix).at(1) << endl;

  // nbat=1;
  // xyA.test_data_in(temp_data_1,temp_targ_1,nbat);

  // nbat=10;
  // xyA.train_data_in(tr_data,tr_targ,nbat);

  // return 0;

  cout << endl;
  cout << "RBM init" << endl;
  xyA.init_RBMs();

  cout << endl;
  cout << "RBM set para" << endl;
  // vector<double> para={0.1,0.1,0.1,0.0002,0.9};  //proper size
  vector<double> para={0.1,0.1,0.1,0.0002};  //wrong size
  unsigned int sizer=strc.size()-1;
  for(unsigned int i=1;i<sizer;i++){
    xyA.Set_RBM_para(para,i);
  }


  cout << endl;
  cout << "RBM training" << endl;
  unsigned int nepoch=10;

  vector<unsigned int> swi;
  unsigned int swit=6;
  swi.push_back(swit);

  vector<vector<double>> swpi;
  para.clear();
  para={0.1,0.1,0.1,0.0002,0.9};
  swpi.push_back(para);

  int trbm;
  trbm=xyA.train_RBMs(nepoch,swi,swpi);
  cout << "trbm " << trbm << endl;

  cout << "print RBM W" << endl;
  double ptype=0.;
  string pre="W_RBM";
  xyA.print_W(ptype,pre);

  cout << endl;
  cout << "MLP training" << endl;

  // nbat=6;
  // xyA.train_data_in(tr_data,tr_targ,nbat);

  // tempv={0.1,3.0,20.,10.,0.1,0.05};  //proper length
  tempv={0.1,3.0,20.,10.,0.1};  //wrong size
  xyA.Set_mlp_para(tempv);

  nepoch=30;
  unsigned int n_one=5;
  double lsleng=3.;
  xyA.train_mlp(nepoch,n_one,lsleng);



  nbat=1;
  xyA.test_data_in(temp_data_1,temp_targ_1,nbat);
  vector<unsigned int> test_nout;
  test_nout=xyA.n_test_miss();
  std::cout << "Test cosmic: totmis= " << test_nout.at(1) << " , out of " << test_nout.at(0) << endl;


  cout << "print MLP W" << endl;
  ptype=1.;
  pre="W_MLP";
  xyA.print_W(ptype,pre);

  return 0;


  cout << endl;
  cout << "n miss calc test" << endl;
  unsigned int nmis;
  tempm.clear();
  tempv={0.,1.};
  tempm.push_back(tempv);
  tempv={1.,0.};
  tempm.push_back(tempv);
  tempv={0.,1.};
  tempm.push_back(tempv);
  tempv={1.,0.};
  tempm.push_back(tempv);
  temp_targ=tempm;

  tempm.clear();
  tempv={0.6,0.5};
  tempm.push_back(tempv);
  tempv={1.,0.};
  tempm.push_back(tempv);
  tempv={0.,1.};
  tempm.push_back(tempv);
  tempv={0.1,0.9};
  tempm.push_back(tempv);

  nmis=xyA.n_miss_match(temp_targ,tempm);
  cout << "n miss= " << nmis << endl;


  return 0;

  cout << endl;
  cout << "Max index test" << endl;
  unsigned int size=temp_targ.size();
  for(unsigned int i=0;i<size;i++){
    tempv=temp_targ.at(i);
    idx=xyA.max_idx(tempv);

    for(unsigned int j=0;j<tempv.size();j++){
      cout << tempv.at(j) << " ";
    }
    cout << " | " << idx << endl;
  }


  cout << endl;
  cout << "Rand W_class test" << endl;
  unsigned int nx=2;
  unsigned int ny=2;

  tempm=xyA.rand_mat(nx,ny);
  show_matrix(tempm);

  tempv=xyA.rand_vec(ny);
  tempm.clear();
  tempm.push_back(tempv);
  show_matrix(tempm);


  string fname="output_test";
  // cout << "Output file name= " << fname << endl;
  /*ofstream outfile(fname);
  outfile << "xy" << endl;
  outfile.close();
  */
  // xydm.print_to_file(b_check,fname);
  // xyA.print_to_file(wbagv.at(1),fname);


  return 0;
}
