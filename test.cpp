
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
  xy_W_bag b_check;

  xy_Autoencoder xyA;


  cout << endl;
  cout << "Read Weights" << endl;
  string prew="W_MLP";
  xyA.read_W(prew);

  /*cout << endl;
  cout << "Print W for test" << endl;
  double pty=1.;
  prew="W_test";
  xyA.print_W(pty,prew);
  */

  ifstream infile;
  string infname;
  double temp_d;
  unsigned int Dx;
  unsigned int Dy;
  unsigned int ndis=5000;  //display reading process

  Dy=2157;

  unsigned int ix;
  unsigned int iy;

  //sim train
  /*infname="xy_PRad_2200MeV_evnets_sim_angmax_gt_5_eclmax_1800_2500_1k_0_373069.txt";
  Dx=1000;
  nbat=10;
  tempv2={0.,1.};
  ix=999;
  iy=718;
  */

  //sim test
  /*infname="xy_PRad_2200MeV_evnets_sim_angmax_gt_5_eclmax_1800_2500_500_374000_563649.txt";
  Dx=500;
  nbat=5;
  tempv2={0.,1.};
  ix=499;
  iy=2154;
  */

  //GEM HyCal match good
  /*infname="xy_PRad_2200MeV_evnets_good_angmax_gt_5_eclmax_1800_2500_500_3016970_4535147.txt";
  Dx=500;
  nbat=5;
  tempv2={0.,1.};
  ix=499;
  iy=86;
  */

  //cosmic file
  infname="xy_PRad_2200MeV_evnets_cosmic_angmax_gt_5_eclmax_1800_2500_107_0_61559.txt";
  Dx=107;
  nbat=1;
  tempv2={1.,0.};
  ix=106;
  iy=875;
  


  cout << endl;
  cout << "Read test file" << endl;
  tempm.clear();
  temp_targ.clear();
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

  cout << tempm.at(ix).at(iy) << endl;
  cout << temp_targ.at(ix).at(0) << " " << temp_targ.at(ix).at(1) << endl;

  xyA.test_data_in(tempm,temp_targ,nbat);


  vector<unsigned int> test_nout;
  test_nout=xyA.n_test_miss();
  std::cout << "Test: totmis= " << test_nout.at(1) << " , out of " << test_nout.at(0) << endl;


  cout << endl;
  cout << "Targ out" << endl;

  xy_W_bag tbd;
  xy_W_bag tbo;
  tbd.W_in(tempm);  //put input data in a bag
  tbo=xyA.calc_target(tbd);  //use Autoencoder to calculate target output

  vector<vector<double>> targ_out;
  targ_out=tbo.Get_W();
  show_matrix(targ_out);

  cout << "missed index" << endl;
  unsigned int ntout=targ_out.size();
  for(unsigned int i=0;i<ntout;i++){
    tempv=targ_out.at(i);
    if(tempv2.at(0)>tempv2.at(1)){
      if(tempv.at(0)<tempv.at(1)){
        cout << i << " ";
      }
    }else{
      if(tempv.at(0)>tempv.at(1)){
        cout << i << " ";
      }
    }
  }
  cout << endl;


  return 0;

}
