
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

#include "stdlib.h"
#include "TROOT.h"
#include "TApplication.h"
#include "Rtypes.h"
#include "math.h"
#include "TFile.h"
#include "TObject.h"
#include "TKey.h"

#include "TChain.h"
#include "TString.h"
#include "TCut.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TEventList.h"
#include "TMath.h"

int main(Int_t argc, char *argv[])
{
  if(argc!=2){
    cout << "No input file when running. Stop." << endl;
    return 0;
  }

  int nrec_out=50;  //show progress every nrec_out events

  TString filename;
  if(argc==2)filename = argv[1];
  cout << "File: " << filename << endl;

  TChain *Tin = new TChain("T","T");
  Tin->AddFile(filename);

  int nev=Tin->GetEntries();
  cout << "Total # of events= " << nev << endl;

  double ID[2157], E[2157];
  UInt_t Nhits,id;
  unsigned int idu;

  double ClAngle[100],angemax;
  double ClEnergy[100],engmax;
  UInt_t Nclusters;
  Int_t EvNb;
  double Good;

  Tin->SetBranchAddress("ID",&ID);
  Tin->SetBranchAddress("Energy",&E);
  Tin->SetBranchAddress("Nhits",&Nhits);
  Tin->SetBranchAddress("Nclusters",&Nclusters);
  Tin->SetBranchAddress("ClAngle",&ClAngle);
  Tin->SetBranchAddress("ClEnergy",&ClEnergy);
  Tin->SetBranchAddress("Good",&Good);
  Tin->SetBranchAddress("EvNb",&EvNb);

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
  cout << "AE1 Read Weights" << endl;
  string prew="W_MLP";
  xyA.read_W(prew);


  xy_Autoencoder xyA2;
  cout << endl;
  cout << "AE2 Read Weights" << endl;
  prew="W_MLP2";
  xyA2.read_W(prew);

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


  //cosmic Nhit>1
  infname="xy_PRad_2200MeV_evnets_cosmic_angmax_gt_5_eclmax_1800_2500_Nhit_gt_1_81_0_61559.txt";
  Dx=81;
  nbat=1;
  tempv2={1.,0.};
  ix=80;
  iy=873;


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


  // return 0;

  // cout << endl;
  // cout << "Targ out together." << endl;

  xy_W_bag tbd;
  xy_W_bag tbo;
  // tbd.W_in(tempm);  //put input data in a bag
  // tbo=xyA.calc_target(tbd);  //use Autoencoder to calculate target output

  vector<vector<double>> targ_out;
  // targ_out=tbo.Get_W();
  // show_matrix(targ_out);


  cout << endl;
  cout << "Targ out" << endl;

  vector<vector<double>> dm_in;
  vector<vector<double>> dm_out;

  // unsigned int sized=tempm.size();
  // sized=2;

  // unsigned int nodis=50;

  targ_out.clear();
  dm_in.clear();
  dm_in.resize(1);
  /*for(unsigned int i=0;i<sized;i++){
    if(i%nodis==0){
      cout << "Calculating ev " << i << " | " << sized << endl;
    }

    dm_in.at(0)=tempm.at(i);
    tbd.W_in(dm_in);
    tbo=xyA.calc_target(tbd);
    dm_out=tbo.Get_W();
    targ_out.push_back(dm_out.at(0));
  }
  */

  // show_matrix(targ_out);

  /*string outname="xy_target_out.txt";
  tbo.W_in(targ_out);
  xyA.print_to_file(tbo,outname);
  */


  double nn_id,nn_id2;

  TFile *xyf = new TFile("PRad_ID.root","RECREATE");
  TTree *Tid = new TTree("Tid","Tid");
  Tid->Branch("nn_id",&nn_id,"nn_id/D");
  Tid->Branch("nn_id2",&nn_id2,"nn_id2/D");
  Tid->Branch("EvNb",&EvNb,"EvNb/I");
  Tid->Branch("Nclusters",&Nclusters,"Nclusters/i");
  Tid->Branch("Nhits",&Nhits,"Nhits/i");
  Tid->Branch("Good",&Good,"Good/D");
  Tid->Branch("engmax",&engmax,"engmax/D");
  Tid->Branch("angemax",&angemax,"angemax/D");


  int emax;
  int ni=0;
  int ncount=0;

  tempv.clear();
  tempv.resize(Dy);
  while(ni<nev){
    Tin->GetEntry(ni);

    ni=ni+1;

    if(Nclusters<1){
      continue;
    }

    emax=TMath::LocMax(Nclusters,&ClEnergy[0]);
    angemax=ClAngle[emax];
    engmax=ClEnergy[emax];

    if(angemax<5.)continue;
    if(engmax<1800.)continue;
    if(engmax>2500.)continue;

    for(int j=0;j<1728;j++){
      if(E[j]>1999.)E[j]=1999.;
      E[j]=E[j]/2000.;
    }


    for(unsigned int j=0;j<2157;j++){
      tempv.at(j)=0.;
    }

    for(UInt_t k=0;k<Nhits;k++){
      id=ID[k];
      idu=id;
      tempv.at(idu)=E[k];
    }

    dm_in.at(0)=tempv;
    tbd.W_in(dm_in);

    tbo=xyA.calc_target(tbd);
    dm_out=tbo.Get_W();
    nn_id=dm_out.at(0).at(1)-dm_out.at(0).at(0);

    tbo=xyA2.calc_target(tbd);
    dm_out=tbo.Get_W();
    nn_id2=dm_out.at(0).at(1)-dm_out.at(0).at(0);

    // cout << nn_id << endl;
    if(nn_id!=nn_id)nn_id=0.;
    if(nn_id2!=nn_id2)nn_id2=0.;

    if(ncount%nrec_out==0){
      cout << "Calculating: " << ni << " | " << nev << " | ncount= " << ncount << " | nn_id= " << nn_id << " | nn_id2= " << nn_id2 << endl;
    }

    Tid->Fill();
    ncount=ncount+1;
  }
  cout << "ncount= " << ncount << endl;
  xyf->Write();

  return 0;


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
