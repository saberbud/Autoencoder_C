#include "xy_Autoencoder.h"
#include "xy_utils.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <cmath>


// constructor
xy_Autoencoder::xy_Autoencoder()
{   
    Weights.clear();
    train_wps.clear();
    test_wps.clear();
    rbmv.clear();
}

void xy_Autoencoder::train_data_in(const std::vector<std::vector<double>> &data, const std::vector<std::vector<double>> &targ, unsigned int &nb)
{
    std::cout << "xy_A:train_data_in" << std::endl;

    unsigned int size=data.size();
    std::cout << "Training data cases= " << size << std::endl;

    if(size!=targ.size()){
      std::cout << "Data and target have different # of cases. STOP." << std::endl;
      return;
    }

    train_data.data_in(data);
    train_data.targ_in(targ);

    train_data.data_shuffle();

    /*std::cout << "Data in train_data: " << std::endl;
    show_matrix(train_data.Get_data());

    std::cout << "Target in train_data: " << std::endl;
    show_matrix(train_data.Get_targ());
    */

    train_data.make_mini_batch(nb);

    std::cout << "mini batch data and target" << std::endl;
    for(unsigned int i=0;i<nb;i++){
      std::cout << std::endl;
      std::cout << "Data:" << std::endl;
      show_matrix(train_data.Get_mini_batch(i).Get_W());
      std::cout << "Target:" << std::endl;
      show_matrix(train_data.Get_mini_batch_t(i).Get_W());
    }

}

void xy_Autoencoder::test_data_in(const std::vector<std::vector<double>> &data, const std::vector<std::vector<double>> &targ, unsigned int &nb)
{
    std::cout << "xy_A:test_data_in" << std::endl;
    
    unsigned int size=data.size();
    std::cout << "Test data cases= " << size << std::endl;

    if(size!=targ.size()){
      std::cout << "Data and target have different # of cases. STOP." << std::endl;
      return;
    }

    test_data.data_in(data);
    test_data.targ_in(targ);

    /*std::cout << "Data in test_data: " << std::endl;
    show_matrix(test_data.Get_data());

    std::cout << "Target in test_data: " << std::endl;
    show_matrix(test_data.Get_targ());
    */

    test_data.make_mini_batch(nb);

    std::cout << "mini batch data and target" << std::endl;
    for(unsigned int i=0;i<nb;i++){
      std::cout << std::endl;
      std::cout << "Data:" << std::endl;
      show_matrix(test_data.Get_mini_batch(i).Get_W());
      std::cout << "Target:" << std::endl;
      show_matrix(test_data.Get_mini_batch_t(i).Get_W());
    }

}

void xy_Autoencoder::Set_net_strucuture(const std::vector<unsigned int> &st)
{
    unsigned int nl=st.size();
    if(nl<2){
      std::cout << "Set net_structure input vector has wrong dimension." << std::endl;
      return;
    }

    std::cout << "# of layers, including data and target_out: " << nl << std::endl;

    Weights.resize(nl);
    // train_wps.resize(nl);
    // test_wps.resize(nl);

    std::cout << "Weights, train_wps, test_wps size: " << Weights.size() << " " << train_wps.size() << " " << test_wps.size() << std::endl;

    unsigned int sizep=nl-1;

    unsigned int nx=st.at(sizep-1);
    unsigned int ny=st.at(sizep);
    Weights.at(sizep).W_in(rand_mat(nx,ny));
    Weights.at(sizep).b_in(rand_vec(ny));

    std::vector<double> tempv;
    std::vector<std::vector<double>> tempm;
    std::cout << "Weight " << sizep << std::endl;
    tempm=Weights.at(sizep).Get_W();
    show_matrix(tempm);
    tempv=Weights.at(sizep).Get_b();
    tempm.clear();
    tempm.push_back(tempv);
    show_matrix(tempm);


    rbmv.resize(sizep);

    for(unsigned int i=1;i<sizep;i++){
      rbmv.at(i).Set_v_size(st.at(i-1));
      rbmv.at(i).Set_h_size(st.at(i));
    }

    for(unsigned int i=1;i<sizep;i++){
      std::cout << "RBM " << i << " vsize= " << rbmv.at(i).GetvSize() << " , hsize= " << rbmv.at(i).GethSize() << std:: endl;
    }

}

void xy_Autoencoder::init_RBMs()
{
    unsigned int sizer=rbmv.size();
    if(sizer<2){
      std::cout << "Structure of RBMs not set yet. Do not init yet." << std::endl;
      return;
    }

    for(unsigned int ir=1;ir<sizer;ir++){
      rbmv.at(ir).init_vh_w();
      rbmv.at(ir).init_v_bias();
      rbmv.at(ir).init_h_bias();
      rbmv.at(ir).init_prods();
      rbmv.at(ir).init_vh_inc();
    }

    std::vector<double> tempv;
    std::vector<std::vector<double>> tempm;

    for(unsigned int ir=1;ir<sizer;ir++){
      std::cout << "RBM " << ir << std::endl;
      tempm=rbmv.at(ir).Get_vishid();
      show_matrix(tempm);

      tempm.clear();
      tempv=rbmv.at(ir).Get_h_bias();
      tempm.push_back(tempv);
      show_matrix(tempm);
      std::cout << std::endl;
    }

}

void xy_Autoencoder::Set_RBM_para(std::vector<double> &para, unsigned int &id)
{
    unsigned int sizer=rbmv.size();
    if((id+1)>sizer){
      std::cout << "Set_RBM_para: index out of RBM set size." << std::endl;
      return;
    }

    rbmv.at(id).Set_train_para(para);
}

int xy_Autoencoder::train_RBMs(unsigned int &nep, std::vector<unsigned int> &sw, std::vector<std::vector<double>> &swp)
{
    unsigned int sizer=rbmv.size();
    unsigned int sizenb=train_data.Get_nbatch();
    train_wps.resize(sizenb);

    std::cout << "para switch times= " << sw.size() << std::endl;
    show_matrix(swp);

    std::cout << "# of epochs= " << nep << std::endl;
    std::cout << "# of RBMs: " << (sizer-1) << " , # of batches= " << sizenb << " , # of wps bags= " << train_wps.size() << std::endl;

    int trbm=0;

    unsigned int nsw;

    // sizer=3;  //for debug or development
    for(unsigned int nr=1;nr<sizer;nr++){ //RBM loop
      nsw=0;
      for(unsigned int ne=0;ne<nep;ne++){  //epoch loop
        if(nsw<sw.size()){
          if(ne==sw.at(nsw)){
            std::cout << "switch paras" << std::endl;
            Set_RBM_para(swp.at(nsw),nr);
            nsw=nsw+1;
          }
        }
        for(unsigned int nb=0;nb<sizenb;nb++){  //mini batch loop
          std::cout << "RBM " << nr << " , epoch " << (ne+1) << " , batch " << nb << std::endl;

          if(nr==1){
            trbm=rbmv.at(nr).rbm_train_one_batch(train_data.Get_mini_batch(nb).Get_W());
          }else{
            trbm=rbmv.at(nr).rbm_train_one_batch(train_wps.at(nb).Get_W());
          }

          // std::cout << "trbm= " << trbm << std::endl;
          if(trbm>0){
            return 1;
          }

          if(ne==(nep-1)){
            train_wps.at(nb).W_in(rbmv.at(nr).Get_poshidprobs());
          }

        }  //mini batch loop
      }  //epoch loop
      Weights.at(nr).W_in(rbmv.at(nr).Get_vishid());
      Weights.at(nr).b_in(rbmv.at(nr).Get_h_bias());
    } //RBM loop


    unsigned int tempi=99;
    std::cout << "RBM weights" << std::endl;
    show_ws(Weights,tempi); 

    std::cout << std::endl;
    std::cout << "RBM probs" << std::endl;
    tempi=0;
    show_ws(train_wps,tempi);
    tempi=99;
    show_ws(train_wps,tempi);

    return 0;
}

void xy_Autoencoder::Set_mlp_para(std::vector<double> &para)
{
    mlp.init_Ps(para);
}

void xy_Autoencoder::train_mlp(unsigned int &nep, unsigned int &n_o, double &leng)
{
    std::cout << "Total epoch= " << nep << " , 1-layer epoch= " << n_o << std::endl;

    mlp.input_weights(Weights);

    std::vector<xy_W_bag> ws_c;
    ws_c=mlp.Get_weights();
    unsigned int sizewc=ws_c.size();
    std::cout << "Weights in mlp before training" << std::endl;
    std::cout << "Weight # in mlp: " << sizewc << std::endl;
    show_ws(ws_c,sizewc);
    std::cout << std::endl;

    double xy_type=0.;  //0: do w_class; 1: do all
    double xy_leng=leng;

    std::vector<std::vector<double>> tempm1;
    std::vector<std::vector<double>> tempm2;
    unsigned int totmis=0;
    unsigned int totdim=0;

    std::vector<unsigned int> tnout;

    xy_W_bag tempbag;
    unsigned int sizenb=train_data.Get_nbatch();
    std::cout << "# of batches: " << sizenb << std::endl;
    for(unsigned int ne=0;ne<nep;ne++){  //epoch loop
      if(ne==n_o)xy_type=1.;

      totmis=0;
      totdim=0;

      if(xy_type<0.5){
        std::cout << "Class layer fine tuning. # search= " << xy_leng << std::endl;
      }else{
        std::cout << "All layers fine tuning. # search= " << xy_leng << std::endl;
      }

      for(unsigned int nb=0;nb<sizenb;nb++){  //mini batch loop
        std::cout << "Epoch: " << ne << " , batch " << nb << std::endl;

        tempbag=train_data.Get_mini_batch(nb);
        mlp.input_data(tempbag);
        tempbag=train_data.Get_mini_batch_t(nb);
        mlp.input_targ(tempbag);

        mlp.eval_probs_base();
        mlp.eval_probs_targ();
        mlp.eval_d_class();
        mlp.eval_d_base();

        mlp.xy_minimize(xy_type,xy_leng);

        ws_c.clear();
        ws_c=mlp.Get_probs();
        sizewc=ws_c.size()-1;
        // show_ws(ws_c,sizewc);
        // sizewc=ws_c.size();
        // show_ws(ws_c,sizewc);

        tempm1=train_data.Get_mini_batch_t(nb).Get_W();
        tempm2=ws_c.at(sizewc).Get_W();
        totmis=totmis+n_miss_match(tempm1,tempm2);
        totdim=totdim+tempm1.size();

      }  //mini batch loop

      std::cout << "totmis= " << totmis << " , out of " << totdim << std::endl;
      tnout=n_test_miss();
      std::cout << "Test: totmis= " << tnout.at(1) << " , out of " << tnout.at(0) << std::endl;
      std::cout << std::endl;
    }  //epoch loop


    std::cout << std::endl;
    std::cout << "Weights in mlp after training" << std::endl;
    ws_c.clear();
    ws_c=mlp.Get_weights();
    sizewc=ws_c.size();
    std::cout << "Weight # in mlp: " << sizewc << std::endl;
    show_ws(ws_c,sizewc);


}

void xy_Autoencoder::set_mlp_w_W()
{
    mlp.input_weights(Weights);
}

std::vector<unsigned int> xy_Autoencoder::n_test_miss(){
    std::vector<unsigned int> P;
    P.clear();

    unsigned int ntot=0;
    unsigned int nmis=0;

    P.push_back(ntot);
    P.push_back(nmis);

    unsigned int sizenb=test_data.Get_nbatch();

    /*std::cout << "test data # batch= " << sizenb << std::endl;
    std::cout << "mini batch data and target" << std::endl;
    for(unsigned int i=0;i<sizenb;i++){
      std::cout << std::endl;
      std::cout << "Data:" << std::endl;
      show_matrix(test_data.Get_mini_batch(i).Get_W());
      std::cout << "Target:" << std::endl;
      show_matrix(test_data.Get_mini_batch_t(i).Get_W());
    }
    */

    xy_W_bag tbd;
    xy_W_bag tb1;
    xy_W_bag tb2;

    std::vector<std::vector<double>> tempm1;
    std::vector<std::vector<double>> tempm2;

    for(unsigned int i=0;i<sizenb;i++){
      tbd=test_data.Get_mini_batch(i);  //data
      tb1=test_data.Get_mini_batch_t(i);  //original target
      tb2=calc_target(tbd);  //calculated target

      // std::cout << "calc target:" << i << std::endl;
      // show_matrix(tb2.Get_W());

      ntot=ntot+tb2.Get_W().size();
      nmis=nmis+n_miss_match(tb1.Get_W(),tb2.Get_W());
    }

    P.at(0)=ntot;
    P.at(1)=nmis;

    return P;
}

xy_W_bag xy_Autoencoder::calc_target(xy_W_bag &A){
    xy_W_bag P;

    mlp.input_data(A);
    mlp.eval_probs_base();
    mlp.eval_probs_targ();

    std::vector<xy_W_bag> tempbs;
    tempbs.clear();
    tempbs=mlp.Get_probs();
    unsigned int id=tempbs.size()-1;
    P=tempbs.at(id);

    return P;
}

unsigned int xy_Autoencoder::n_miss_match(const std::vector<std::vector<double>> &A, const std::vector<std::vector<double>> &B)
{
    unsigned int nm=0;

    unsigned int size=A.size();
    if(size!=B.size()){
      std::cout << "n_miss_match: input matrices have different dimension." << std::endl;
      return nm;
    }

    unsigned int ta;
    unsigned int tb;

    for(unsigned int i=0;i<size;i++){
      ta=max_idx(A.at(i));
      tb=max_idx(B.at(i));

      if(ta!=tb){
        nm=nm+1;
      }
    }

    return nm;
}

void xy_Autoencoder::show_ws(std::vector<xy_W_bag> &A, unsigned int &id)
{
    unsigned int size=A.size();

    std::vector<double> tempv;
    std::vector<std::vector<double>> tempm;

    if(id<size){
      std::cout << "W " << id << std::endl;
      tempm=A.at(id).Get_W();
      show_matrix(tempm);
      std::cout << "b " << id << std::endl;
      tempm.clear();
      tempv=A.at(id).Get_b();
      tempm.push_back(tempv);
      show_matrix(tempm);
      return;
    }

    for(unsigned int i=1;i<size;i++){
      std::cout << "W " << i << std::endl;
      tempm=A.at(i).Get_W();
      show_matrix(tempm);
      std::cout << "b " << i << std::endl;
      tempm.clear();
      tempv=A.at(i).Get_b();
      tempm.push_back(tempv);
      show_matrix(tempm);
    }

}

void xy_Autoencoder::Set_one_W(xy_W_bag &A, unsigned int &id)
{
    unsigned int sizew=Weights.size();
    unsigned int sizep=sizew-1;

    if(id>sizew || id>sizep){
      std::cout << "Set_one_W input id out of Weights range." << std::endl;
      return;
    }

    Weights.at(id)=A;

}

unsigned int xy_Autoencoder::max_idx(const std::vector<double> &A)
{
    unsigned int id=0;

    unsigned int size=A.size();
    for(unsigned int i=1;i<size;i++){
      if(A.at(i)>A.at(i-1))id=i;
    }

    return id;
}

std::vector<double> xy_Autoencoder::rand_vec(unsigned int &ny)
{
    std::vector<double> P;
    P.clear();

    std::mt19937 rng;
    rng.seed(std::random_device()());
    std::normal_distribution<double> norm_dist(0., 1.);

    double temp;
    for(unsigned int i=0;i<ny;i++){
      temp=norm_dist(rng);
      temp=temp*0.1;
      P.push_back(temp);
    }

    return P;
}

std::vector<std::vector<double>> xy_Autoencoder::rand_mat(unsigned int &nx, unsigned int &ny)
{
    std::vector<std::vector<double>> P;
    P.clear();

    std::mt19937 rng;
    rng.seed(std::random_device()());
    std::normal_distribution<double> norm_dist(0., 1.);

    std::vector<double> tempv;
    double temp;

    for(unsigned int j=0;j<nx;j++){
      tempv.clear();
      for(unsigned int i=0;i<ny;i++){
        temp=norm_dist(rng);
        temp=temp*0.1;
        tempv.push_back(temp);
      }
      P.push_back(tempv);
    }

    return P;
}

void xy_Autoencoder::print_W(double &type, std::string &pref)
{
    std::vector<xy_W_bag> tbv;

    if(type<0.5){
      tbv=Weights;
    }else{
      tbv=mlp.Get_weights();
    }

    unsigned int sizew=tbv.size();

    std::cout << "W number file:" << std::endl;
    std::ofstream outfile(pref);
    outfile << sizew << std::endl;
    outfile.close();

    std::string fname;
    char buffer[50];
    for(unsigned int i=1;i<sizew;i++){
      std::sprintf(buffer, "_%d", i);
      fname=buffer;
      fname=pref+fname;
      // fname=pref+std::Form("_%d",i);
      // std::cout << "fname: " << fname << std::endl;
      print_to_file(tbv.at(i),fname);
    }

}

void xy_Autoencoder::print_to_file(xy_W_bag &A, std::string &name)
{
    std::cout << "print to file: " << name << std::endl;


    std::vector<double> tempv;
    std::vector<std::vector<double>> tempm;

    tempm=A.Get_W();
    tempv=A.Get_b();

    unsigned int nx=tempm.size();
    unsigned int ny;
    if(nx<1){
      ny=0;
    }else{
      ny=tempm.at(0).size();
    }

    std::ofstream outfile(name);
    outfile << nx << " " << ny << std::endl;

    for(unsigned int i=0;i<nx;i++){
      if(ny!=tempm.at(i).size()){
        std::cout << "row length is not const. Stop." << std::endl;
        return;
      }
      for(unsigned int j=0;j<ny;j++){
        outfile << std::setprecision(11) << tempm.at(i).at(j) << " ";
      }
      outfile << std::endl;
    }

    nx=1;
    ny=tempv.size();
    outfile << nx << " " << ny << std::endl;
    for(unsigned int j=0;j<ny;j++){
      outfile << std::setprecision(11) << tempv.at(j) << " ";
    }
    outfile << std::endl;

    outfile.close();

}
