#include "xy_MLP.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>


// constructor
xy_MLP::xy_MLP()
{
    // place holder
}

xy_W_bag xy_MLP::calc_prob(xy_W_bag &A, xy_W_bag &B)
{
    xy_W_bag P;

    std::vector<std::vector<double>> A_w;
    A_w=A.Get_W();

    std::vector<std::vector<double>> B_w;
    B_w=B.Get_W();

    std::vector<double> B_b;
    B_b=B.Get_b();

    unsigned int nx=A_w.size();
    if(nx<1){
      std::cout << "calc_prob input A matrix size=0." << std::endl;
      return P;
    }

    unsigned int ny=A_w.at(0).size();
    for(unsigned int i=0;i<nx;i++){
      if(ny!=A_w.at(i).size()){
        std::cout << "calc_prob: A matrix row size not same." << std::endl;
        return P;
      }
    }
    // std::cout << "calc_prob A matrix size= " << nx << " X " << ny << std::endl;

    unsigned int nx2=B_w.size();
    if(nx2<1){
      std::cout << "calc_prob: input B matrix size=0." << std::endl;
      return P;
    }

    unsigned int ny2=B_b.size();
    for(unsigned int j=0;j<nx2;j++){
      if(ny2!=B_w.at(j).size()){
        std::cout << "calc_prob: B matrix row size not same." << std::endl;
        return P;
      }
    }
    // std::cout << "calc_prob B matrix size= " << nx2 << " X " << ny2 << std::endl;

    if(ny!=nx2){
      std::cout << "calc_prob: input A B matrix dimension mismatch" << std::endl;
      return P;
    }

    std::vector<double> tempv;
    std::vector<std::vector<double>> tempm;
    tempm.clear();
    double d,vh,b,S,Sp;
    for(unsigned int i=0;i<nx;i++){
      tempv.clear();
      for(unsigned int k=0;k<ny2;k++){
        S=0.;
        for(unsigned int j=0;j<ny;j++){
          d=A_w.at(i).at(j);
          vh=B_w.at(j).at(k);
          S=S+d*vh;
        }
        b=B_b.at(k);
        S=S+b;
        Sp=1./(1.+std::exp(-S));
        tempv.push_back(Sp);
      }
      tempm.push_back(tempv);
    }


    // P=A;
    P.W_in(tempm);

    return P;
}

xy_W_bag xy_MLP::calc_targ(xy_W_bag &A, xy_W_bag &B)
{
    xy_W_bag P;

    std::vector<std::vector<double>> A_w;
    A_w=A.Get_W();

    std::vector<std::vector<double>> B_w;
    B_w=B.Get_W();

    std::vector<double> B_b;
    B_b=B.Get_b();

    unsigned int nx=A_w.size();
    if(nx<1){
      std::cout << "calc_prob input A matrix size=0." << std::endl;
      return P;
    }

    unsigned int ny=A_w.at(0).size();
    for(unsigned int i=0;i<nx;i++){
      if(ny!=A_w.at(i).size()){
        std::cout << "calc_prob: A matrix row size not same." << std::endl;
        return P;
      }
    }


    unsigned int nx2=B_w.size();
    if(nx2<1){
      std::cout << "calc_prob: input B matrix size=0." << std::endl;
      return P;
    }

    unsigned int ny2=B_b.size();
    for(unsigned int j=0;j<nx2;j++){
      if(ny2!=B_w.at(j).size()){
        std::cout << "calc_prob: B matrix row size not same." << std::endl;
        return P;
      }
    }


    if(ny!=nx2){
      std::cout << "calc_prob: input A B matrix dimension mismatch" << std::endl;
      return P;
    }

    std::vector<double> tempv;
    std::vector<std::vector<double>> tempm;
    tempm.clear();
    double d,vh,b,S,Sp;
    for(unsigned int i=0;i<nx;i++){
      tempv.clear();
      for(unsigned int k=0;k<ny2;k++){
        S=0.;
        for(unsigned int j=0;j<ny;j++){
          d=A_w.at(i).at(j);
          vh=B_w.at(j).at(k);
          S=S+d*vh;
        }
        b=B_b.at(k);
        S=S+b;
        Sp=std::exp(S);
        tempv.push_back(Sp);
      }
      tempm.push_back(tempv);
    }

    for(unsigned int i=0;i<nx;i++){
      S=0.;
      for(unsigned int k=0;k<ny2;k++){
        d=tempm.at(i).at(k);
        S=S+d;
      }
      if(S==0.)continue;

      for(unsigned int k=0;k<ny2;k++){
        d=tempm.at(i).at(k);
        d=d/S;
        tempm.at(i).at(k)=d;
      }
    }

    P.W_in(tempm);

    return P;
}

double xy_MLP::calc_err(xy_W_bag &A, xy_W_bag &B)
{
    double P=0.;

    std::vector<std::vector<double>> A_w;
    A_w=A.Get_W();

    std::vector<std::vector<double>> B_w;
    B_w=B.Get_W();

    unsigned int nx=A_w.size();
    if(nx<1){
      std::cout << "calc_err input A matrix size=0." << std::endl;
      return P;
    }

    unsigned int ny=A_w.at(0).size();
    for(unsigned int i=0;i<nx;i++){
      if(ny!=A_w.at(i).size()){
        std::cout << "calc_err: A matrix row size not same." << std::endl;
        return P;
      }
    }

    unsigned int nx2=B_w.size();
    if(nx2<1){
      std::cout << "calc_err input B matrix size=0." << std::endl;
      return P;
    }

    unsigned int ny2=B_w.at(0).size();
    for(unsigned int i=0;i<nx2;i++){
      if(ny2!=B_w.at(i).size()){
        std::cout << "calc_err: B matrix row size not same." << std::endl;
        return P;
      }
    }

    if(nx!=nx2 || ny!=ny2){
      std::cout << "calc_err: matrices A and B have different dimensions." << std::endl;
      return P;
    }


    double d,vh,S;
    S=0.;
    for(unsigned int i=0;i<nx;i++){
      for(unsigned int j=0;j<ny;j++){
        d=A_w.at(i).at(j);
        vh=B_w.at(i).at(j);
        if(vh>0.){
          vh=std::log(vh);
        }else{
          vh=0.;
        }
        S=S+d*vh;
      }
    }

    P=-S;

    return P;
}

std::vector<xy_W_bag> xy_MLP::calc_d_class(std::vector<xy_W_bag> &A)
{
    std::vector<xy_W_bag> P;
    P.clear();
    P.resize(2);

    unsigned int size=A.size();
    if(size!=3){
      std::cout << "calc_d_class: input bag # NOT 3" << std::endl;
      return P;
    }

    std::vector<std::vector<double>> t_out;
    t_out=A.at(0).Get_W();

    unsigned int nx_t_out=t_out.size();
    if(nx_t_out<1){
      std::cout << "calc_d_class:input bag 1 is null." << std::endl;
      return P;
    }

    unsigned int ny_t_out=t_out.at(0).size();
    for(unsigned int i=0;i<nx_t_out;i++){
      if(ny_t_out!=t_out.at(i).size()){
        std::cout << "calc_d_class: 1st matrix row size not same." << std::endl;
        return P;
      }
    }
    // std::cout << "calc_d_class t_out size= " << nx_t_out << " X " << ny_t_out << std::endl;


    std::vector<std::vector<double>> t_org;
    t_org=A.at(1).Get_W();

    unsigned int nx_t_org=t_org.size();
    if(nx_t_org<1){
      std::cout << "calc_d_class:input bag 2 is null." << std::endl;
      return P; 
    } 
    
    unsigned int ny_t_org=t_org.at(0).size();
    for(unsigned int i=0;i<nx_t_org;i++){
      if(ny_t_org!=t_org.at(i).size()){
        std::cout << "calc_d_class: 1st matrix row size not same." << std::endl;
        return P; 
      } 
    } 
    // std::cout << "calc_d_class t_org size= " << nx_t_org << " X " << ny_t_org << std::endl;

    if(nx_t_out!=nx_t_org || ny_t_out!=ny_t_org){
      std::cout << "calc_d_class: target matrices have different dimension." << std::endl;
      return P;
    }


    std::vector<std::vector<double>> wp;
    wp=A.at(2).Get_W();

    unsigned int nx_wp=wp.size();
    if(nx_wp<1){
      std::cout << "calc_d_class:input bag 3 is null." << std::endl;
      return P;
    }
    
    unsigned int ny_wp=wp.at(0).size();
    for(unsigned int i=0;i<nx_wp;i++){
      if(ny_wp!=wp.at(i).size()){
        std::cout << "calc_d_class: 1st matrix row size not same." << std::endl;
        return P;
      }
    }
    // std::cout << "calc_d_class wp size= " << nx_wp << " X " << ny_wp << std::endl;

    if(nx_wp!=nx_t_out){
      std::cout << "nx in target and wprob matrices different." << std::endl;
      return P;
    }


    std::vector<std::vector<double>> IX;
    IX.clear();
    IX=t_out;
    for(unsigned int i=0;i<nx_t_out;i++){
      for(unsigned int j=0;j<ny_t_out;j++){
        IX.at(i).at(j)=IX.at(i).at(j)-t_org.at(i).at(j);
      }
    }

    P.at(0).W_in(IX);

    std::vector<double> tempv;
    std::vector<std::vector<double>> tempm;
    tempm.clear();
    double temp_d,temp_b,temp_S;
    for(unsigned int k=0;k<ny_wp;k++){
      tempv.clear();
      for(unsigned int j=0;j<ny_t_out;j++){
        temp_S=0.;
        for(unsigned int i=0;i<nx_t_out;i++){
          temp_d=wp.at(i).at(k);
          temp_b=IX.at(i).at(j);
          temp_S=temp_S+temp_d*temp_b;
        }
        tempv.push_back(temp_S);
      }
      tempm.push_back(tempv);
    }

    P.at(1).W_in(tempm);

    tempv.clear();
    for(unsigned int j=0;j<ny_t_out;j++){
      temp_S=0.;
      for(unsigned int i=0;i<nx_t_out;i++){
        temp_b=IX.at(i).at(j);
        temp_S=temp_S+temp_b;
      }
      tempv.push_back(temp_S);
    }

    P.at(1).b_in(tempv);
    // P.at(1).Wb_comb();


    // P=A;

    return P;
}

std::vector<xy_W_bag> xy_MLP::calc_d_base(std::vector<xy_W_bag> &A)
{
    std::vector<xy_W_bag> P;
    P.clear();
    P.resize(2);

    unsigned int size=A.size();
    // std::cout << "Input bags: " << size << std::endl;
    if(size!=4){
      std::cout << "calc_d_base: input bag # NOT 4" << std::endl;
      return P;
    }

    std::vector<std::vector<double>> IX;
    IX=A.at(0).Get_W();

    unsigned int nx_IX=IX.size();
    if(nx_IX<1){
      std::cout << "calc_d_base:input bag 1 is null." << std::endl;
      return P;
    }

    unsigned int ny_IX=IX.at(0).size();
    for(unsigned int i=0;i<nx_IX;i++){
      if(ny_IX!=IX.at(i).size()){
        std::cout << "calc_d_base: 1st matrix row size not same." << std::endl;
        return P;
      }
    }
    // std::cout << "calc_d_base IX size= " << nx_IX << " X " << ny_IX << std::endl;


    std::vector<std::vector<double>> wcl;
    wcl=A.at(1).Get_W();

    unsigned int nx_wcl=wcl.size();
    if(nx_wcl<1){
      std::cout << "calc_d_base:input bag 2 is null." << std::endl;
      return P;
    }

    unsigned int ny_wcl=wcl.at(0).size();
    for(unsigned int i=0;i<nx_wcl;i++){
      if(ny_wcl!=wcl.at(i).size()){
        std::cout << "calc_d_base: 2nd matrix row size not same." << std::endl;
        return P;
      }
    }
    // std::cout << "calc_d_base wcl size= " << nx_wcl << " X " << ny_wcl << std::endl;


    std::vector<std::vector<double>> wpa;
    wpa=A.at(2).Get_W();

    unsigned int nx_wpa=wpa.size();
    if(nx_wpa<1){
      std::cout << "calc_d_base:input bag 3 is null." << std::endl;
      return P;
    }

    unsigned int ny_wpa=wpa.at(0).size();
    for(unsigned int i=0;i<nx_wpa;i++){
      if(ny_wpa!=wpa.at(i).size()){
        std::cout << "calc_d_base: 3rd matrix row size not same." << std::endl;
        return P;
      }
    }
    // std::cout << "calc_d_base w_prob a size= " << nx_wpa << " X " << ny_wpa << std::endl;


    std::vector<std::vector<double>> wpb;
    wpb=A.at(3).Get_W();

    unsigned int nx_wpb=wpb.size();
    if(nx_wpb<1){
      std::cout << "calc_d_base:input bag 4 is null." << std::endl;
      return P;
    }

    unsigned int ny_wpb=wpb.at(0).size();
    for(unsigned int i=0;i<nx_wpb;i++){
      if(ny_wpb!=wpb.at(i).size()){
        std::cout << "calc_d_base: 4th matrix row size not same." << std::endl;
        return P;
      }
    }
    // std::cout << "calc_d_base w_prob b size= " << nx_wpb << " X " << ny_wpb << std::endl;

    if(ny_IX!=ny_wcl || nx_wcl!=ny_wpa || nx_IX!=nx_wpa){
      std::cout << "calc_d_base: IX, Wcl, Wpa matrix product dimension mismatch." << std::endl;
      return P;
    }

    if(nx_IX!=nx_wpb){
      std::cout << "calc_d_base: IX, Wpb matrix product dimension mismatch." << std::endl;
      return P;
    }

    std::vector<double> tempv;
    std::vector<std::vector<double>> tempm;
    tempm.clear();
    double temp_d,temp_b,temp_S;
    for(unsigned int i=0;i<nx_IX;i++){
      tempv.clear();
      for(unsigned int k=0;k<nx_wcl;k++){
        temp_S=0.;
        for(unsigned int j=0;j<ny_IX;j++){
          temp_d=IX.at(i).at(j);
          temp_b=wcl.at(k).at(j);
          temp_S=temp_S+temp_d*temp_b;
        }
        temp_d=wpa.at(i).at(k);
        temp_b=1.-temp_d;
        temp_S=temp_S*temp_d*temp_b;
        tempv.push_back(temp_S);
      }
      tempm.push_back(tempv);
    }

    std::vector<std::vector<double>> IXN;
    IXN.clear();
    IXN=tempm;

    P.at(0).W_in(IXN);

    tempm.clear();
    for(unsigned int j=0;j<ny_wpb;j++){
      tempv.clear();
      for(unsigned int k=0;k<ny_wpa;k++){
        temp_S=0.;
        for(unsigned int i=0;i<nx_wpb;i++){
          temp_d=wpb.at(i).at(j);
          temp_b=IXN.at(i).at(k);
          temp_S=temp_S+temp_d*temp_b;
        }
        tempv.push_back(temp_S);
      }
      tempm.push_back(tempv);
    }

    P.at(1).W_in(tempm);

    tempv.clear();
    for(unsigned int k=0;k<ny_wpa;k++){
      temp_S=0.;
      for(unsigned int i=0;i<nx_wpb;i++){
        temp_b=IXN.at(i).at(k);
        temp_S=temp_S+temp_b;
      }
      tempv.push_back(temp_S);
    }

    P.at(1).b_in(tempv);

    return P;
}

void xy_MLP::init_Ps(const std::vector<double> &A)
{
    P_init=1;

    unsigned int size=A.size();

    if(size!=6){
      std::cout << "Init para: input dimension != 6, use default." << std::endl;

      P_INT=0.1; P_EXT=3.0; P_MAX=20.; P_RATIO=10.; P_SIG=0.1; P_RHO=0.05;
    }else{
      P_INT=A.at(0); P_EXT=A.at(1); P_MAX=A.at(2); P_RATIO=A.at(3); P_SIG=A.at(4); P_RHO=A.at(5);
    }

    std::cout << "P_INT, EXT, MAX, RATIO, SIG, RHO= " << P_INT << " " << P_EXT << " " << P_MAX << " " << P_RATIO << " " << P_SIG << " " << P_RHO << std::endl;

}

void xy_MLP::input_weights(std::vector<xy_W_bag> &A)
{
    unsigned int size=A.size();
    if(size<2){
      std::cout << "input_weights: bag # <2." << std::endl;
      return;
    }

    ws=A;

}

void xy_MLP::input_data(xy_W_bag &A)
{
    unsigned int size=ws.size();
    if(size<2){
      std::cout << "input weights first, so dimension defined." << std::endl;
      return;
    }

    wps.resize(size);

    wps.at(0)=A;

}

void xy_MLP::input_targ(xy_W_bag &A)
{
    targorg=A;
}

void xy_MLP::eval_probs_base()
{
    unsigned int size=ws.size();
    if(size<2){
      std::cout << "Weights not (properly) initialized." << std::endl;
      return;
    }

    unsigned int sizep=wps.size();
    if(size!=sizep){
      std::cout << "Data not been input properly." << std::endl;
      return;
    }

    sizep=size-1;
    for(unsigned int i=1;i<sizep;i++){
      wps.at(i)=calc_prob(wps.at(i-1),ws.at(i));
    }


}

void xy_MLP::eval_probs_targ()
{
    unsigned int size=ws.size();
    if(size<2){
      std::cout << "Weights not (properly) initialized." << std::endl;
      return;
    }

    unsigned int sizep=wps.size();
    if(size!=sizep){
      std::cout << "Data not been input properly." << std::endl;
      return;
    }

    sizep=size-1;

    unsigned int sizem;
    sizem=wps.at(sizep-1).Get_W().size();
    if(sizem<1){
      std::cout << "Penultimate prob not calculated." << std::endl;
      return;
    }

    wps.at(sizep)=calc_targ(wps.at(sizep-1),ws.at(sizep));

}

void xy_MLP::eval_d_class()
{
    unsigned int size=ws.size();
    if(size<2){
      std::cout << "Weights not (properly) initialized." << std::endl;
      return;
    }

    unsigned int sizep=wps.size();
    if(size!=sizep){
      std::cout << "Data not been input properly." << std::endl;
      return;
    }

    sizep=size-1;

    unsigned int sizet;
    sizet=targorg.Get_W().size();
    if(sizet<1){
      std::cout << "Original target not defined." << std::endl;
      return;
    }

    std::vector<xy_W_bag> bag_in;
    bag_in.push_back(wps.at(sizep));  //targetout
    bag_in.push_back(targorg);  //target original
    bag_in.push_back(wps.at(sizep-1));  //wprobs

    std::vector<xy_W_bag> bag_out;
    bag_out=calc_d_class(bag_in);

    IX_xy=bag_out.at(0);

    ws.at(sizep).d_W_in(bag_out.at(1).Get_W());
    ws.at(sizep).d_b_in(bag_out.at(1).Get_b());

}

void xy_MLP::eval_d_base()
{
    unsigned int size=ws.size();
    if(size<2){
      std::cout << "Weights not (properly) initialized." << std::endl;
      return;
    }

    unsigned int sizep=wps.size();
    if(size!=sizep){
      std::cout << "Data not been input properly." << std::endl;
      return;
    }

    sizep=size-1;

    unsigned int sizet;
    sizet=IX_xy.Get_W().size();
    if(sizet<1){
      std::cout << "IX from eval_d_class not defined." << std::endl;
      return;
    }

    std::vector<xy_W_bag> bag_in;
    std::vector<xy_W_bag> bag_out;
    unsigned int idx;
    for(unsigned int i=1;i<sizep;i++){
      idx=sizep-i;

      bag_in.clear();
      bag_in.push_back(IX_xy);  //IX_class
      bag_in.push_back(ws.at(idx+1));  //w
      bag_in.push_back(wps.at(idx));  //w3probs
      bag_in.push_back(wps.at(idx-1));  //w2probs

      bag_out.clear();
      bag_out=calc_d_base(bag_in);

      IX_xy=bag_out.at(0);

      ws.at(idx).d_W_in(bag_out.at(1).Get_W());
      ws.at(idx).d_b_in(bag_out.at(1).Get_b());
    }


}

void xy_MLP::eval_err()
{
    unsigned int size=ws.size();
    if(size<2){
      std::cout << "Weights not (properly) initialized." << std::endl;
      return;
    }

    unsigned int sizep=wps.size();
    if(size!=sizep){
      std::cout << "Data not been input properly." << std::endl;
      return;
    }

    sizep=size-1;

    unsigned int sizet;
    sizet=targorg.Get_W().size();
    if(sizet<1){
      std::cout << "Original target not defined." << std::endl;
      return;
    }

    unsigned int sizeto;
    sizeto=wps.at(sizep).Get_W().size();
    if(sizeto!=sizet){
      std::cout << "Calculated target not found." << std::endl;
      return;
    }

    error=calc_err(targorg,wps.at(sizep));
    // std::cout << "error= " << error << std::endl;

}

int xy_MLP::check_run()
{
    int P=0;

    if(std::isinf(error) || std::isnan(error)){
      return 1;
    }

    unsigned int size=ws.size();

    std::vector<std::vector<double>> dW;
    std::vector<double> db;

    unsigned int nx;
    unsigned int ny;

    for(unsigned int i=1;i<size;i++){
      dW.clear();
      db.clear();

      dW=ws.at(i).Get_d_W();
      db=ws.at(i).Get_d_b();

      nx=dW.size();
      ny=db.size();

      if(nx<1 || ny<1){
        return 1;
      }

      for(unsigned int j=0;j<nx;j++){
        if(ny!=dW.at(j).size()){
          return 1;
        }
        for(unsigned int k=0;k<ny;k++){
          if(std::isinf(dW.at(j).at(k)) || std::isnan(dW.at(j).at(k))){
            return 1;
          }
        }
      }

      for(unsigned int k=0;k<ny;k++){
        if(std::isinf(db.at(k)) || std::isnan(db.at(k))){
          return 1;
        }
      }

    }



    return P;
}

std::vector<xy_W_bag> xy_MLP::add_d_W(std::vector<xy_W_bag> &A, double &x)
{
    std::vector<xy_W_bag> P;
    P.clear();

    unsigned int size=A.size();
    if(size<2){
      std::cout << "add_d_W: input bag # 0" << std::endl;
      return P;
    }

    P.resize(size);

    std::vector<std::vector<double>> W;
    std::vector<std::vector<double>> dW;
    std::vector<double> b;
    std::vector<double> db;

    unsigned int nx;
    unsigned int ny;

    for(unsigned int i=1;i<size;i++){
      W.clear();
      dW.clear();
      b.clear();
      db.clear();

      W=A.at(i).Get_W();
      b=A.at(i).Get_b();

      nx=W.size();
      ny=b.size();

      dW=A.at(i).Get_d_W();
      db=A.at(i).Get_d_b();

      if(nx!=dW.size() || ny!=db.size()){
        std::cout << "add_d_W: W dW, b db do not have same size." << std::endl;
        return P;
      }
      // std::cout << i << " nx= " << nx << " ny= " << ny << std::endl;

      for(unsigned int j=0;j<nx;j++){
        if(ny!=W.at(j).size()){
          std::cout << "add_d_W: W ny is not constant." << std::endl;
          return P;
        }
        for(unsigned int k=0;k<ny;k++){
          W.at(j).at(k)=W.at(j).at(k)+x*dW.at(j).at(k);
        }
      }
      P.at(i).W_in(W);

      for(unsigned int k=0;k<ny;k++){
        b.at(k)=b.at(k)+x*db.at(k);
      }
      P.at(i).b_in(b);

    }




    return P;
}

std::vector<xy_W_bag> xy_MLP::dW_to_W(std::vector<xy_W_bag> &A)
{
    std::vector<xy_W_bag> P;
    P.clear();

    unsigned int size=A.size();

    P.resize(size);

    for(unsigned int i=1;i<size;i++){
      P.at(i).W_in(A.at(i).Get_d_W());
      P.at(i).b_in(A.at(i).Get_d_b());
    }



    return P;
}

std::vector<xy_W_bag> xy_MLP::W_to_dW(std::vector<xy_W_bag> &A, std::vector<xy_W_bag> &B)
{
    std::vector<xy_W_bag> P;
    P.clear();

    unsigned int size=A.size();
    if(size!=B.size()){
      std::cout << "W_to_dW: A B have different dimension." << std::endl;
      return P;
    }

    P=B;

    for(unsigned int i=1;i<size;i++){
      P.at(i).d_W_in(A.at(i).Get_W());
      P.at(i).d_b_in(A.at(i).Get_b());
    }


    return P;
}

std::vector<xy_W_bag> xy_MLP::mat_sum(double &fa, std::vector<xy_W_bag> &A, double &fb, std::vector<xy_W_bag> &B)
{
    std::vector<xy_W_bag> P;
    P.clear();

    unsigned int size=A.size();
    if(size!=B.size()){
      std::cout << "mat_sum: A B have different dimension." << std::endl;
      return P;
    }

    P.resize(size);

    std::vector<std::vector<double>> W;
    std::vector<std::vector<double>> dW;
    std::vector<double> b;
    std::vector<double> db;

    unsigned int nx;
    unsigned int ny;


    for(unsigned int i=1;i<size;i++){
      W.clear();
      dW.clear();
      b.clear();
      db.clear();

      W=A.at(i).Get_W();
      b=A.at(i).Get_b();

      nx=W.size();
      ny=b.size();

      dW=B.at(i).Get_W();
      db=B.at(i).Get_b();

      if(nx!=dW.size() || ny!=db.size()){
        std::cout << "mat_sum: A B matrices do not have same size." << std::endl;
        return P;
      }

      for(unsigned int j=0;j<nx;j++){
        if(ny!=W.at(j).size()){
          std::cout << "mat_sum: W ny is not constant." << std::endl;
          return P;
        }
        for(unsigned int k=0;k<ny;k++){
          W.at(j).at(k)=fa*W.at(j).at(k)+fb*dW.at(j).at(k);
        }
      }
      P.at(i).W_in(W);

      for(unsigned int k=0;k<ny;k++){
        b.at(k)=fa*b.at(k)+fb*db.at(k);
      }
      P.at(i).b_in(b);

    }




    return P;
}

double xy_MLP::s_prod(std::vector<xy_W_bag> &A, std::vector<xy_W_bag> &B)
{
    double P;
    P=0.;

    unsigned int size=A.size();
    if(size!=B.size()){
      std::cout << "s_prod: A B have different dimension." << std::endl;
      return P;
    }

    std::vector<std::vector<double>> W;
    std::vector<std::vector<double>> dW;
    std::vector<double> b;
    std::vector<double> db;

    unsigned int nx;
    unsigned int ny;

    for(unsigned int i=1;i<size;i++){
      W.clear();
      dW.clear();
      b.clear();
      db.clear();

      W=A.at(i).Get_W();
      b=A.at(i).Get_b();

      nx=W.size();
      ny=b.size();

      dW=B.at(i).Get_W();
      db=B.at(i).Get_b();

      if(nx!=dW.size() || ny!=db.size()){
        std::cout << "s_prod: A B matrices do not have same size." << std::endl;
        return P;
      }

      for(unsigned int j=0;j<nx;j++){
        if(ny!=W.at(j).size()){
          std::cout << "s_prod: W ny is not constant." << std::endl;
          P=0.;
          return P;
        }
        for(unsigned int k=0;k<ny;k++){
          P=P+W.at(j).at(k)*dW.at(j).at(k);
        }
      }

      for(unsigned int k=0;k<ny;k++){
        P=P+b.at(k)*db.at(k);
      }

    }




    return P;
}

void xy_MLP::xy_minimize(double &f_type, double &f_leng)
{
    std::cout << "minimize: Type= " << f_type << " , length= " << f_leng << std::endl;
    // std::cout << "P_INT, EXT, MAX, RATIO, SIG, RHO= " << P_INT << " " << P_EXT << " " << P_MAX << " " << P_RATIO << " " << P_SIG << " " << P_RHO << std::endl;

    std::vector<xy_W_bag> ws0;
    ws0=ws;  //backup original weights

    std::vector<xy_W_bag> wsorg;
    wsorg=ws;  //backup original weights

    std::vector<xy_W_bag> wps0;
    wps0=wps;  //backup original probs

    unsigned int sizep;
    sizep=ws.size()-1;
    // std::cout << "penultimate layer: " << sizep << std::endl;

    double f_i=0.;
    double ls_failed=0.;

    if(f_type>0.5)eval_probs_base();
    eval_probs_targ();

    eval_d_class();
    if(f_type>0.5)eval_d_base();

    eval_err();
    double f0=error;

    // std::cout << "xy1| f0= " << f0 << std::endl;

    std::vector<xy_W_bag> s;
    s.clear();
    if(f_type>0.5){
      s=ws;
    }else{
      s.resize(2);
      s.at(1)=ws.at(sizep);
    }
    // std::cout << "s size= " << s.size() << std::endl;

    std::vector<xy_W_bag> tempbag;
    tempbag.clear();
    tempbag=dW_to_W(s);
    // std::cout << "tempbag size= " << tempbag.size() << std::endl;

    // cbag.clear();
    // cbag=s;
    // cbag=tempbag;

    std::vector<xy_W_bag> df0;
    df0=tempbag;

    double temp_d=-1.;
    double temp_b=0.;

    s=mat_sum(temp_d,tempbag,temp_b,tempbag);  //-df0

    tempbag.clear();

    // cbag.clear();
    // cbag=s;

    double d0;

    d0=s_prod(s,s);
    d0=-d0;

    // std::cout << "d0= " << d0 << std::endl;

    double x3;
    x3=1./(1.-d0);

    // std::cout << "x3= " << x3 << std::endl;


    double F0;
    std::vector<xy_W_bag> dF0;

    std::vector<xy_W_bag> w_up;  //update ws

    double x2,f2,d2,f3,d3;
    std::vector<xy_W_bag> df3;
    df3.clear();

    double x1,d1,f1;
    double A,B;

    double x4,d4,f4;

    double done_ext=0.;
    double success=0.;
    double M=P_MAX;
    double done_int=0.;

    double stop_ls=0.;

    int cxy;  //check run error

    while(f_i<f_leng && stop_ls<0.5){   //number of line searches loop
      f_i=f_i+1;

      // std::cout << std::endl;
      std::cout << "minimize: f_i= " << f_i << " , ls_failed= " << ls_failed << std::endl;

      wsorg=ws;  //weights updated in ws

      ws0=ws;
      F0=f0;
      dF0=df0;

      M=P_MAX;

      // std::cout << "F0= " << F0 << std::endl;
      // cbag.clear();
      // cbag=dF0;

      done_ext=0.;
      while(done_ext<0.5){  //extrapolation loop
        x2=0.;f2=f0;d2=d0;f3=f0;
        df3=df0;

        success=0.;
        // std::cout << "xy_EXT: done_ext= " << done_ext << " , success= " << success << std::endl;

        while(success<0.5 && M>0.5){  //try error loop
          M=M-1.;
          // std::cout << "try error loop: success= " << success << " , M= " << M << std::endl;

          tempbag.clear();
          if(f_type>0.5){
            tempbag=ws;
          }else{
            tempbag.resize(2);
            tempbag.at(1)=ws.at(sizep);
          }

          w_up=W_to_dW(s,tempbag);
          tempbag.clear();

          tempbag=add_d_W(w_up,x3);
          w_up=tempbag;
          tempbag.clear();

          // cbag.clear();
          // cbag=w_up;

          if(f_type>0.5){
            ws=w_up;
          }else{
            ws.at(sizep)=w_up.at(1);
          }

          if(f_type>0.5)eval_probs_base();
          eval_probs_targ();

          eval_d_class();
          if(f_type>0.5)eval_d_base();

          eval_err();
          f3=error;

          // std::cout << "f3= " << f3 << std::endl;

          df3.clear();
          tempbag.clear();
          if(f_type>0.5){
            tempbag=ws;
          }else{
            tempbag.resize(2);
            tempbag.at(1)=ws.at(sizep);
          }
          df3=dW_to_W(tempbag);
          tempbag.clear();

          // cbag.clear();
          // cbag=df3;

          success=1.;

          cxy=check_run();
          // std::cout << "Check error: " << cxy << std::endl;

          if(cxy>0. || f3==0.){
            // std::cout << "check error>0. Reduce x3 size." << std::endl;
            x3=(x2+x3)/2.0;
            success=0.;
            ws=wsorg;
          }

        }  //try error loop

        std::cout << "Try err loop in EXT finished: M, x2, f2, d2, f3= " << M << " , " << x2 << " , " << f2 << " , " << d2 << " , " << f3 << std::endl;


        if(f3<F0){  //update or not
          F0=f3;
          dF0=df3;
          ws0=ws;
        }
        ws=wsorg;

        d3=s_prod(df3,s);
        // std::cout << "d3= " << d3 << std::endl;


        // std::cout << "done_ext= " << done_ext << std::endl;
        if(d3>P_SIG*d0 || f3>(f0+x3*P_RHO*d0) || M<1.){
          // std::cout << "Done with extrapolation." << std::endl;
          done_ext=1.;
          continue;
        }
        // std::cout << "done_ext= " << done_ext << std::endl;

        // std::cout << "Not done yet, test. done_ext= " << done_ext << std::endl;

        x1=x2;f1=f2;d1=d2;
        x2=x3;f2=f3;d2=d3;
        A = 6.*(f1-f2)+3.*(d2+d1)*(x2-x1);
        B = 3.*(f2-f1)-(2.*d1+d2)*(x2-x1);
        x3 = x1-d1*(x2-x1)*(x2-x1)/(B+std::sqrt(B*B-A*d1*(x2-x1)));

        // std::cout << "x1,x2,x3= " << x1 << " , " << x2 << " , " << x3 << std::endl;
        // std::cout << "A, B= " << A << " , " << B << std::endl;
        // std::cout << "d1, d2= " << d1 << " , " << d2 << std::endl;
        // std::cout << "f1,f2,f3= " << f1 << " , " << f2 << " , " << f3 << std::endl;

        if(std::isinf(x3) || std::isnan(x3) || x3<0.){
          x3 = x2*P_EXT;
        }else if(x3 > x2*P_EXT){
          x3 = x2*P_EXT;
        }else if(x3<(x2+P_INT*(x2-x1))){
          x3 = x2+P_INT*(x2-x1);
        }

        // std::cout << "x3= " << x3 << std::endl;

      }  //extrapolation loop

      std::cout << "EXT ends" << std::endl;

      // cbag.clear();
      /*if(f_type>0.5){
        cbag=ws;
      }else{
        cbag.resize(2);
        cbag.at(1)=ws.at(sizep);
      }
      */
      /*if(f_type>0.5){
        cbag=ws0;
      }else{
        cbag.resize(2);
        cbag.at(1)=ws0.at(sizep);
      }
      */

      // d3=-d3;  //just a test, do not do this
      // std::cout << "d3= " << d3 << " " << std::abs(d3) << " " << std::fabs(d3) << std::endl;
      // std::cout << "-P_SIG*d0= " << (-P_SIG*d0) << std::endl;
      // std::cout << "f3= " << f3 << std::endl;
      // std::cout << "f0+x3*P_RHO*d0= " << (f0+x3*P_RHO*d0) << std::endl;

      done_int=0.;
      if(std::fabs(d3)>(-P_SIG*d0) || f3>(f0+x3*P_RHO*d0)){
        done_int=0.;
      }else{
        done_int=1.;
      }
      if(M<1.)done_int=1.;


      while(done_int<0.5){  //interpolation loop
        // std::cout << "Interpolating. done_int= " << done_int << " | M= " << M << std::endl;

        x4 = x3; f4 = f3; d4 = d3;
        if(d3>0. || f3>(f0+x3*P_RHO*d0)){
          x4 = x3; f4 = f3; d4 = d3;
        }else{
          x2 = x3; f2 = f3; d2 = d3;
        }
        // std::cout << "x4, f4, d4= " << x4 << " " << f4 << " " << d4 << std::endl;

        if(f4>f0){
          x3 = x2-(0.5*d2*(x4-x2)*(x4-x2))/(f4-f2-d2*(x4-x2));
        }else{
          A = 6*(f2-f4)/(x4-x2)+3*(d4+d2);
          B = 3*(f4-f2)-(2*d2+d4)*(x4-x2);
          x3 = x2+(std::sqrt(B*B-A*d2*(x4-x2)*(x4-x2))-B)/A;
        }
        // std::cout << "x3= " << x3 << std::endl;

        if(std::isinf(x3) || std::isnan(x3)){
          x3 = (x2+x4)/2;
        }
        // std::cout << "INT: x3, x2, x4= " << x3 << " " << x2 << " " << x4 << std::endl;

        if(x3<(x4-P_INT*(x4-x2))){
          temp_d=x3;
        }else{
          temp_d=(x4-P_INT*(x4-x2));
        }

        if(temp_d>(x2+P_INT*(x4-x2))){
          x3=temp_d;
        }else{
          x3=(x2+P_INT*(x4-x2));
        }
        std::cout << "INT:x3= " << x3 << std::endl;


        tempbag.clear();
        if(f_type>0.5){
          tempbag=ws;
        }else{
          tempbag.resize(2);
          tempbag.at(1)=ws.at(sizep);
        }

        w_up=W_to_dW(s,tempbag);
        tempbag.clear();

        tempbag=add_d_W(w_up,x3);
        w_up=tempbag;
        tempbag.clear();

        // cbag.clear();
        // cbag=s;
        // cbag=tempbag;
        // cbag=w_up;

        if(f_type>0.5){
          ws=w_up;
        }else{
          ws.at(sizep)=w_up.at(1);
        }

        if(f_type>0.5)eval_probs_base();
        eval_probs_targ();

        eval_d_class();
        if(f_type>0.5)eval_d_base();

        eval_err();
        f3=error;

        // std::cout << "f3= " << f3 << std::endl;


        df3.clear();
        tempbag.clear();
        if(f_type>0.5){
          tempbag=ws;
        }else{
          tempbag.resize(2);
          tempbag.at(1)=ws.at(sizep);
        }
        df3=dW_to_W(tempbag);
        tempbag.clear();

        // cbag.clear();
        // cbag=df3;

        if(f3<F0){  //update or not
          F0=f3;
          dF0=df3;
          ws0=ws;
        }
        ws=wsorg;

        M=M-1.;

        d3=s_prod(df3,s);
        std::cout << "INT: d3= " << d3 << " , M= " << M << std::endl;

        if(std::fabs(d3)>(-P_SIG*d0) || f3>(f0+x3*P_RHO*d0)){
          done_int=0.;
          // std::cout << "Not done yet. M= " << M << std::endl;
        }else{
          done_int=1.;
        }

        if(M<1.){
          done_int=1.;
          // std::cout << "Run out of M. Done." << std::endl;
        }

        // done_int=1.;
      }  //interpolation loop

      std::cout << "INT done." << std::endl;

      // cbag.clear();
      /*if(f_type>0.5){
        cbag=ws;
      }else{
        cbag.resize(2);
        cbag.at(1)=ws.at(sizep);
      }
      */
      /*if(f_type>0.5){
        cbag=ws0;
      }else{
        cbag.resize(2);
        cbag.at(1)=ws0.at(sizep);
      }
      */

      if(std::fabs(d3)<(-P_SIG*d0) && f3<(f0+x3*P_RHO*d0)){  //ls not failed
      // if(1.<0.){  //test else output
        std::cout << "ls not failed. Value= " << f3 << std::endl;

        tempbag.clear();
        if(f_type>0.5){
          tempbag=ws;
        }else{
          tempbag.resize(2);
          tempbag.at(1)=ws.at(sizep);
        }

        w_up=W_to_dW(s,tempbag);
        tempbag.clear();

        tempbag=add_d_W(w_up,x3);
        w_up=tempbag;
        tempbag.clear();

        if(f_type>0.5){
          ws=w_up;
        }else{
          ws.at(sizep)=w_up.at(1);
        }

        f0=f3;

        temp_d=s_prod(df3,df3)-s_prod(df0,df3);
        temp_b=s_prod(df0,df0);
        // std::cout << "temp_d= " << temp_d << " , temp_b= " << temp_b << std::endl;

        temp_d=temp_d/temp_b;
        temp_b=-1.;
        tempbag.clear();
        tempbag=mat_sum(temp_d,s,temp_b,df3);
        s.clear();
        s=tempbag;

        cbag.clear();
        cbag=s;

        df0=df3;
        d3=d0;
        d0=s_prod(df0,s);

        // std::cout << "d3= " << d3 << " , d0= " << d0 << std::endl;

        if(d0>0.){
          temp_d=-1.;
          temp_b=0.;
          s=mat_sum(temp_d,df0,temp_b,df0);
          d0=s_prod(s,s);
          d0=-1.*d0;
        }

        // std::cout << "d3= " << d3 << " , d0= " << d0 << std::endl;

        // cbag.clear();
        // cbag=s;

        if(P_RATIO<(d3/d0)){
          x3=x3*P_RATIO;
        }else{
          x3=x3*d3/d0;
        }

        // std::cout << "x3= " << x3 << std::endl;

        ls_failed=0.;
      }else{  //ls failed
        std::cout << "ls failed." << std::endl;

        ws=ws0;
        f0=F0;
        df0=dF0;

        //break condition
        if(ls_failed>0.5){
          stop_ls=1.0;
          continue;
        }

        temp_d=-1.;
        temp_b=0.;
        s=mat_sum(temp_d,df0,temp_b,df0);
        d0=s_prod(s,s);
        d0=-1.*d0;
        x3=1./(1.-d0);

        // std::cout << "f0= " << f0 << " , d0= " << d0 << " , x3= " << x3 << std::endl;

        // cbag.clear();
        // cbag=s;

        ls_failed=1.;
      }  //ls faled or not

      cbag.clear();
      if(f_type>0.5){
        cbag=ws;
      }else{
        cbag.resize(2);
        cbag.at(1)=ws.at(sizep);
      }
      

      // std::cout << "line search end" << std::endl;
    }  //line search loop




}
