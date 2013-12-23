#include "mps_op.h"
#include <fstream>

using std::ofstream;
using std::ifstream;

void save_site(const MPS<Quantum>& mps, int site, const char *filename){
  char name[50];
  sprintf(name,"%s/%d.mps",filename,site);
  ofstream fout(name);
  boost::archive::binary_oarchive oar(fout);
  oar << mps[site];
}

void load_site(MPS<Quantum>& mps, int site ,const char *filename){
  char name[50];
  sprintf(name,"%s/%d.mps",filename,site);
  ifstream fin(name);
  boost::archive::binary_iarchive iar(fin);
  iar >> mps[site];
}

double norm_on_disk(MPS<Quantum> &mps, const char* filename) {
  QSDArray<2> E;
  int L = mps.size();
  load_site(mps, 0, filename);
  QSDcontract(1.0,mps[0],shape(0,1),mps[0].conjugate(),shape(0,1),0.0,E);
  mps[0].clear();
  // intermediate
  QSDArray<3> I;
  for(int i = 1; i < L; ++i){
    load_site(mps, i, filename);
    //construct intermediate, i.e. past X to E
    QSDcontract(1.0,E,shape(0),mps[i],shape(0),0.0,I);
    //clear structure of E
    E.clear();
    //construct E for site i by contracting I with Y
    QSDcontract(1.0,I,shape(0,1),mps[i].conjugate(),shape(0,1),0.0,E);
    mps[i].clear();
    I.clear();
    //bad style: if no blocks remain, return zero
    if(E.begin() == E.end()) {
      return 0.0;
    }
  }
  return sqrt((*(E.find(shape(0,0))->second))(0,0));
}


void normalize_on_disk(MPS<Quantum>& mps, const char* filename) {
  double norm = norm_on_disk(mps, filename);
  int L = mps.size();  
  double alpha = pow(1./norm, 1./(double)L);
  for (int i = 0; i < L; ++i) {
    load_site(mps, i, filename);
    QSDscal(alpha, mps[i]);
    save_site(mps, i, filename);
    mps[i].clear();
  }
}

void compress_on_disk(MPS<Quantum>& mps,const MPS_DIRECTION &dir,int D, const char *filename, bool store){
  int L = mps.size();//length of the chain
  double acc_norm = 1.;

  if(dir == MPS_DIRECTION::Left) {
    SDArray<1> S;//singular values
    QSDArray<2> V;//V^T
    QSDArray<3> U;//U --> unitary left normalized matrix
    load_site(mps, 0, filename);  // load first site
    for(int i = 0;i < L - 1;++i){
      cout << "Compress site " << i << endl;
      //redistribute the norm over the chain: for stability reasons
      // note this is not complete, the sites on the left are not affected
      // need final renormalization
      double nrm = sqrt(QSDdotc(mps[i],mps[i]));
      acc_norm *= pow(nrm, 1./(double)L);
      QSDscal(acc_norm / nrm, mps[i]);
      
      //then svd
      QSDgesvd(RightArrow,mps[i],S,U,V,D);
      //copy unitary to mpx
      QSDcopy(U,mps[i]);
      if (store) {
        save_site(mps, i, filename);
        mps[i].clear();
      }

      //paste S and V together
      SDdidm(S,V);
      // now read next site
      load_site(mps, i+1, filename);
      //and multiply with mpx on the next site
      U = mps[i + 1];
      //when compressing dimensions will change, so reset:
      mps[i + 1].clear();
      QSDcontract(1.0,V,shape(1),U,shape(0),0.0,mps[i + 1]);
    }
    cout << "Compress site " << L-1 << endl;    
    double nrm = sqrt(QSDdotc(mps[L-1],mps[L-1]));
    acc_norm *= pow(nrm, 1./(double)L);
    QSDscal(acc_norm/nrm,mps[L-1]);
    if (store) {
      save_site(mps, L-1, filename);
      mps[L-1].clear();
    }
  } else {//right
    SDArray<1> S;//singular values
    QSDArray<3> V;//V^T --> unitary right normalized matrix
    QSDArray<2> U;//U
    load_site(mps, L-1, filename);  // load first site
    for(int i = L - 1;i > 0;--i){
      cout << "Compress site " << i << endl;      
      //redistribute the norm over the chain: for stability reasons
      double nrm = sqrt(QSDdotc(mps[i],mps[i]));
      acc_norm *= pow(nrm, 1./(double)L);      
      QSDscal(acc_norm/nrm,mps[i]);
      //then SVD: 
      QSDgesvd(RightArrow,mps[i],S,U,V,D);
      //copy unitary to mpx
      QSDcopy(V,mps[i]);
      if (store) {
        save_site(mps, i, filename);
        mps[i].clear();
      }

      //paste U and S together
      SDdimd(U,S);
      // now read next site
      load_site(mps, i-1, filename);
      //and multiply with mpx on the next site
      V = mps[i - 1];
      //when compressing dimensions will change, so reset:
      mps[i - 1].clear();
      QSDcontract(1.0,V,shape(2),U,shape(0),0.0,mps[i - 1]);

    }
    cout << "Compress site " << 0 << endl;    
    double nrm = sqrt(QSDdotc(mps[0],mps[0]));
    acc_norm *= pow(nrm, 1./(double)L);    
    QSDscal(acc_norm/nrm,mps[0]);
    if (store) {
      save_site(mps, 0, filename);
      mps[0].clear();
    }
  }
  // now normalize all the sites
  if (store) {
    normalize_on_disk(mps, filename);
  } else {
    normalize(mps);
  }
}

tuple<SDArray<1>, Qshapes<Quantum>> Schmidt_on_disk(MPS<Quantum>& mps, int site, const char* filename) {
  vector<double> coef;
  int L = mps.size();
  if (site < 0) {
    site = L/2;
  }
  
  load_site(mps, 0, filename);  // load first site    
  for(int i = 0;i < site;++i){
    SDArray<1> S;//singular values
    QSDArray<2> V;//V^T
    QSDArray<3> U;//U --> unitary left normalized matrix
    
    //then svd
    QSDgesvd(RightArrow,mps[i],S,U,V,0);
    //copy unitary to mps
    QSDcopy(U,mps[i]);
    save_site(mps, i, filename);
    mps[i].clear();

    //paste S and V together
    SDdidm(S,V);
    // now read next site
    load_site(mps, i+1, filename);
    //and multiply with mpx on the next site
    U = mps[i + 1];
    //when compressing dimensions will change, so reset:
    mps[i + 1].clear();
    QSDcontract(1.0,V,shape(1),U,shape(0),0.0,mps[i + 1]);
  }

  // save matrix[site]
  save_site(mps, site, filename);
  mps[site].clear();

  load_site(mps, L-1, filename);  // load first site
  for(int i = L - 1;i > site;--i){
    SDArray<1> S;//singular values
    QSDArray<3> V;//V^T --> unitary right normalized matrix
    QSDArray<2> U;//U
    //then SVD: 
    QSDgesvd(RightArrow,mps[i],S,U,V,0);
    //copy unitary to mpx
    QSDcopy(V,mps[i]);
    save_site(mps, i, filename);
    mps[i].clear();

    //paste U and S together
    SDdimd(U,S);
    // now read next site
    load_site(mps, i-1, filename);
    //and multiply with mpx on the next site
    V = mps[i - 1];
    //when compressing dimensions will change, so reset:
    mps[i - 1].clear();
    QSDcontract(1.0,V,shape(2),U,shape(0),0.0,mps[i - 1]);
  }
  SDArray<1> S;//singular values
  QSDArray<3> V;//V^T
  QSDArray<2> U;//U --> unitary left normalized matrix
  QSDgesvd(RightArrow,mps[site],S,U,V,0);

  // save matrix[site]
  save_site(mps, site, filename);
  mps[site].clear();
  return std::make_tuple(S, U.qshape()[0]);
}




