#include "mps_op.h"
#include <fstream>

using std::ofstream;
using std::ifstream;

MPS<Quantum> allocate(int nsites, int M, int M_total, int M_min, bool dim) {
  // physical quantum numbers and dimensions
  Qshapes<Quantum> qp;
  physical(2, qp);
  cout << qp << endl;
  Dshapes dp(qp.size(), 1);

  // right quantum numbers and dimensions
  vector<Qshapes<Quantum>> qr(nsites);
  vector<Dshapes> dr(nsites);
  qr[0] = qp;
  dr[0] = dp;

  // for these sites
  for (int i = 1; i < nsites-1; i++) {
    qr[i] = qr[i-1] * qp;

    for (auto dim_last: dr[i-1]) {
      for (auto dim_phys: dp) {
        dr[i].push_back(dim_last*dim_phys);
      }
    }
    
    // sort before merge
    Quantum sortq;
    int sortd;
    for (int j = 0; j < qr[i].size(); ++j) {
      for (int k = j+1; k < qr[i].size(); ++k) {
        if (qr[i][k] < qr[i][j]) {
          sortq = qr[i][j];
          sortd = dr[i][j];
          qr[i][j] = qr[i][k];
          dr[i][j] = dr[i][k];
          qr[i][k] = sortq;
          dr[i][k] = sortd;
        }
      }
    }
    
    // now merge
    int j = 0;
    while (j < qr[i].size()) {
      int k = j+1;
      while (k < qr[i].size()) {
        if (qr[i][k] == qr[i][j]) {
          qr[i].erase(qr[i].begin()+k);
          dr[i][j] += dr[i][k];
          dr[i].erase(dr[i].begin()+k);
        } else {
          ++k;
        }
        if (dr[i][j] > M && M_total <= 0) {
          dr[i][j] = M;
        }
      }
      ++j;
    }

    if (M_total > 0) {
      int d_total = 0;
      for (int j = 0; j < dr[i].size(); j++) {
        d_total += dr[i][j];
      }
      if (d_total > M_total * dr[i].size()) {
        double sh =  (double) M_total * dr[i].size() / d_total;
        for (int j = 0; j < dr[i].size(); j++) {
          int nd = dr[i][j];
          dr[i][j] = (int)((double)dr[i][j] * sh);
          if (dr[i][j] < M_min) {
            dr[i][j] = (nd > M_min)?M_min:nd;
          }
        }
      }
    }
  }

  qr[nsites-1] = Qshapes<Quantum>(1, Quantum::zero());
  dr[nsites-1] = Dshapes(1, 1);

  if (dim) {
    Qshapes<Quantum> tmpq;
    Dshapes tmpd;    
    for (int i = nsites-2; i >=0; --i) {
      tmpq.clear();
      for(int j = 0;j < qr[i+1].size();++j) {
        for(int k = qp.size() - 1;k >= 0;--k) {
          tmpq.push_back(qr[i + 1][j] * (-qp[k]));
        }
      }
      tmpd.clear();
      for(int j = 0;j < dr[i+1].size();++j) {
        for(int k = dp.size() - 1;k >= 0;--k) {
          tmpd.push_back(dr[i+1][j]*dp[k]);
        }
      }
      // sort
      Quantum sortq;
      int sortd;
      for (int j = 0; j < tmpq.size(); ++j) {
        for (int k = j+1; k < tmpq.size(); ++k) {
          if (tmpq[k] < tmpq[j]) {
            sortq = tmpq[j];
            sortd = tmpd[j];            
            tmpq[j] = tmpq[k];
            tmpd[j] = tmpd[k];            
            tmpq[k] = sortq;
            tmpd[k] = sortd;
          }
        }
      }
      // merge
      int tmpd_total = 0;
      int j = 0;
      while (j < tmpq.size()) {
        int k = j+1;
        while (k < tmpq.size()) {
          if (tmpq[k] == tmpq[j]) {
            tmpq.erase(tmpq.begin()+k);
            tmpd[j] += tmpd[k];
            tmpd.erase(tmpd.begin()+k);
          } else {
            ++k;
          }
          if (tmpd[j] > M && M_total <= 0) {
            tmpd[j] = M;
          }
        }
        ++j;
      }

      if (M_total > 0) {
        int d_total = 0;
        for (int j = 0; j < tmpd.size(); j++) {
          d_total += tmpd[j];
        }
        if (d_total > M_total * tmpd.size()) {
          double sh =  (double) M_total * tmpd.size() / d_total;
          for (int j = 0; j < tmpd.size(); j++) {
            int nd = tmpd[j];
            tmpd[j] = (int)((double)tmpd[j] * sh);
            if (tmpd[j] < M_min) {
              tmpd[j] = (nd > M_min)?M_min:nd;
            }
          }
        }
      }

      // keep only blocks present in both qr[i] and tmpq
      for (int k = 0; k < qr[i].size(); ++k) {
        auto it = std::find(tmpq.begin(), tmpq.end(), qr[i][k]);
        if (it == tmpq.end()) {
          qr[i].erase(qr[i].begin() + k);
          dr[i].erase(dr[i].begin() + k);
          --k;
        } else {
          if (dr[i][k] > tmpd[it-tmpq.begin()]) {
            dr[i][k] = tmpd[it-tmpq.begin()];
          }
        }
      }
    }
  } else {
    Qshapes<Quantum> tmpq;
    for (int i = nsites-2; i >=0; --i) {
      tmpq.clear();
      tmpq = qr[i+1] * (-qp);
      // keep only blocks present in both qr[i] and tmpq
      for (int k = 0; k < qr[i].size(); ++k) {
        if (std::find(tmpq.begin(), tmpq.end(), qr[i][k]) == tmpq.end()) {
          qr[i].erase(qr[i].begin() + k);
          dr[i].erase(dr[i].begin() + k);
          --k;
        }
      }
    }
  }

  if (M_total > 0) {
    for (int i = nsites-2; i >=0; --i) {
      int d_total = 0;
      for (int j = 0; j < dr[i].size(); j++) {
        d_total += dr[i][j];
      }
      if (d_total > M_total * dr[i].size()) {
        double sh =  (double) M_total * dr[i].size() / d_total;
        for (int j = 0; j < dr[i].size(); j++) {
          int nd = dr[i][j];
          dr[i][j] = (int)((double)dr[i][j] * sh);
          if (dr[i][j] < M_min) {
            dr[i][j] = (nd > M_min)?M_min:nd;
          }
        }
      }
    }
  }

  // MPS declaration
  MPS<Quantum> A(nsites);

  // allocate the tensors
  TVector<Qshapes<Quantum>, 3> qshape;
  TVector<Dshapes, 3> dshape;
  // site 0
  Qshapes<Quantum> ql(1, Quantum::zero());
  Dshapes dl(ql.size(), 1);
  qshape = make_array(ql, qp, -qr[0]);
  dshape = make_array(dl, dp, dr[0]);
  A[0].resize(Quantum::zero(), qshape, dshape, false);
  // other sites
  # pragma omp parallel shared(nsites, qr, dr, A) private(qshape, dshape, ql, dl)
  {
  # pragma omp for schedule(guided, 1) 
  for (int i = 1; i < nsites; ++i) {
    ql = qr[i-1];
    dl = dr[i-1];
    qshape = make_array(ql, qp, -qr[i]);
    dshape = make_array(dl, dp, dr[i]);
    A[i].resize(Quantum::zero(), qshape, dshape, false);
  }
  }

  return A;
}

map<TVector<Quantum, 3>, TArray<double, 3>> allocate_blocks(const QSDArray<3,Quantum>& site) {
  map<TVector<Quantum, 3>, TArray<double, 3>> blocks;
  int idx_lim = 1;
  for (auto j:site.qshape()) {
    idx_lim *= j.size();
  }
  auto dshape = site.dshape();
  for (int j = 0; j < idx_lim; ++j) {
    auto index = site.index(j);
    auto qindex =  site.qindex(index);
    // nonzero
    if (qindex[0] * qindex[1] * qindex[2] == Quantum::zero()) {
      IVector<3> d;
      d[0] = dshape[0][index[0]];
      d[1] = dshape[1][index[1]];
      d[2] = dshape[2][index[2]];
      blocks[qindex] = TArray<double, 3>(d);
    }
  }
  return blocks;
}

void insert_blocks(QSDArray<3,Quantum>& site, map<TVector<Quantum, 3>, TArray<double, 3>> blocks) {
  int idx_lim = 1;
  for (auto j:site.qshape()) {
    idx_lim *= j.size();
  }
  for (int j = 0; j < idx_lim; j++) {
    auto index = site.index(j);
    auto qindex =  site.qindex(index);
    if (qindex[0] * qindex[1] * qindex[2] == Quantum::zero()) {
      site.insert(index, blocks[qindex]);
      // delete the blocks immediately
      blocks.erase(qindex);
    }
  }
}

const TVector<Quantum, 3> make_qindex(int lq, int rq) {
  TVector<Quantum, 3> vec;
  vec[0] = Quantum(-lq);
  vec[1] = Quantum(lq-rq);
  vec[2] = Quantum(rq);
  return vec;
}

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




