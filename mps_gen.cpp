#include "schmidt.h"
#include <cassert>
#include <omp.h>
/*
QSDArray<3, Quantum> CoupledBasis::generate() {
  QSDArray<3, Quantum> A;
  // build qshape, dshape basically transform from the data we already have
  TVector<Qshapes<Quantum>, 3> qshape;
  TVector<Dshapes, 3> dshape;
  for (int i = 0; i < ql.size(); ++i) {
    qshape[0].push_back(Quantum(ql[i]));
    dshape[0].push_back(dl[i]);
  }
  for (int i = 0; i < qp.size(); ++i) {
    qshape[1].push_back(Quantum(qp[i]));
    dshape[1].push_back(dp[i]);
  }
  for (int i = 0; i < qr.size(); ++i) {
    qshape[2].push_back(-Quantum(qr[i]));
    dshape[2].push_back(dr[i]);
  }

  A.resize(Quantum::zero(), qshape, dshape, false);
  //cout << A.qshape() << endl;
  //cout << "nblock = " << block.size() << endl;
  # pragma omp parallel for shared(A, cout) schedule(dynamic, 1)
  for (int i = 0; i < block.size(); ++i) {
    //cout << block[i] << endl;

    // electron numbers, since we use the Schmidt basis of the "left behind" sites, right-hand sites, spin quantum number reverses
    int nla = (nsites-ql[block[i][0]])/2-lc;
    int nlb = (nsites+ql[block[i][0]])/2-lc;
    int nra = (nsites-qr[block[i][2]]-1)/2-rc;
    int nrb = (nsites+qr[block[i][2]]-1)/2-rc;
    Spin s;
    if (qp[block[i][1]] == -1) {
      s = Spin::down;
    } else {
      s = Spin::up;
    }
    auto iter_l = lbasis -> iterator(nla, nlb);
    auto iter_r = rbasis -> iterator(nra, nrb);
    int phys = qp[block[i][1]];
    # pragma omp critical
    A.reserve(block[i]);
    // now fill in data
    DArray<3> dense;
    dense.reference(*(A.find(block[i]) -> second));
    for (int j = 0; j < iter_l.size(); ++j) {
      for (int k = 0; k < iter_r.size(); ++k) {
        dense(j, 0, k) = overlap(iter_l.get_pair(j), iter_r.get_pair(k), s, nla, nlb, nra, nrb);
      }
    }
    # pragma omp critical    
    cout << nla << " " << nlb << " " << nra << " " << nrb << endl;    
  }
  return std::move(A);
}
*/
