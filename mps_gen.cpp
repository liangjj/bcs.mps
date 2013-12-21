#include "schmidt.h"
#include <cassert>

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
  
  for (int i = 0; i < block.size(); ++i) {
    A.reserve(block[i]);
    // now fill in data
    DArray<3> dense;
    dense.reference(*(A.find(block[i]) -> second));
    // electron numbers
    int nla = (nsites+ql[block[i][0]])/2-lc;
    int nlb = (nsites-ql[block[i][0]])/2-lc;
    int nra = (nsites+qr[block[i][2]]-1)/2-rc;
    int nrb = (nsites-qr[block[i][2]]-1)/2-rc;

    SymmetricMatrix sa(nla, nla), sb(nlb, nlb);

    auto iter_l = lbasis -> iterator(nla, nlb);
    auto iter_r = rbasis -> iterator(nra, nrb);
    int phys = qp[block[i][1]];
    assert(iter_l.size() == dense.shape()[0]);
    assert(iter_r.size() == dense.shape()[2]);

    for (int j = 0; j < iter_l.size(); ++j) {
      for (int k = 0; k < iter_r.size(); ++k) {
      }
    }
    cout << dense << endl;      
  }
  return std::move(A);
}
