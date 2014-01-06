#include "schmidt.h"
#include <cassert>
#include <omp.h>

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
  //cout << "qshape = " << A.qshape() << endl;
  //cout << "dshape = " << A.dshape() << endl;
  //cout << "nblock = " << block.size() << endl;
  
  // now generate these blocks
  for (int i = 0; i < block.size(); ++i) {
    cout << block[i] << endl;
    // get information of each block
    Spin s = (qp[block[i][1]] == -1) ? (Spin::down) : (Spin::up);
    int nl = nsites - ql[block[i][0]] - lc;
    int nr = nsites - 1 - qr[block[i][2]] - rc;
    auto iter_l = lbasis -> iterator(nl);
    auto iter_r = rbasis -> iterator(nr);
    
    A.reserve(block[i]);
    // now fill in data    
    DArray<3> dense;
    dense.reference(*(A.find(block[i]) -> second));
    
    for (int j = 0; j < iter_l.size(); ++j) {
      for (int k = 0; k < iter_r.size(); ++k) {
        dense(j, 0, k) = overlap(iter_l.get_config(j), iter_r.get_config(k), s, nl, nr);
      }
    }
    /*
    // electron numbers, since we use the Schmidt basis of the "left behind" sites, right-hand sites, spin quantum number reverses

    # pragma omp parallel for default(shared) schedule(static)    
    */
  }
  return std::move(A);
}
