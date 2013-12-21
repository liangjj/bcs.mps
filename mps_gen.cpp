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
    qshape[2].push_back(Quantum(qr[i]));
    dshape[2].push_back(dr[i]);
  }

  cout << qshape << endl;
  cout << dshape << endl;
  //A.resize(Quantum::zero());
  return std::move(A);
}
