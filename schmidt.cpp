#include "schmidt.h"
#include <cassert>
#include <algorithm>    // std::find
#include <math.h>       /* sqrt */
#include <iomanip>

using std::endl;

SchmidtBasis::SchmidtBasis(const Matrix& nat_orbs, const vector<double>& occ, double thr1p, double thrnp): thr(thrnp) {
  assert(nat_orbs.Ncols() == occ.size());
  vector<int> idx_c, idx_a;
  for (auto it = occ.begin(); it != occ.end(); ++it) {
    if (1. - *it < thr1p) {
      idx_c.push_back((int)(it-occ.begin()));
    } else if (*it > thr1p) {
      idx_a.push_back((int)(it-occ.begin()));
    }
  }
  core.ReSize(nat_orbs.Nrows(), idx_c.size());
  active.ReSize(nat_orbs.Nrows(), idx_a.size());
  for (int i = 0; i < idx_c.size(); ++i) {
    core.Column(i+1) << nat_orbs.Column(idx_c[i]+1);
  }
  for (int i = 0; i < idx_a.size(); ++i) {
    active.Column(i+1) << nat_orbs.Column(idx_a[i]+1);
    weight.push_back(occ.at(idx_a[i]));
  }
}

std::ostream& operator <<(std::ostream& os, const SchmidtBasis& basis) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(10);
  os << "Core Orbitals (" << basis.ncore() << ")\n";
  //os << basis.core << endl;
  os << "Active Space (" << basis.nactive() << ")\n";
  //os << basis.active << endl;
  for (int i = 0; i < basis.weight.size(); ++i) {
    cout << basis.weight[i] << "  ";
  }
  cout << endl;
  cout << basis.thr << endl;
  return os;
}
