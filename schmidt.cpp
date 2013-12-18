#include "schmidt.h"
#include <cassert>
#include <algorithm>    // std::find
#include <math.h>       /* sqrt */
#include <iomanip>
#include <algorithm>    // std::reverse

using std::endl;

uint choose(int iN, int iR){
    if (iR < 0 || iR > iN) {
        return 0;
    }
    if (iR > iN/2) {
      return choose(iN, iN-iR);
    }
    uint iComb = 1;
    int i = 0;
    while (i < iR) {
        ++i;
        iComb *= iN - i + 1;
        iComb /= i;
    }
    return iComb;
}

ActiveSpaceIterator::ActiveSpaceIterator(int _nsites, int _nocc): nsites(_nsites), nocc(_nocc), ptr(0) {
  max = choose(nsites, nocc);
}

uint ActiveSpaceIterator::addr(const vector<bool>& bits) const {
  int occ = 0;
  uint address = 0;
  for (int i = 0; i < nsites; ++i) {
    if (bits[i]) {
      address += choose(i, ++occ);
    }
  }
  return address;
}

vector<bool> ActiveSpaceIterator::bits(uint address) const {
  vector<bool> temp_bits;
  int _nocc = nocc;
  for (int i = nsites-1; i >= 0; --i) {
    if (address >= choose(i, _nocc)) {
      temp_bits.push_back(true);
      address -= choose(i, _nocc--);
    } else {
      temp_bits.push_back(false);
    }
  }
  std::reverse(temp_bits.begin(), temp_bits.end());
  return std::move(temp_bits);
}

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
  addr_ptr.resize(active.Ncols());
  std::fill(addr_ptr.begin(), addr_ptr.end(), 0);
}

std::ostream& operator <<(std::ostream& os, const SchmidtBasis& basis) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(10);
  os << "Core Orbitals (" << basis.ncore() << ")\n";
  //os << basis.core << endl;
  os << "Active Space (" << basis.nactive() << ") with weights:\n";
  for (int i = 0; i < basis.weight.size(); ++i) {
    os << basis.weight[i] << "  ";
  }
  os << endl;
  //os << basis.active << endl;
  return os;
}

CoupledBasis::CoupledBasis(const SchmidtBasis& _lbasis, const SchmidtBasis& _rbasis): lbasis(&_lbasis),  rbasis(&_rbasis) {
  assert(lbasis->nsites() == rbasis->nsites()+1);
  nsites = lbasis->nsites();
  lc = lbasis->ncore();
  la = lbasis->nactive();
  rc = rbasis->ncore();
  ra = rbasis->nactive();
  contract1p();
}

void CoupledBasis::contract1p() {
  cc = lbasis->get_core().Rows(2, nsites).t() * rbasis->get_core();
  ac = lbasis->get_active().Rows(2, nsites).t() * rbasis->get_core();
  ca = lbasis->get_core().Rows(2, nsites).t() * rbasis->get_active();
  aa = lbasis->get_active().Rows(2, nsites).t() * rbasis->get_active();
  cs = lbasis->get_core().Row(1).t();
  as = lbasis->get_active().Row(1).t();
}

CoupledBasis::~CoupledBasis() {
  lbasis = nullptr;
  rbasis = nullptr;
}

