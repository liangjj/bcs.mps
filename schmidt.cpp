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

ActiveSpaceIterator::ActiveSpaceIterator(int _nsites, int _noccA, int _noccB, const SchmidtBasis* const _basis): nsites(_nsites), noccA(_noccA), noccB(_noccB), basis(_basis) {
  maxA = choose(nsites, noccA);
  maxB = choose(nsites, noccB);
  assert(maxA != 0 && maxB != 0);
  double threshold = basis -> get_thr();
  // now do one iteration to store the possible indices for each spin
  for (int i = 0; i < maxA; ++i) {
    double w = basis -> get_weight(bits(i, Spin::up));
    if (w > threshold) {
      weightA.insert(std::pair<uint, double>(i, w));
    }
  }
  ptrA = weightA.cbegin();
  for (int i = 0; i < maxB; ++i) {
    double w = basis -> get_weight(bits(i, Spin::down));
    if (w > threshold) {
      weightB.insert(std::pair<uint, double>(i, w));
    }
  }
  ptrB = weightB.cbegin();
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

vector<bool> ActiveSpaceIterator::bits(uint address, Spin s) const {
  vector<bool> temp_bits;
  int _nocc;
  if (s == Spin::up) {
    _nocc = noccA;
  } else {
    _nocc = noccB;
  }
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

std::pair<vector<bool>, vector<bool>> ActiveSpaceIterator::fetch() const {
  auto bitsA = bits(ptrA->first, Spin::up);
  auto bitsB = bits(ptrB->first, Spin::down);
  return std::pair<vector<bool>, vector<bool>>(bitsA, bitsB);
}

void ActiveSpaceIterator::find() {
  double threshold = basis -> get_thr();
  while (!end()) {
    if (ptrA->second * ptrB->second > threshold) {
      break;
    }
    next();
  }
}

void ActiveSpaceIterator::next() {
  ++ptrB;
  if (ptrB == weightB.cend()) {
    ptrB = weightB.cbegin();
    ++ptrA;
  }
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

double SchmidtBasis::get_weight(const vector<bool>& bits) const {
  double w = 1.;
  for (int i = 0; i < bits.size(); ++i) {
    if (bits[i]) {
      w *= weight[i];
    } else {
      w *= 1.-weight[i];
    }
  }
  return w;
}

ActiveSpaceIterator SchmidtBasis::iterator(int noccA, int noccB) { // noccs of active space
  vector<int> occs = {noccA, noccB};
  auto it = iterators.find(occs);
  if (it == iterators.end()) {
    ActiveSpaceIterator asi(nactive(), noccA, noccB, this);
    iterators.insert(std::pair<vector<int>, ActiveSpaceIterator>(occs, std::move(asi)));
  }
  return std::move(iterators.at(occs));  
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

