#include "schmidt.h"
#include <cassert>
#include <algorithm>    // std::find, std::reverse
#include <math.h>       /* sqrt */
#include <iomanip>

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

ActiveSpaceIterator::ActiveSpaceIterator(int _nsites, int _noccA, int _noccB, const SchmidtBasis* _basis): nsites(_nsites), noccA(_noccA), noccB(_noccB), basis(_basis) {
  // maximum indices
  uint maxA = choose(nsites, noccA);
  uint maxB = choose(nsites, noccB);
  assert(maxA != 0 && maxB != 0);
  double threshold = basis -> get_thr();
  // first do one iteration to find the possible indices for each spin
  map<uint, double> weightA, weightB;
  for (uint i = 0; i < maxA; ++i) {
    double w = basis -> get_weight(bits(i, Spin::up));
    if (w > threshold) {
      weightA.insert(std::pair<uint, double>(i, w));
    }
  }
  for (uint i = 0; i < maxB; ++i) {
    double w = basis -> get_weight(bits(i, Spin::down));
    if (w > threshold) {
      weightB.insert(std::pair<uint, double>(i, w));
    }
  }

  // the second stage: find possible pairs
  for (auto itA = weightA.cbegin(); itA != weightA.cend(); ++itA) {
    for (auto itB = weightB.cbegin(); itB != weightB.cend(); ++itB) {
      if (itA->second * itB->second > threshold) {
        list.push_back(std::pair<uint, uint>(itA->first, itB->first));
      }
    }
  }
}

std::pair<vector<bool>, vector<bool>> ActiveSpaceIterator::get_pair(int i) const {
  auto bitsA = bits(list[i].first, Spin::up);
  auto bitsB = bits(list[i].second, Spin::down);
  return std::pair<vector<bool>, vector<bool>>(bitsA, bitsB);
}

SchmidtBasis::SchmidtBasis(const Matrix& nat_orbs, const vector<double>& occ, double thr1p, double thrnp): thr(thrnp) {
  assert(nat_orbs.Ncols() == occ.size());
  vector<int> idx_c, idx_a;
  for (int i = 0; i < occ.size(); ++i) {
    if (1-occ[i] < thr1p) {
      idx_c.push_back(i);
    } else if (occ[i] > thr1p) {
      idx_a.push_back(i);
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

CoupledBasis::CoupledBasis(SchmidtBasis& _lbasis, SchmidtBasis& _rbasis): lbasis(&_lbasis),  rbasis(&_rbasis) {
  assert(lbasis->nsites() == rbasis->nsites()+1);
  nsites = lbasis->nsites();
  lc = lbasis->ncore();
  la = lbasis->nactive();
  rc = rbasis->ncore();
  ra = rbasis->nactive();
  // contract one-particle orbitals
  contract1p();
  // calculate possible quantum numbers
  quantum_number();
  // calculate dimensions
  dimensions();
  // make possible block list
  for (int i = 0; i < ql.size(); ++i) {
    for (int j = 0; j < qp.size(); ++j) {
      for (int k = 0; k < qr.size(); ++k) {
        if (ql[i]+qp[j] == qr[k]) {
          IVector<3> temp = {i, j, k};
          block.push_back(temp);
        }
      }
    }
  }
}

void CoupledBasis::contract1p() {
  cc = lbasis->get_core().Rows(2, nsites).t() * rbasis->get_core();
  ac = lbasis->get_active().Rows(2, nsites).t() * rbasis->get_core();
  ca = lbasis->get_core().Rows(2, nsites).t() * rbasis->get_active();
  aa = lbasis->get_active().Rows(2, nsites).t() * rbasis->get_active();
  cs = lbasis->get_core().Row(1).t();
  as = lbasis->get_active().Row(1).t();
}

void CoupledBasis::quantum_number() {
  qp = {-1, 1};
  // nelec is the same as nsites (left)
  int nelec_l = nsites - 2*lc;
  int nelec_r = nsites - 2*rc - 1;
  
  for (int ka = 0; ka <= la; ++ka) {
    if (nelec_l-ka >= 0 && nelec_l-ka <= la) {
      ql.push_back(2*ka - nelec_l);
    }
  }
  for (int ka = 0; ka <= ra; ++ka) {
    if (nelec_r-ka >= 0 && nelec_r-ka <= ra) {
      qr.push_back(2*ka - nelec_r);
    }
  }
}

void CoupledBasis::dimensions() {
  dp = {1, 1};
  for (int q:ql) {
    int neleca = (nsites-q)/2 - lc;
    int nelecb = (nsites+q)/2 - lc;
    auto it = lbasis -> iterator(neleca, nelecb);
    dl.push_back(it.size());
  }
  for (int q:qr) {
    int neleca = (nsites-q-1)/2 - rc;
    int nelecb = (nsites+q-1)/2 - rc;
    auto it = rbasis -> iterator(neleca, nelecb);
    dr.push_back(it.size());
  }
  // now if any dimension is 0, delete that quantum number and dimension
  for (int i = 0; i < dl.size(); ++i) {
    if (dl[i] == 0) {
      dl.erase(dl.begin()+i);
      ql.erase(ql.begin()+i);
      i -= 1;
    }
  }
  for (int i = 0; i < dr.size(); ++i) {
    if (dr[i] == 0) {
      dr.erase(dr.begin()+i);
      qr.erase(qr.begin()+i);
      i -= 1;
    }
  }
}

double CoupledBasis::overlap(const std::pair<vector<bool>, vector<bool>> left, const std::pair<vector<bool>, vector<bool>> right, Spin s, int nla, int nlb, int nra, int nrb) const {
  Matrix mat_a(nla+lc, nla+lc), mat_b(nlb+lc, nlb+lc);
  int sa = 1-int(s);
  int sb = int(s);
  if (cc.Storage()) { // core-core part
    mat_a.SubMatrix(nla+1, nla+lc, nla+lc-rc+1, nla+lc) = cc;
    mat_b.SubMatrix(nlb+1, nlb+lc, nlb+lc-rc+1, nlb+lc) = cc;
  }
  if (nla) {   // active-core part (alpha)
    int count = 0;
    for (int i = 0; i < la; ++i) {
      if (left.first[i]) {
        mat_a.SubMatrix(count+1, count+1, nla+lc-rc+1, nla+lc) = ac.Row(i+1);
        if (sa) {
          mat_a(count+1, 1) = as(i+1, 1);
        }
        ++count;
      }
    }   
  }
  if (nra) {   // core-active part (alpha)
    int count = 0;
    for (int i = 0; i < ra; ++i) {
      if (right.first[i]) {
        mat_a.SubMatrix(nla+1, nla+lc, count+sa+1, count+sa+1) = ca.Column(i+1);
        ++count;
      }
    }
  }
  if (nla && nra) {   // active-active (alpha)
    int count_l = 0;
    for (int i = 0; i < la; ++i) {
      if (left.first[i]) {
        int count_r = 0;
        for (int j = 0; j < ra; ++j) {
          if (right.first[j]) {
            mat_a(count_l+1, count_r+sa+1) = aa(i+1, j+1);
            ++count_r;
          }
        }
        ++count_l;        
      }
    }
  }
  if (nlb) {  // active-core (beta)
    int count = 0;
    for (int i = 0; i < la; ++i) {
      if (left.second[i]) {
        mat_b.SubMatrix(count+1, count+1, nlb+lc-rc+1, nlb+lc) = ac.Row(i+1);
        if (sb) {
          mat_b(count+1, 1) = as(i+1, 1);
        }
        ++count;
      }
    }
  }
  if (nrb) { // core-active (beta)
    int count = 0;
    for (int i = 0; i < ra; ++i) {
      if (right.second[i]) {
        mat_b.SubMatrix(nlb+1, nlb+lc, count+sb+1, count+sb+1) = ca.Column(i+1);
        ++count;
      }
    }
  }
  if (nlb && nrb) {   // active-active (beta)
    int count_l = 0;
    for (int i = 0; i < la; ++i) {
      if (left.second[i]) {
        int count_r = 0;
        for (int j = 0; j < ra; ++j) {
          if (right.second[j]) {
            mat_b(count_l+1, count_r+sb+1) = aa(i+1, j+1);
            ++count_r;
          }
        }
        ++count_l;        
      }
    }
  }
  if (sa) {  // onsite (alpha)
    mat_a.SubMatrix(nla+1, nla+lc, 1, 1) = cs;
  } else {  // onsite (beta)
    mat_b.SubMatrix(nlb+1, nlb+lc, 1, 1) = cs;
  }
  double sign = 1.;
  if (sb && (nla+lc)%2 == 1) {
    sign = -1.;
  }
  double detA = (nla+lc == 0) ? 1.: mat_a.Determinant();
  double detB = (nlb+lc == 0) ? 1.: mat_b.Determinant();

  return sign * detA * detB;
}

CoupledBasis::~CoupledBasis() {
  lbasis = nullptr;
  rbasis = nullptr;
}

