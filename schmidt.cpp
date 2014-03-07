#include "schmidt.h"
#include <cassert>
#include <algorithm>    // std::find, std::reverse
#include <math.h>       /* sqrt */
#include <iomanip>
#include <boost/mpi.hpp>

using std::endl;
namespace mpi=boost::mpi;

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

vector<bool> ActiveSpaceIterator::bits(uint address, int nocc, bool half) const {
  vector<bool> temp_bits;
  int _nocc, _nsites = nsites;
  if (nocc >= 0) {
    _nocc = nocc;
  } else {
    _nocc = nex;
  }
  if (half) {
    _nsites /= 2;
  }

  for (int i = _nsites-1; i >= 0; --i) {
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

void ActiveSpaceIterator::initialize(int _nsites, int _nex, const SchmidtBasis* _basis) {
  nsites = _nsites;
  nex = _nex;
  basis = _basis;
  initialized = true;
  vector<int> nfirst, nsecond; // two halves
  for (int i = 0; i <= nex; ++i) {
    if (i <= nsites/2 && nex-i <= nsites/2) {
      nfirst.push_back(i);
      nsecond.push_back(nex-i);
    }
  }

  double threshold = basis -> get_thr();
  for (int i = 0; i < nfirst.size(); ++i) {
    // maximum indices    
    uint max_first = choose(nsites/2, nfirst[i]);
    uint max_second = choose(nsites/2, nsecond[i]);
    assert(max_first != 0 && max_second != 0);
    map<vector<bool>, double> weight_first, weight_second;
    // first do one iteration to find the possible indices for each spin    
    for (uint j = 0; j < max_first; ++j) {
      vector<bool> bit_rep = bits(j, nfirst[i], true);
      double w = basis -> get_weight(bit_rep, 0);
      if (w > threshold) {
        weight_first.insert(std::pair<vector<bool>, double>(bit_rep, w));
      }
    }
    for (uint j = 0; j < max_second; ++j) {
      vector<bool> bit_rep = bits(j, nsecond[i], true);      
      double w = basis -> get_weight(bit_rep, nsites/2);
      if (w > threshold) {
        weight_second.insert(std::pair<vector<bool>, double>(bit_rep, w));
      }
    }
    // the second stage: find possible pairs
    for (auto it1 = weight_first.cbegin(); it1 != weight_first.cend(); ++it1) {
      for (auto it2 = weight_second.cbegin(); it2 != weight_second.cend(); ++it2) {
        if (it1->second * it2->second > threshold) {
          vector<bool> merge = it1->first;
          merge.insert(merge.end(), it2->first.begin(), it2->first.end());
          list.push_back(merge);
        }
      }
    }
  }
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
  assert(weight.size() % 2 == 0);
  for (int i = 0; i < weight.size(); ++i) {
    assert(fabs(weight[i] + weight[idx_a.size()-1-i]-1.) < 1e-12);
  }
  // making iterators and compute dimensions
  dimensions();
  broadcast_iterators();
}

void SchmidtBasis::dimensions() {
  for (int i = ncore(); i <= ncore() + nactive(); ++i) {
    if (i % 2 == 0) {
      quantums.push_back(nsites()-i);
    }
  }
  mpi::communicator world;
  dims.resize(quantums.size(), 0);
  for (int i = 0; i < quantums.size(); ++i) {
    int nex = q2nex(quantums[i]);
    if (i % world.size() == world.rank()) {
      auto it = iterator(nex);
      dims[i] = it.size();
    } else {
      iterator(nex, false);
    }
    iterator_map.insert(std::pair<int, int>(nex, i % world.size()));
  }
  
  for (int i = 0; i < dims.size(); ++i) {
    broadcast(world, dims[i], iterator_on_rank(q2nex(quantums[i])));
    if (dims[i] == 0) {
      dims.erase(dims.begin()+i);
      quantums.erase(quantums.begin()+i);
      i -= 1;
    }
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
  os << "Quantum Numbers\t" << basis.quantums << endl;
  os << "Dimensions\t" << basis.dims << endl;
  //os << basis.active << endl;
  return os;
}

double SchmidtBasis::get_weight(const vector<bool>& bits, int shift) const {
  double w = 1.;
  for (int i = 0; i < bits.size(); ++i) {
    if (bits[i]) {
      w *= weight[shift+i];
    } else {
      w *= 1.-weight[shift+i];
    }
  }
  return w;
}

ActiveSpaceIterator SchmidtBasis::iterator(int nex, bool init) { // occupation number of active space
  auto it = iterators.find(nex);
  if (it == iterators.end()) {
    ActiveSpaceIterator asi;
    if (init) {
      asi.initialize(nactive(), nex, this);
    }
    iterators.insert(std::pair<int, ActiveSpaceIterator>(nex, std::move(asi)));
  }
  return std::move(iterators.at(nex));
}

void SchmidtBasis::broadcast_iterators() {
  mpi::communicator world;
  for (int i = 0; i < quantums.size(); ++i) {
    broadcast(world, iterators[q2nex(quantums[i])], iterator_on_rank(q2nex(quantums[i])));
    iterators[q2nex(quantums[i])].set_basis(this);
  }
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
  // possible quantum numbers
  ql = lbasis->quantums;
  qr = rbasis->quantums;
  qp = {-1, 1};
  // dimensions
  dl = lbasis->dims;
  dr = rbasis->dims;
  dp = {1, 1};
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
  cc = lbasis->get_core().Rows(2, nsites).t() * rbasis->get_core().Rows(1, nsites-1)
      + lbasis->get_core().Rows(nsites+2, nsites*2).t() * rbasis->get_core().Rows(nsites, (nsites-1) * 2);
  ac = lbasis->get_active().Rows(2, nsites).t() * rbasis->get_core().Rows(1, nsites-1)
      + lbasis->get_active().Rows(nsites+2, nsites*2).t() * rbasis->get_core().Rows(nsites, (nsites-1) * 2);
  ca = lbasis->get_core().Rows(2, nsites).t() * rbasis->get_active().Rows(1, nsites-1)
      + lbasis->get_core().Rows(nsites+2, nsites*2).t() * rbasis->get_active().Rows(nsites, (nsites-1) * 2);
  aa = lbasis->get_active().Rows(2, nsites).t() * rbasis->get_active().Rows(1, nsites-1)
      + lbasis->get_active().Rows(nsites+2, nsites*2).t() * rbasis->get_active().Rows(nsites, (nsites-1) * 2);
  cs.ReSize(lc, 2);
  as.ReSize(la, 2);
  cs.Column(1) << lbasis->get_core().Row(nsites+1).t();
  cs.Column(2) << lbasis->get_core().Row(1).t();
  as.Column(1) << lbasis->get_active().Row(nsites+1).t();
  as.Column(2) << lbasis->get_active().Row(1).t();
}

double CoupledBasis::overlap(const vector<bool>& left, const vector<bool>& right, Spin s, int nl, int nr) const {
  // we assume the number of orbitals are the same, to enhance performance
  // if the numbers are different, the overlap is 0, but we will get error here
  int ns = 2-int(s)*2;
  Matrix mat(nl+lc, nl+lc);
  if (cc.Storage()) { // core-core part
    mat.SubMatrix(nl+1, nl+lc, nl+lc-rc+1 , nl+lc) = cc;
  }
  if (nl) { // active-core part
    int count = 0;
    for (int i = 0; i < la; ++i) {
      if (left[i]) {
        mat.SubMatrix(count+1, count+1, nl+lc-rc+1, nl+lc) = ac.Row(i+1);
        if (ns) {
          mat(count+1, 1) = as(i+1, 1);
          mat(count+1, 2) = as(i+1, 2);
        }
        ++count;        
      }
    }
  }
  if (nr) {
    int count = 0;
    for (int i = 0; i < ra; ++i) {  // core-active part
      if (right[i]) {
        mat.SubMatrix(nl+1, nl+lc, count+ns+1, count+ns+1) = ca.Column(i+1);
        ++count;        
      }
    }
  }
  if (nl && nr) {  // active-active part
    int count_l = 0;
    for (int i = 0; i < la; ++i) {
      if (left[i]) {
        int count_r = 0;
        for (int j = 0; j < ra; ++j) {
          if (right[j]) {
            mat(count_l+1, count_r+ns+1) = aa(i+1, j+1);
            ++count_r;
          }
        }
        ++count_l;
      }
    }
  }
  if (ns) {
    mat.SubMatrix(nl+1, nl+lc, 1, 2) = cs;
  }

  int sign = ((nr+rc)%2 == 0) ? 1:-1;
  double det = (nl+lc == 0) ? 1.: mat.Determinant();
  
  return sign * det;
}

CoupledBasis::~CoupledBasis() {
  lbasis = nullptr;
  rbasis = nullptr;
}

