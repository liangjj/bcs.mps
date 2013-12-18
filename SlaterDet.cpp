#include "SlaterDet.h"
#include <cassert>
#include <algorithm>    // std::find
#include <math.h>       /* sqrt */
#include <iomanip>

SlaterDet::SlaterDet(Matrix*& orb_up, Matrix*& orb_down, 
    vector<int>& isites, double ifactor): factor(ifactor), sites(isites) {
  assert((*orb_up).Nrows() == sites.size());
  assert((*orb_down).Nrows() == sites.size());
  orb[0] = orb_up;
  orb[1] = orb_down;
}

int SlaterDet::nelec(Spin s) const {
  return orb[int(s)]->Ncols();
}

double SlaterDet::get_factor() const {
  return factor;
}

int SlaterDet::nspin() const {
  return nelec(Spin::up) - nelec(Spin::down);
}

string SlaterDet::list_sites() const {
  string str = "";
  for (auto site: sites) {
    str +=  std::to_string(static_cast<long long>(site)) + "  ";
  }
  return str;
}

std::ostream& operator << (std::ostream& os, const SlaterDet& det) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(10);
  os << "Factor=" << det.factor << "\n";
  os << "Sites:\n";
  os << det.list_sites();
  os << "\n";
  os << "Number of particles: alpha " << det.nelec(Spin::up) << " beta " << det.nelec(Spin::down) << "\n";
  os << "----------\n" << *det.orb[int(Spin::up)]<< "\n";
  os <<  "----------\n" << *det.orb[int(Spin::down)]<< "\n";
  return os;
}

double operator*(const SlaterDet& det1, const SlaterDet& det2) {
  if (det1.nelec(Spin::up) != det2.nelec(Spin::up) || det1.nelec(Spin::down) != det2.nelec(Spin::down)) {
    return 0.;
  }
  double result = 1.;
  Spin ud[2] = {Spin::up, Spin::down};
  for (auto s:ud) {
    if (det1.nelec(s) != 0) {
      Matrix temp = det1.orb[int(s)]->t() * (*(det2.orb[int(s)]));
      result *= temp.Determinant();
    }
  }
  return result;
}


tuple<int, double, ColumnVector*, Matrix*> SlaterDet::make_emb(const Matrix& sorb, int site) {
  // projected overlap matrix
  assert(site == 1);
  SymmetricMatrix S;
  S << sorb.Row(site).t() * sorb.Row(site);
  DiagonalMatrix ew;
  Matrix ev;
  Jacobi(S, ew, ev);
  Matrix sorbR(sorb * ev);
  //cout << sorbR << endl;
  for (int i = 1; i < sorbR.Ncols(); ++i) {
    sorbR.Column(i) /= sorbR.Column(i).NormFrobenius();
  }
  double bath_factor = sorbR(site, sorbR.Ncols());
  int sign = 1-2*(ev.Determinant()<0); 
  if (sorb.Ncols() % 2 == 0) { // move the last column (bath) to the first
    sign *= -1;
  }
  //double norm_bath = sorbR.Column(sorbR.Ncols()).NormFrobenius();
  double norm_bath = sqrt(1.-bath_factor * bath_factor);
  //cout << fabs(bath_factor * bath_factor + norm_bath * norm_bath - 1.) << endl;
  Matrix *env = new Matrix(sorbR.Nrows()-1, sorb.Ncols()-1);
  env -> Rows(1, site-1) << sorbR.SubMatrix(1, site-1, 1, sorbR.Ncols()-1);
  env -> Rows(site, env -> Nrows()) << sorbR.SubMatrix(site+1, sorbR.Nrows(), 1, sorbR.Ncols()-1);
  // Matrix* env = new Matrix(sorbR.Columns(1, sorbR.Ncols() - 1));
  ColumnVector* bath = nullptr;
  if (norm_bath > 3e-6) {
    bath = new ColumnVector(sorbR.Nrows()-1);
    bath -> Rows(1, site-1) << sorbR.SubMatrix(1, site-1, sorbR.Ncols(), sorbR.Ncols());
    bath -> Rows(site, bath -> Nrows()) << sorbR.SubMatrix(site+1, sorbR.Nrows(),  sorbR.Ncols(), sorbR.Ncols());
    //if (fabs(norm_bath-bath->NormFrobenius())/(norm_bath+bath->NormFrobenius()) > 1e-6) {
    //  cout << norm_bath << "\t" << bath->NormFrobenius() << endl;
    //}
    *bath /= norm_bath;
  }
  return std::forward_as_tuple(sign, bath_factor, bath, env); // sign, factor, bath, env
}

vector<SlaterDet*> SlaterDet::Schmidt(int site, bool detail) {
  int site_idx = std::find(sites.begin(), sites.end(), site) - sites.begin() + 1;
  //for (int i = 0; i < sites.size(); ++i) {
  //  cout << sites[i] << "\t";
  //}
  //cout << "\n";
  assert(site_idx <= sites.size());
  vector<SlaterDet*> vec(2, nullptr);
  if (nelec(Spin::up) == 0 && nelec(Spin::down) == 0) { // no electron in both spin
  } else if (nelec(Spin::up) == 0) { // no electron in spin up
    auto emb = make_emb(*orb[int(Spin::down)], site_idx);
    Matrix* empty = new Matrix(sites.size()-1, 0);
    vector<int> rest_sites = sites;
    rest_sites.erase(std::find(rest_sites.begin(), rest_sites.end(), site));
    vec[1] = new SlaterDet(empty, get<3>(emb), rest_sites, factor * get<1>(emb) * get<0>(emb));
    delete get<2>(emb);
  } else if (nelec(Spin::down) == 0) { // no electron in spin down
    auto emb = make_emb(*orb[int(Spin::up)], site_idx);
    Matrix* empty = new Matrix(sites.size()-1, 0);
    vector<int> rest_sites = sites;
    rest_sites.erase(std::find(rest_sites.begin(), rest_sites.end(), site));
    vec[0] = new SlaterDet(get<3>(emb), empty, rest_sites, factor * get<1>(emb) * get<0>(emb));
    delete get<2>(emb);
  } else { // both spin have electrons
    auto emb_u = make_emb(*orb[int(Spin::up)], site_idx);
    auto emb_d = make_emb(*orb[int(Spin::down)], site_idx);
    double com_factor = factor * get<0>(emb_u) * get<0>(emb_d);
    if (detail) {
      cout << com_factor << endl;
      cout << factor << endl;
      cout << get<1>(emb_u) << endl;
      cout << get<1>(emb_d) << endl;
    }
    int sign = 1;
    if (nelec(Spin::up)%2 == 1) {
      sign = -1;
    }
    vector<int> rest_sites = sites;
    rest_sites.erase(std::find(rest_sites.begin(), rest_sites.end(), site));

    if ((get<2>(emb_d)) != nullptr) { // no bath in spin down
      Matrix* orb_d = new Matrix(rest_sites.size(), nelec(Spin::down));
      orb_d -> Column(1) << *get<2>(emb_d);
      orb_d -> Columns(2, nelec(Spin::down)) << *get<3>(emb_d);
      if (sqrt(1.-get<1>(emb_d) * get<1>(emb_d)) == sqrt(1.-get<1>(emb_d) * get<1>(emb_d))) {
        vec[0] = new SlaterDet(get<3>(emb_u), orb_d, rest_sites, com_factor * get<1>(emb_u) * sqrt(1.-get<1>(emb_d) * get<1>(emb_d)));
      }
      delete get<2>(emb_d);
    }
    if ((get<2>(emb_u)) != nullptr) { // no bath in spin up
      Matrix* orb_u = new Matrix(rest_sites.size(), nelec(Spin::up));
      orb_u -> Column(1) << *get<2>(emb_u);
      orb_u -> Columns(2, nelec(Spin::up)) << *get<3>(emb_u);
      if (sqrt(1.- get<1>(emb_u) * get<1>(emb_u)) == sqrt(1.- get<1>(emb_u) * get<1>(emb_u))) {
        vec[1] = new SlaterDet(orb_u, get<3>(emb_d), rest_sites, com_factor * sign * get<1>(emb_d) * sqrt(1.- get<1>(emb_u) * get<1>(emb_u)));
      }
      delete get<2>(emb_u);
    }
  }
  return vec;
}

SlaterDet::~SlaterDet() {
  for (int i = 0; i < 2; i++) {
    if (orb[i] != nullptr) {
      delete orb[i];
      orb[i] = nullptr;
    }
  }
}


void WfnContainer::push_back(SlaterDet*& det, int origin_spin, int origin_idx) {
  int spin = det->nspin();
  if (data.find(spin) == data.end()) {
    data[spin] = new vector<SlaterDet*>;
    block[spin] = new vector<int>;
    l_idx[spin] = new vector<int>;    
    factor[spin] = new vector<double>;
    keys.push_back(spin);
  }
  data[spin] -> push_back(det);
  block[spin] -> push_back(origin_spin);
  l_idx[spin] -> push_back(origin_idx);  
  sum += det -> get_factor() * det -> get_factor();
  factor[spin] -> push_back(fabs(det -> get_factor()));
  det = nullptr;
}

const vector<SlaterDet*>& WfnContainer::operator[](int spin) {
  return *(data[spin]);
}

const vector<int>& WfnContainer::qp() const {
  return keys;
}

int WfnContainer::size(int spin) const {
  return data.find(spin) -> second -> size();
}

std::ostream& operator <<(std::ostream& os, const WfnContainer& wfns) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(10);
  for (int spin: wfns.qp()) {
    os << "spin=" << spin << " : no. of wavefunctions " << wfns.size(spin) << "\n";
    auto map_it = wfns.data.find(spin);
    //for (auto it = map_it -> second -> begin(); it != map_it -> second -> end(); ++it) {
    //  cout << **it << "\n";
    //}
  }
  return os;
}

void WfnContainer::scale(SlaterDet& det, double s) {
  det.factor /= s;
}

void WfnContainer::renormalize() {
  // first discard states with two low norms
  for (int spin: qp()) {
    // do this only when there are more than M states associated with the quantum number
    if (size(spin) > M) {
      // find threshold of factor
      std::nth_element(factor[spin] -> begin(), factor[spin] -> end() - M, factor[spin] -> end());
      double thr = 1e-10;
      if (thr < *(factor[spin]->end() - M)) {
        thr = *(factor[spin]->end() - M);
      }
      auto& tmpd = *data[spin];
      auto& tmpb = *block[spin];
      auto& tmpl = *l_idx[spin];

      for (int i = 0; i < tmpd.size(); ++i) {
        if (fabs(tmpd[i] -> get_factor()) < thr+ 1e-10) {
          sum -= tmpd[i] -> get_factor() * tmpd[i] -> get_factor();
          delete tmpd[i];
          tmpd.erase(tmpd.begin()+i);
          tmpb.erase(tmpb.begin()+i);
          tmpl.erase(tmpl.begin()+i);
          --i;
        }
      }
    }
  }
  // now normalize all states
  double norm = sqrt(sum);
  //cout << sum << endl;
  sum = 0.;
  for (int spin: qp()) {
    for (int i = 0; i < data[spin] -> size(); ++i) {
      scale(*(data[spin] -> at(i)), norm);
      factor[spin] -> at(i) = fabs(data[spin] -> at(i) -> get_factor());
      sum += factor[spin] -> at(i) * factor[spin] -> at(i);
    }
  }
  //cout << sum << endl;
}

void WfnContainer::renormalize(map<int, int>& dims) {
  // first discard states with two low norms
  for (int spin: qp()) {
    // do this only when there are more than M states associated with the quantum number
    int sdim = dims[spin];
    if (size(spin) > sdim) {
      // find threshold of factor
      std::nth_element(factor[spin] -> begin(), factor[spin] -> end() - sdim, factor[spin] -> end());
      double thr = 1e-10;
      if (thr < *(factor[spin] -> end()-sdim)) {
        thr = *(factor[spin] -> end()-sdim);
      }
      auto& tmpd = *data[spin];
      auto& tmpb = *block[spin];
      auto& tmpl = *l_idx[spin];

      for (int i = 0; i < tmpd.size(); ++i) {
        if (fabs(tmpd[i] -> get_factor()) < thr+ 1e-10) {
          sum -= tmpd[i] -> get_factor() * tmpd[i] -> get_factor();
          delete tmpd[i];
          tmpd.erase(tmpd.begin()+i);
          tmpb.erase(tmpb.begin()+i);
          tmpl.erase(tmpl.begin()+i);
          --i;
        }
      }
    }
  }
  // now normalize all states
  double norm = sqrt(sum);
  sum = 0;
  for (int spin: qp()) {
    for (int i = 0; i < data[spin] -> size(); ++i) {
      scale(*(data[spin] -> at(i)), norm);
      factor[spin] -> at(i) = fabs(data[spin] -> at(i) -> get_factor());
      sum += factor[spin] -> at(i) * factor[spin] -> at(i);
    }
  }
}

void WfnContainer::clear() {
  for (auto it = data.begin(); it != data.end(); ++it) {
    for (auto it_vec = it -> second -> begin(); it_vec != it -> second -> end(); ++it_vec) {
      if (*it_vec != nullptr) {
        delete (*it_vec);
        *it_vec = nullptr;
      }
    }
    if (it -> second != nullptr) {
      delete it -> second;
      it -> second = nullptr;
    }
  }
  for (auto it = block.begin(); it != block.end(); ++it) {
    delete it -> second;
    it -> second = nullptr;
  }
  for (auto it = factor.begin(); it != factor.end(); ++it) {
    delete it -> second;
    it -> second = nullptr;
  }
  sum = 0;
  data.clear();
  keys.clear();
  block.clear();
  factor.clear();
}

const vector<int>& WfnContainer::lq(int spin) const {
  return *(block.find(spin) -> second);
}

const vector<int>& WfnContainer::ld(int spin) const {
  return *(l_idx.find(spin) -> second);
}

WfnContainer::~WfnContainer() {
  for (auto it = data.begin(); it != data.end(); ++it) {
    for (auto it_vec = it -> second -> begin(); it_vec != it -> second -> end(); ++it_vec) {
      if (*it_vec != nullptr) {
        delete (*it_vec);
        *it_vec = nullptr;
      }
    }
    if (it -> second != nullptr) {
      delete it -> second;
      it -> second = nullptr;
    }
  }
  for (auto it = l_idx.begin(); it != l_idx.end(); ++it) {
    delete it -> second;
    it -> second = nullptr;
  }
  for (auto it = block.begin(); it != block.end(); ++it) {
    delete it -> second;
    it -> second = nullptr;
  }
  for (auto it = factor.begin(); it != factor.end(); ++it) {
    delete it -> second;
    it -> second = nullptr;
  }
}

