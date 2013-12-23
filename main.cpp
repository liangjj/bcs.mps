#include <cmath>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <string>
#include <sstream>      // std::istringstream

#include "include.h"
#include "schmidt.h"
#include "mps_op.h"
#include "newmat10/newmatap.h"
#include "newmat10/newmatio.h"
#include "boost/filesystem.hpp"

using std::cout;
using std::endl;
using std::string;
using std::istringstream;
using std::ifstream;

namespace btas { typedef SpinQuantum Quantum; }; // Defined as default quantum number class

const string mps_dir = "/home/boxiao/mps/mps.out";
const string temp_dir = "/scratch/gpfs/boxiao/MPSTemp";

using namespace btas;
using namespace mpsxx;

void permute(Matrix& orbs, vector<int>& order) {
  Matrix orbs_i = orbs;
  for (int i = 0; i < order.size(); ++i) {
    orbs.Row(i+1) = orbs_i.Row(order[i]);
  }
}

Matrix read_input(char* file) {
  string line;
  ifstream in(file);
  istringstream is(line);
  std::getline(in, line);  // the first line is an explanation of the file.
    
  int nsites, norbs;
  std::getline(in, line);
  is.str(line);
  is >> nsites >> norbs;
  std::getline(in, line);
  
  vector<int> order;
  if (line.compare("default") == 0) {
    for (int i = 0; i < nsites; ++i) {
      order.push_back(i+1);
    }
  } else {
    int temp;
    is.clear();
    is.str(line);
    for (int i = 0; i < nsites; ++i) {
      is >> temp;
      order.push_back(temp);
    }
    assert(order.size() == nsites);
  }
  Matrix orbs(nsites, norbs);
  double temp;
  for (int i = 0; i < norbs; i++) {
    for (int j = 0; j < nsites; j++) {
      in >> temp;
      orbs(j+1, i+1) = temp;
    }
  }
  permute(orbs, order);
  return std::move(orbs);
}

string mktmpdir(const string& prefix) {
  char* temp = new char[prefix.size() + 11];
  std::strcpy(temp, (prefix + "/tmpXXXXXX").c_str());
  mkdtemp(temp);
  // create this folder
  cout << "MPS Temporary Directory " << temp << endl;
  return string(temp);
}

void banner() {
  cout << "-----------------------------------------------------------------------\n";
  cout << "                           G P S - M P S                               \n";
  cout << "   (Gutzwiller Projection of Single Slater Determinants through MPS)   \n\n";
  cout << "                           Bo-Xiao Zheng                               \n";
  cout << "-----------------------------------------------------------------------\n\n";
}

void natural_orbs(const SymmetricMatrix& rdm, vector<double>& occs, Matrix& orbs) {
  occs.clear();
  DiagonalMatrix D;
  Jacobi(rdm, D, orbs);
  for (int i = 0; i < D.Nrows(); ++i) {
    occs.push_back(D(i+1, i+1));
  }
}

int main(int argc, char* argv[]){
  cout.setf(std::ios::fixed, std::ios::floatfield);
  cout.precision(10);

  assert(argc > 1);
  double thr1p = pow(10., -7.);
  double thrnp = pow(10., -8.);
  int M = 0;
  if (argc > 2) {
    thrnp = pow(10., -atoi(argv[2]));
  }
  if (argc > 3) {
    M = atoi(argv[3]);
  }
  banner();
  string mps_temp = mktmpdir(temp_dir);
  // read input file
  Matrix coefs = read_input(argv[1]);
  int nsites = coefs.Nrows();
  int norbs = coefs.Ncols();
  // density matrix
  SymmetricMatrix rdm;
  rdm << coefs * coefs.t();
  vector<double> occs(norbs, 1.);
  SchmidtBasis lbasis(coefs, occs, thr1p, thrnp);
  // prepare MPS
  MPS<Quantum> A(nsites);

  for (int site = 0; site < nsites; ++site) {
    // first build right basis
    SymmetricMatrix prdm;
    prdm << rdm.SubMatrix(site+2, nsites, site+2, nsites);
    Matrix natorbs;
    natural_orbs(prdm, occs, natorbs);
    SchmidtBasis rbasis(natorbs, occs, thr1p, thrnp);
    cout << "Site: " << site+1 << endl;
    cout << "Left\n" << lbasis << endl;
    cout << "Right\n" << rbasis << endl;    
    CoupledBasis basis_pair(lbasis, rbasis);
    // do some thing
    A[site] = basis_pair.generate();
    save_site(A, site, mps_temp.c_str());
    A[site].clear();
    lbasis = std::move(rbasis);
  }

  compress_on_disk(A, MPS_DIRECTION::Right, M, mps_temp.c_str(), true);  
  cout << "\nnow calculate entanglement entropy\n";
  auto tup = Schmidt_on_disk(A, -1, mps_temp.c_str());

  SDArray<1> sc = std::get<0>(tup);
  Qshapes<Quantum> sq = std::get<1>(tup);
  auto iter = sc.begin();  
  map<int, vector<double>> coef;
  for (int i = 0; i < sq.size(); i++) {
    int sp = sq[i].gSz();
    coef[sp];
    for (auto it_d = iter -> second -> begin(); it_d != iter -> second -> end(); ++it_d) {
      coef[sp].push_back(-2.*log(*it_d));
    }
    std::sort(coef[sp].begin(), coef[sp].end());
    ++iter;    
  }

  cout << "Entanglement Spectra:\n";
  for (auto it = coef.begin(); it != coef.end(); it++) {
    cout << "Section S=" << it -> first << endl;
    for (auto it_d = it -> second.begin(); it_d < it -> second.end(); ++it_d) {
      cout << *it_d << "\t";
    }
    cout << endl;
  }

  boost::filesystem::path to_remove(mps_temp);
  boost::filesystem::remove_all(to_remove);
  return 0;
  /*
  ofstream ofs((mps_dir+"/es").c_str());
  ofs.setf(std::ios::fixed, std::ios::floatfield);
  ofs.precision(10);
  for (auto it = coef.begin(); it != coef.end(); it++) {
    ofs << "Section S=" << it -> first;
    for (auto it_d = it -> second.begin(); it_d < it -> second.end(); ++it_d) {
      ofs << *it_d << "\t";
    }
    ofs << endl;
  }
  ofs.close();
  */
}
