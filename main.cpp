#include <cmath>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <string>
#include <sstream>      // std::istringstream
#include <omp.h>

#include "include.h"
#include "SlaterDet.h"
#include "newmat10/newmatap.h"
#include "newmat10/newmatio.h"
#include "mps_op.h"
#include "boost/filesystem.hpp"

using std::cout;
using std::endl;
using std::string;
using std::istringstream;
using std::ifstream;

namespace btas { typedef SpinQuantum Quantum; }; // Defined as default quantum number class

const string mps_dir = "/home/boxiao/mps/mps.out";

using namespace btas;
using namespace mpsxx;

/*
void is_canonical(QSDArray<3,Quantum>& mat, MPS_DIRECTION dir) {
    enum{p,q,r,s};
    QSDArray<2,Quantum> loc;
    if (dir == MPS_DIRECTION::Left) {
      QSDindexed_contract (1.0, mat, shape(p,q,r), mat.conjugate(), shape(p,q,s), 0.0, loc, shape(r,s));
    } else {
      QSDindexed_contract (1.0, mat, shape(p,q,r), mat.conjugate(), shape(s,q,r), 0.0, loc, shape(p,s));
    }
    cout << loc;
    // need more to finish
}
*/

int read_input(char* file, vector<int>& order, Matrix*& orbs) {
  string line;
  ifstream in(file);
  istringstream is(line);

  std::getline(in, line);  // the first line is an explanation of the file.
    
  int nsites, norbs;
  std::getline(in, line);
  is.str(line);
  is >> nsites >> norbs;
  std::getline(in, line);
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

  orbs = new Matrix(nsites, norbs);
  double temp;
  for (int i = 0; i < norbs; i++) {
    for (int j = 0; j < nsites; j++) {
      in >> temp;
      (*orbs)(j+1, i+1) = temp;
    }
  }
  return 0;
}

string mktmpdir(string prefix) {
  char* temp = new char[prefix.size() + 11];
  std::strcpy(temp, (prefix + "/tmpXXXXXX").c_str());
  mkdtemp(temp);
  // create this folder
  cout << temp << endl;
  return string(temp);
}

void permute(Matrix*& orbs, vector<int>& order) {
  Matrix* orbs_i = orbs;
  orbs = new Matrix();
  orbs -> ReSize(*orbs_i);
  for (int i = 0; i < order.size(); ++i) {
    orbs -> Row(i+1) = orbs_i -> Row(order[i]);
  }
  delete orbs_i;
}

int main(int argc, char* argv[]){
  cout.setf(std::ios::fixed, std::ios::floatfield);
  cout.precision(10);

  assert(argc > 2);
  int M = atoi(argv[2]);
  cout << "M = " << M << "\n";
  string mps_temp = mktmpdir("/scratch/gpfs/boxiao/MPSTemp");
  // read input file
  // include eigenvectors (only RHF: number of sites, number of electrons per spin)
  // site order
  vector<int> order;
  Matrix* orbs;
  read_input(argv[1], order, orbs);
  permute(orbs, order);
  // creat the first Slater determinant object (whole system)
  vector<int> sites;
  for (int i = 0; i < (*orbs).Nrows(); ++i) {
    sites.push_back(i+1);
  }
  int nsites = sites.size();

  Matrix* orbs1 = new Matrix(*orbs);
  SlaterDet* Det = new SlaterDet(orbs, orbs1, sites, 1.);
  WfnContainer *wfns;
  wfns = new WfnContainer(M);
  wfns -> push_back(Det, 0, 0); // wfns takes over the slater determinant and Det repoints to null.

  // prepare MPS
  MPS<Quantum> A = allocate(nsites, M, false);
  for (int i = 0; i < nsites-1; ++i) {
    int site = i+1;
    cout << "site " << i << "  index " << order[i] << endl;
    WfnContainer* new_wfns = new WfnContainer(M);
    
    # pragma omp parallel default(none) shared(wfns, new_wfns, site)
    {
    # pragma omp for schedule(guided, 1)
    for (int ispin = 0; ispin < wfns -> qp().size(); ++ispin) {
      int spin = wfns -> qp()[ispin];
    //for (int spin: wfns -> qp()) { // spin is left index
      for (int j = 0; j < (*wfns)[spin].size(); ++j) {
        SlaterDet* det = (*wfns)[spin][j];
        auto wfn_ud = det -> Schmidt(site);
        if (wfn_ud[0] != nullptr) {
          # pragma omp critical
          new_wfns -> push_back(wfn_ud[0], spin, j);
        }
        if (wfn_ud[1] != nullptr) {
          # pragma omp critical          
          new_wfns -> push_back(wfn_ud[1], spin, j);
        }
      }
    }
    }
    delete wfns;
    new_wfns -> renormalize(); // will 1. fix dimension if exceed M 2. normalize the states
    //cout << *new_wfns << endl;
    wfns = new_wfns;

    // allocate the blocks
    // first get the size of index space
    auto blocks = allocate_blocks(A[i]);
    // now fill in numbers
    for (int ispin = 0; ispin < wfns -> qp().size(); ++ispin) {
      int spin = wfns -> qp()[ispin];
      for (int j = 0; j < wfns -> size(spin); ++j) {
        (blocks[make_qindex(wfns -> lq(spin)[j], spin)])(wfns -> ld(spin)[j], 0, j) = 1.;
      }
    }
    // set blocks
    insert_blocks(A[i], blocks);
    save_site(A, i, mps_temp.c_str());
    A[i].clear();
  }

  // a little bit different for the last site
  auto blocks = allocate_blocks(A[nsites-1]);
  // now fill in numbers directly
  cout << "site " << nsites-1 << "  index " << order[nsites-1] << "(last)" << endl;
  for (int spin: wfns -> qp()) {
    for (int j = 0; j < wfns -> size(spin); ++j) {
      (blocks[make_qindex(spin, 0)])(j, 0, 0) = (*wfns)[spin][j] -> get_factor();
    }   
  }
  delete wfns;  
  // set blocks 
  insert_blocks(A[nsites-1], blocks);
  save_site(A, nsites-1, mps_temp.c_str());
  A[nsites-1].clear();
  // now right-canonicalize it
  cout << "\ncompress the state\n";
  compress_on_disk(A, MPS_DIRECTION::Right, M, mps_temp.c_str(), true);

  cout << "\nnow calculate entanglement entropy\n";
  int cut = -1;
  if (argc > 3) {
    cut = atoi(argv[3]);
  }
  auto tup = Schmidt_on_disk(A, cut, mps_temp.c_str());

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
  boost::filesystem::path to_remove(mps_temp);
  boost::filesystem::remove_all(to_remove);
  return 0;
}
