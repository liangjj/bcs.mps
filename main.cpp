#include <cmath>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <string>
#include <sstream>      // std::istringstream
#include <boost/mpi.hpp>

#include "include.h"
#include "schmidt.h"
#include "mps_op.h"
#include "newmat10/newmatap.h"
#include "newmat10/newmatio.h"
#include "newmat10/newmatutils.h"

#define BOOST_NO_CXX11_SCOPED_ENUMS
#include "boost/filesystem.hpp"
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/string.hpp>
#include <boost/algorithm/string.hpp>

using std::cout;
using std::endl;
using std::string;
using std::istringstream;
using std::ifstream;
namespace mpi = boost::mpi;
using boost::trim;
using boost::is_any_of;

namespace btas { typedef SpinQuantum Quantum; }; // Defined as default quantum number class

const string mps_dir = "/home/boxiao/mps/mps.out";
const string temp_dir = "/scratch/boxiao/MPSTemp";

using namespace btas;
using namespace mpsxx;

void permute(Matrix& orbs, const vector<int>& order) {
  Matrix orbs_i = orbs;
  for (int i = 0; i < order.size(); ++i) {
    orbs.Row(i+1) = orbs_i.Row(order[i]);
  }
}

void read_config(string file, double& thr1p, double& thrnp, int& M, bool& calc_spectra, bool& savemps) {
  string line;
  ifstream in(file.c_str());

  vector<string> tokens;
  while(std::getline(in, line)) {
    if (line.find('!') != string::npos) {
      boost::erase_tail(line, line.size() - line.find('!'));
    }
    trim(line);
    boost::split(tokens, line, is_any_of("\t ="), boost::token_compress_on);
    if (tokens[0].size() == 0) {
      continue;
    }
    if (boost::iequals(tokens[0], "thr1p")) {
      thr1p = atof(tokens[1].c_str());
    } else if (boost::iequals(tokens[0], "thrnp")) {
      thrnp = atof(tokens[1].c_str());
    } else if (boost::iequals(tokens[0], "M")) {
      M = atoi(tokens[1].c_str());
    } else if (boost::iequals(tokens[0], "nospectra")) {
      calc_spectra = false;
    } else if (boost::iequals(tokens[0], "savemps")) {
      savemps = true;
    } else {
      cout << "\nUnrecognized option in config file:" << endl;
      cout << "\t" << tokens[0] << endl;
      abort();
    }
  }
}

Matrix read_orbitals(string file) {
  string line;
  ifstream in(file.c_str());
  istringstream is(line);
  std::getline(in, line);  // the first line is an explanation of the file.
    
  int nsites;
  std::getline(in, line);
  is.str(line);
  is >> nsites;
  std::getline(in, line);
  
  vector<int> order;
  
  // input site order
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

  Matrix u(nsites, nsites), v(nsites, nsites); // u and v
  double temp;
  for (int i = 0; i < nsites; i++) {
    for (int j = 0; j < nsites; j++) {
      in >> temp;
      u(j+1, i+1) = temp;
    }
    for (int j = 0; j < nsites; j++) {
      in >> temp;
      v(j+1, i+1) = temp;
    }
  }
  permute(u, order);
  permute(v, order);
  Matrix orbs(nsites*2, nsites);
  orbs.Rows(1, nsites) = u;
  orbs.Rows(nsites+1, nsites*2) = v;
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

void natural_orbs(const SymmetricMatrix& rho, const SymmetricMatrix& kappa, vector<double>& occs, Matrix& orbs) {
  occs.clear();
  int nsites = rho.Nrows();
  
  DiagonalMatrix D, eye(nsites);
  Matrix V;
  eye = 1.;

  SymmetricMatrix grdm(nsites*2);
  grdm.SubMatrix(1, nsites, 1, nsites) = eye - rho;
  grdm.SubMatrix(nsites+1, nsites*2, 1, nsites) = kappa;
  grdm.SubMatrix(1, nsites, nsites+1, nsites*2) = kappa;
  grdm.SubMatrix(nsites+1, nsites*2, nsites+1, nsites*2) = rho;

  Jacobi(grdm, D, orbs);
  for (int i = 0; i < D.Nrows(); ++i) {
    occs.push_back(D(i+1, i+1));
  }
}

int main(int argc, char* argv[]){
  mpi::environment env(argc, argv);
  cout.setf(std::ios::fixed, std::ios::floatfield);
  cout.precision(10);

  assert(argc > 1);

  // set defaults
  double thr1p = 1e-7, thrnp = 1e-8;
  int M = 0;
  bool calc_spectra = true, savemps = false;
  
  string mps_temp;
  Matrix coefs;

  mpi::communicator world;  
  if (world.rank() == 0) {
    banner();
    mps_temp = mktmpdir(temp_dir);    
    // read input file
    read_config(string(argv[1]) + "/config.in", thr1p, thrnp, M, calc_spectra, savemps);
    coefs = read_orbitals(string(argv[1]) + "/orbitals.in"); // coefficients: upper-half u, lower-half v
  }
  broadcast(world, mps_temp, 0);
  broadcast(world, thr1p, 0);
  broadcast(world, thrnp, 0);
  broadcast(world, M, 0);
  broadcast(world, coefs, 0);
  broadcast(world, savemps, 0);
  
  int nsites = coefs.Nrows()/2;
  int norbs = nsites;
  // density matrix
  SymmetricMatrix rho, kappa;
  rho << coefs.Rows(nsites+1, nsites*2) * coefs.Rows(nsites+1, nsites*2).t(); // \rho = VV^t
  kappa << coefs.Rows(1, nsites) * coefs.Rows(nsites+1, nsites*2).t();  // \kappa = UV^t
  vector<double> occs(nsites, 1.);
  SchmidtBasis lbasis(coefs, occs, thr1p, thrnp);
  // prepare MPS
  MPS<Quantum> A(nsites);

  for (int site = 0; site < nsites; ++site) {
    // first build right basis
    SymmetricMatrix prho, pkappa;
    prho << rho.SubMatrix(site+2, nsites, site+2, nsites);
    pkappa << kappa.SubMatrix(site+2, nsites, site+2, nsites);
    Matrix natorbs;
    natural_orbs(prho, pkappa, occs, natorbs);
    SchmidtBasis rbasis(natorbs, occs, thr1p, thrnp);

    if (world.rank() == 0) {    
      cout << "Site: " << site+1 << endl;
      cout << "Left\n" << lbasis << endl;
      cout << "Right\n" << rbasis << endl;
    }

    CoupledBasis basis_pair(lbasis, rbasis);
    A[site] = basis_pair.generate();
    if (world.rank() == 0) {
      save_site(A, site, mps_temp.c_str());
      A[site].clear();
    }
    lbasis = std::move(rbasis);
  }

  if (world.rank() == 0) {
    compress_on_disk(A, MPS_DIRECTION::Right, M, mps_temp.c_str(), true);

    if (calc_spectra) {
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

      ofstream fspec((string(argv[1]) + "/spectra.out").c_str());
      for (auto it = coef.begin(); it != coef.end(); it++) {
        fspec << "Section S=" << it -> first << endl;
        for (auto it_d = it -> second.begin(); it_d < it -> second.end(); ++it_d) {
          fspec << *it_d << "\t";
        }
        fspec << endl;
      }
      fspec.close();
    }
  }
  
  world.barrier();
  boost::filesystem::path mps_tmp_store(mps_temp);
  if (savemps) {
    if (world.rank() == 0) {
      cout << "\nnow save mps files" << endl;
    }
    boost::filesystem::path mps_store(string(argv[1]) + "/mps.out");
    boost::filesystem::create_directory(mps_store);
    for (int i = 0; i < nsites; ++i) {
      if (world.rank() == i % world.size()) {
        string filename = std::to_string(i) + ".mps";
        boost::filesystem::path p1(mps_temp + "/" + filename);
        boost::filesystem::path p2(string(argv[1]) + "/mps.out/" + filename);
        boost::filesystem::copy_file(p1, p2, boost::filesystem::copy_option::overwrite_if_exists);
      }
    }
  }

  world.barrier();

  if (world.rank() == 0) {
    boost::filesystem::remove_all(mps_tmp_store);
  }
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
