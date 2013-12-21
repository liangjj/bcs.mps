#ifndef SCHMIDT
#define SCHMIDT

#include "include.h"
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "newmat10/newmatap.h"
#include "newmat10/newmatio.h"

using std::vector;
using std::string;
using std::map;

enum class Spin: int {up = 0, down = 1};

// forward declaration
uint choose(int, int);
class SchmidtBasis;

class ActiveSpaceIterator {
private:
  // number of sites, occupations alpha and beta
  int nsites, noccA, noccB;
  // pointer to the parent Schmidt basis
  const SchmidtBasis* basis;
  // stores all possible indices, alpha and beta spin
  vector<pair<uint, uint>> list;

  // private functions: internal conversion
  uint addr(const vector<bool>&) const;  // bits -> address
  vector<bool> bits(uint, Spin) const;         // address to bits

public:
  // constructor
  ActiveSpaceIterator(int, int, int, const SchmidtBasis*);

  // iterator behaviors
  std::pair<vector<bool>, vector<bool>> get_pair(int i) const;
  int size() const {
    return list.size();
  }

  // destructor
  ~ActiveSpaceIterator() {  basis = nullptr;}
};

class SchmidtBasis {
private:
  // matrices of core and active orbitals
  Matrix core, active;
  // weight of each active orbital
  vector<double> weight;
  // threshold to discard the Slater determinant
  double thr;
  // store iterators
  map<vector<int>, ActiveSpaceIterator> iterators;

public:
  // constructors
  SchmidtBasis(const Matrix& _core, const Matrix& _active): core(_core), active(_active), weight(_active.Ncols(), 0.5), thr(0.) {};
  SchmidtBasis(const Matrix&, const vector<double>&, double, double);

  // cout
  friend std::ostream& operator <<(std::ostream&, const SchmidtBasis&);

  // get properties
  int ncore() const {  return core.Ncols();}
  int nactive() const {  return active.Ncols();}
  int nsites() const { return core.Nrows();}
  double get_thr() const {  return thr;}
  double get_weight(const vector<bool>&) const;

  // iterators
  ActiveSpaceIterator iterator(int, int);

  // member acess
  const Matrix get_core() const {  return std::move(core);}
  const Matrix get_active() const {  return std::move(active);}
  const ColumnVector get_core(int n) const {  return std::move(core.Column(n));}
  const ColumnVector get_active(int n) const {  return std::move(active.Row(n));}

  // destructor
  ~SchmidtBasis() {}
};

class CoupledBasis {
private:
  // left and right basis of the site
  SchmidtBasis *lbasis, *rbasis;
  // Matrix elements between orbitals
  Matrix cc, ac, ca, aa, cs, as; // c - core, a - active, s - site
  // number of sites (left)
  int nsites;
  // left/right core/active space size
  int lc, rc, la, ra; 
  vector<int> ql, qp, qr;
  vector<int> dl, dp, dr;
  vector<vector<int>> block;

  // member function contract 1-particle part
  void contract1p();
  // calculate possible quantum number on left, right and physical indices
  void quantum_number();
  void dimensions();

public:
  // constructor
  CoupledBasis(SchmidtBasis&, SchmidtBasis&);
  
  // generate MPS site
  TVector<Quantum, 3> generate();

  // destructor
  ~CoupledBasis();
};
#endif
