#ifndef SCHMIDT
#define SCHMIDT

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "newmat10/newmatap.h"
#include "newmat10/newmatio.h"

using std::vector;
using std::string;

enum class Spin: int {up = 0, down = 1};

uint choose(int, int);

class ActiveSpaceIterator {
private:
  int nsites, nocc;
  uint max, ptr;

  uint addr(const vector<bool>&) const;
  vector<bool> bits(uint) const;
public:
  ActiveSpaceIterator(int, int);
  void reset() {
    ptr = 0;
  }
  vector<bool> next();
  bool end();
};

class SchmidtBasis {
private:
  Matrix core, active;
  vector<double> weight;
  double thr;
  vector<uint> addr_ptr;
public:
  // constructors
  SchmidtBasis(const Matrix& _core, const Matrix& _active): core(_core), active(_active), weight(_active.Ncols(), 0.5), thr(0.), addr_ptr(_active.Ncols(), 0) {};
  SchmidtBasis(const Matrix&, const vector<double>&, double, double);
  // cout
  friend std::ostream& operator <<(std::ostream&, const SchmidtBasis&);
  // get properties
  int ncore() const {
    return core.Ncols();
  }
  int nactive() const {
    return active.Ncols();
  }
  int nsites() const {
   return core.Nrows();
  }
  double get_thr() const {
    return thr;
  }
  // member acess
  const Matrix get_core() const {
    return std::move(core);
  }
  const Matrix get_active() const {
    return std::move(active);
  }
  const ColumnVector get_core(int n) const {
    return std::move(core.Column(n));
  }
  const ColumnVector get_active(int n) const {
    return std::move(active.Row(n));
  }
  // destructor
  ~SchmidtBasis() {}
};

class CoupledBasis {
private:
  const SchmidtBasis *lbasis, *rbasis;
  Matrix cc, ac, ca, aa, cs, as; // c - core, a - active, s - site
  int nsites;
  int lc, rc, la, ra; 
  void contract1p();
public:
  CoupledBasis(const SchmidtBasis&, const SchmidtBasis&);
  ~CoupledBasis();
};
#endif
