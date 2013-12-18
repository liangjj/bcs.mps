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

class SchmidtBasis {
private:
  Matrix core, active;
  vector<double> weight;
  double thr;
public:
  SchmidtBasis(const Matrix& _core, const Matrix& _active): core(_core), active(_active), weight(_active.Ncols(), 0.5), thr(0.) {};
  SchmidtBasis(const Matrix&, const vector<double>&, double, double);
  friend std::ostream& operator <<(std::ostream&, const SchmidtBasis&);
  int ncore() const {
    return core.Ncols();
  }
  int nactive() const {
    return  active.Ncols();
  }
};


#endif
