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
public:
  SchmidtBasis(const Matrix& _core, const Matrix& _active): core(_core), active(_active) {};
  SchmidtBasis(const Matrix&, const vector<double>&, double);
  friend std::ostream& operator <<(std::ostream&, const SchmidtBasis&);

};


#endif
