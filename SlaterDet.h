#ifndef SLATER
#define SLATER

#include <iostream>
#include <tuple>
#include <vector>
#include <string>
#include <map>
#include "newmat10/newmatap.h"
#include "newmat10/newmatio.h"

using std::vector;
using std::string;
using std::tuple;
using std::map;
using std::get;

enum class Spin: int {up = 0, down = 1};

class SlaterDet;

class WfnContainer {
  map<int, vector<SlaterDet*>*> data;
  map<int, vector<int>*> block;
  map<int, vector<int>*> l_idx;  
  map<int, vector<double>*> factor;
  vector<int> keys;
  double sum;
  int M;
public:
  WfnContainer(int iM = 0): M(iM), sum(0) {}
  void push_back(SlaterDet*&, int, int); // return index
  const vector<SlaterDet*>& operator[](int); 
  const vector<int>& qp() const;
  int size(int) const;
  friend std::ostream& operator <<(std::ostream&, const WfnContainer&);
  void scale(SlaterDet&, double);
  void renormalize();
  void renormalize(map<int, int>&);
  void clear();
  const vector<int>& lq(int spin) const;
  const vector<int>& ld(int spin) const;  
  ~WfnContainer();
};

class SlaterDet {
  double factor;
  const vector<int> sites; // site index should start from 1
  Matrix* orb[2];
  tuple<int, double, ColumnVector*, Matrix*> make_emb(const Matrix&, int site);
public:
  SlaterDet(Matrix*&, Matrix*&, vector<int>&, double);
  int nelec(Spin) const;
  int nspin() const;
  double get_factor() const;
  string list_sites() const;
  friend std::ostream& operator <<(std::ostream&, const SlaterDet&);
  friend void WfnContainer::scale(SlaterDet&, double);
  friend double operator*(const SlaterDet& det1, const SlaterDet& det2);
  vector<SlaterDet*> Schmidt(int site, bool detail = false);
  ~SlaterDet();
};

#endif
