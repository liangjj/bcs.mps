#ifndef SCHMIDT
#define SCHMIDT

#include "include.h"
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "newmat10/newmatap.h"
#include "newmat10/newmatio.h"
#include "newmat10/newmatutils.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

using std::vector;
using std::string;
using std::map;

using namespace btas;
using namespace mpsxx;

enum class Spin: int {up = 0, down = 1};

// forward declaration
uint choose(int, int);
class SchmidtBasis;

class ActiveSpaceIterator {
private:
  // number of sites, number of excititations needed
  int nsites, nex;
  // pointer to the parent Schmidt basis
  const SchmidtBasis* basis;
  // stores all possible indices, first and second half (temp)
  vector<vector<bool>> list;

  // private functions: internal conversion
  uint addr(const vector<bool>&) const;  // bits -> address
  vector<bool> bits(uint address, int nocc = -1, bool half = false) const;  // address to bits
  
  bool initialized;
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version) {
    ar & nsites;
    ar & nex;
    ar & list;
    ar & initialized;
  }
public:
  // constructor
  void initialize(int, int, const SchmidtBasis*); 
  ActiveSpaceIterator(int _nsites, int _nex, const SchmidtBasis* _basis) {  initialize(_nsites, _nex, _basis);}
  ActiveSpaceIterator(): initialized(false) {};

  // iterator behaviors
  vector<bool> get_config(int i) const {  return std::move(list[i]);}
  int size() const {  return list.size();}
  void set_basis(const SchmidtBasis* _basis) {  basis = _basis;}
  bool is_initialized() {  return initialized;}

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
  map<int, ActiveSpaceIterator> iterators;
  map<int, int> iterator_map;

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version) {
    ar & core;
    ar & active;
    ar & weight;
    ar & thr;
    ar & iterators & iterator_map;
    ar & dims & quantums;

    for (auto iter = iterators.begin(); iter != iterators.end(); ++iter) {
      iter -> second.set_basis(this);
    }
  }
public:
  // constructors
  SchmidtBasis(const Matrix& _core, const Matrix& _active): core(_core), active(_active), weight(_active.Ncols(), 0.5), thr(0.) {};
  SchmidtBasis(const Matrix&, const vector<double>&, double, double);
  void dimensions();

  // cout
  friend std::ostream& operator <<(std::ostream&, const SchmidtBasis&);
  vector<int> dims, quantums;  

  // get properties
  int ncore() const {  return core.Ncols();}
  int nactive() const {  return active.Ncols();}
  int nsites() const {  return core.Nrows()/2;}
  double get_thr() const {  return thr;}
  double get_weight(const vector<bool>&, int shift = 0) const;

  // used for mpi
  void broadcast_iterators();
  int iterator_on_rank(int nex) {return iterator_map[nex];}
  int q2nex(int q) const {  return -q + nsites() - ncore();} 

  // iterators
  ActiveSpaceIterator iterator(int nex, bool init = true);

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
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version) {
    ar & (*lbasis) & (*rbasis);
    ar & cc & ac & ca & aa & cs & as;
    ar & nsites & lc & rc & la & ra;
    ar & ql & qp & qr;
    ar & dl & dp & dr;
    ar & block;
  }
  // left and right basis of the site
  SchmidtBasis *lbasis, *rbasis;
  // Matrix elements between orbitals
  Matrix cc, ac, ca, aa, cs, as; // c - core, a - active, s - site (two columns a_{i\up}^\dagger a_{i\down})
  // number of sites (left)
  int nsites;
  // left/right core/active space size
  int lc, rc, la, ra; 
  vector<int> ql, qp, qr;
  vector<int> dl, dp, dr;
  vector<IVector<3>> block;

  // member function contract 1-particle part
  void contract1p();
  // compute matrix elements: given active space bits left/right, and the onsite state
  // calculat overlap
  // sign convention : <l_vac|<l_core|<l_act|s_i>|r_act>r_core>|r_vac>, resulted basis becomes |s_1,s_2,...>
  double overlap(const vector<bool>&, const vector<bool>&, Spin, int, int, Matrix&) const;

public:
  // constructor
  CoupledBasis(SchmidtBasis&, SchmidtBasis&);
  
  // generate MPS site
  QSDArray<3, Quantum> generate();

  // destructor
  ~CoupledBasis();
};
#endif
