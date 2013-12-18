#ifndef MPS_OP
#define MPS_OP

#include "SpinQuantum.h"
#include "include.h"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <tuple>

using std::map;
using std::tuple;

MPS<Quantum> allocate(int nsites, int M, int M_total = -1, int M_min = -1, bool dim = false);

map<TVector<Quantum, 3>, TArray<double, 3>> allocate_blocks(const QSDArray<3,Quantum>& );

void insert_blocks(QSDArray<3,Quantum>&, map<TVector<Quantum, 3>, TArray<double, 3>>);

const TVector<Quantum, 3> make_qindex(int, int);

void save_site(const MPS<Quantum>&, int, const char *);

void load_site(MPS<Quantum>&, int,const char *);

double norm_on_disk(MPS<Quantum>&, const char*);

void normalize_on_disk(MPS<Quantum>&, const char*);

void compress_on_disk(MPS<Quantum>& mps,const MPS_DIRECTION& dir, int D, const char* filename, bool store = false);

tuple<SDArray<1>, Qshapes<Quantum>> Schmidt_on_disk(MPS<Quantum>& mps, int site, const char* filename);
#endif
