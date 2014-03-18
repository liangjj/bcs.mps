#include "schmidt.h"
#include <cassert>
#include <boost/mpi.hpp>

namespace mpi = boost::mpi;

QSDArray<3, Quantum> CoupledBasis::generate() {
  QSDArray<3, Quantum> A;
  // build qshape, dshape basically transform from the data we already have
  TVector<Qshapes<Quantum>, 3> qshape;
  TVector<Dshapes, 3> dshape;
  for (int i = 0; i < ql.size(); ++i) {
    qshape[0].push_back(Quantum(ql[i]));
    dshape[0].push_back(dl[i]);
  }
  for (int i = 0; i < qp.size(); ++i) {
    qshape[1].push_back(Quantum(qp[i]));
    dshape[1].push_back(dp[i]);
  }
  for (int i = 0; i < qr.size(); ++i) {
    qshape[2].push_back(-Quantum(qr[i]));
    dshape[2].push_back(dr[i]);
  }

  A.resize(Quantum::zero(), qshape, dshape, false);
  mpi::communicator world;

  if (world.rank() == 0) {
    //cout << "qshape = " << A.qshape() << endl;
    //cout << "dshape = " << A.dshape() << endl;
    cout << "nblock = " << block.size() << endl;
  }
  
  // now generate these blocks
  for (int i = 0; i < block.size(); ++i) {
    if (world.rank() == 0) {
      cout << "block " << i;
    }
    // get information of each block
    Spin s = (qp[block[i][1]] == -1) ? (Spin::down) : (Spin::up);
    int nl = nsites - ql[block[i][0]] - lc;
    int nr = nsites - 1 - qr[block[i][2]] - rc;
    auto iter_l = lbasis -> iterator(nl);
    auto iter_r = rbasis -> iterator(nr);

    // total size of the block
    size_t size = dl[block[i][0]];
    size *=  dr[block[i][2]];

    // declarations
    IVector<3> stride;
    DArray<3> dense;

    if (world.rank() == 0) {
      A.reserve(block[i]);
      dense.reference(*(A.find(block[i]) -> second));
      stride = dense.stride();
      cout << "\t" << "nelements = " << stride[0] << endl;
    }
    broadcast(world, stride, 0);

    vector<size_t> size_procs(world.size(), size / world.size());
    vector<size_t> shift_procs(world.size(), 0);    
    for (int j = 0; j < size % world.size(); ++j) {
      size_procs[j] += 1;
    }
    for (int j = 1; j < world.size(); ++j) {
      shift_procs[j] = shift_procs[j-1] + size_procs[j-1];
    }
    
    vector<double> my_array;
    vector<double>::iterator my_it;

    if (world.rank() == 0) {
      my_it = dense.begin();
    } else {
      my_array.resize(size_procs[world.rank()]);
      my_it = my_array.begin();
    }

    Matrix workspace;
    for (size_t j = shift_procs[world.rank()]; j <  shift_procs[world.rank()] + size_procs[world.rank()]; ++j) {
      int idx_l = j / stride[0];
      int idx_r = (j % stride[1]) / stride[2];
      *my_it = overlap(iter_l.get_config(idx_l), iter_r.get_config(idx_r), s, nl, nr, workspace);
      ++my_it;
    }

    if (world.rank() == 0) {
      MPI_Request req[world.size()-1];
      for (int j = 1; j < world.size(); ++j) {
        MPI_Irecv(dense.data() + shift_procs[j], size_procs[j], MPI_DOUBLE, j, j, world, &req[j-1]);
      }
      MPI_Waitall(world.size()-1, req, MPI_STATUS_IGNORE);
    } else {
      MPI_Request req;
      MPI_Isend(my_array.data(), size_procs[world.rank()], MPI_DOUBLE, 0, world.rank(), world, &req);
      MPI_Wait(&req, MPI_STATUS_IGNORE);
    }
  }
  cout << endl;
  return std::move(A);
}
