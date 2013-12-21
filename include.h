//nog enkele definities:
#ifndef INCLUDE_H_
#define INCLUDE_H_

#include "SpinQuantum.h"
namespace btas { typedef SpinQuantum Quantum; };

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <btas/blas_cxx_interface.h>

#include <btas/TVector.h>

#include <btas/DENSE/DArray.h>
#include <btas/QSPARSE/QSDArray.h>
#include <btas/QSPARSE/QSDcontract.h>

#include <MPSblas.h>

#include "SpinHamiltonian.h"

#endif
