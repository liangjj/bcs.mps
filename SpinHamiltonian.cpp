#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::ostream;

#include "include.h"

/**
 * @param d local dimension
 * @param qp Qshapes object containing the local quantumnumbers on output, input is destroyed
 */
template<class Q>
void physical(int d,Qshapes<Q> &qp){

   qp.clear();

   int m = -d + 1;

   while(m < d){

      qp.push_back(Q(m));

      m += 2;

   }

}

template void physical<Quantum>(int d,Qshapes<Quantum> &);

