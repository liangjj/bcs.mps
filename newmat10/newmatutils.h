/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef NEWMATUTILS_HEADER_H
#define NEWMATUTILS_HEADER_H
#include "newmat.h"
#include <fstream>

namespace boost {
namespace serialization {

template<class Archive>
void serialize(Archive & ar, SymmetricMatrix& a, const unsigned int version)
{
  int dim = a.Nrows();
  ar & dim;
  if(dim != a.Nrows())
    a.ReSize(dim);
  for(int i=0;i<a.Storage();++i)
    ar & a.Store()[i];
}

template<class Archive>
void serialize(Archive & ar, DiagonalMatrix& a, const unsigned int version)
{
  int dim = a.Nrows();
  ar & dim;
  if(dim != a.Nrows())
    a.ReSize(dim);
  for(int i=0;i<a.Storage();++i)
    ar & a.Store()[i];
}

template<class Archive>
void serialize(Archive & ar, Matrix& a, const unsigned int version)
{
  int Nrs = a.Nrows();
  int Ncs = a.Ncols();
  ar & Nrs;
  ar & Ncs;

  if(a.Nrows() != Nrs && a.Ncols() != Ncs)
    {
      a.ReSize(Nrs, Ncs);
    }
  for(int i=0;i<a.Storage();++i) {
    ar & a.Store()[i];
  }
}



} }

#endif
