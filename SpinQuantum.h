#ifndef SPINQUANTUM_H
#define SPINQUANTUM_H

#include <iostream>
#include <iomanip>
#include <boost/serialization/serialization.hpp>
/**
 * class which defines the quantumnumbers of the system, in this case spinquantumnumbers
 */
class SpinQuantum
{

   public:

      /**
       * empty constructor, set the spin to zero
       */
      SpinQuantum(){

         Sz = 0;

      }

      /**
       * constructor with input, set the spin to input
       * @param Sz_i input quantumnumber
       */
      SpinQuantum(int Sz_i){

         Sz = Sz_i;

      }

      /**
       * copy constructor
       */
      SpinQuantum(const SpinQuantum &qn_copy){

         Sz = qn_copy.gSz();

      }

      /**
       * @return the quantumnumber
       */
      int gSz() const {

         return Sz;

      }

      /**
       *  equality operator overloaded
       * @param qn_i input
       * @return true if input == *this
       */
      inline bool operator==(const SpinQuantum &qn_i) const { 

         return (Sz == qn_i.gSz());

      }

      /**
       * inequality operator overloaded
       * @param qn_i input
       * @return true if input != *this
       */
      inline bool operator!=(const SpinQuantum& qn_i) const {

         return Sz != qn_i.gSz(); 

      }

      /**
       * < comparison operator overloaded
       * @param qn_i input
       * @return true if *this < input
       */
      inline bool operator<(const SpinQuantum& qn_i) const {

         return Sz < qn_i.gSz(); 

      }

      /**
       * > comparison operator overloaded
       * @param qn_i input
       * @return true if *this > input
       */
      inline bool operator>(const SpinQuantum& qn_i) const { 

         return Sz > qn_i.gSz(); 

      }

      /**
       * operator acting on quantumnumbers
       * @param qn_i input
       * @return new SpinQuantum object with Sz = *this + input
       */
      inline SpinQuantum operator*(const SpinQuantum& qn_i) const {

         return SpinQuantum(Sz + qn_i.gSz());

      }

      /**
       * overload the + operator: basically makes a copy of the input object
       * @param input object q
       */
      friend SpinQuantum operator+ (const SpinQuantum& q) {

         return SpinQuantum(q.gSz()); 

      }

      /**
       * overload the - operator: basically makes a copy of the input object with negative sign
       * @param input object q
       */
      friend SpinQuantum operator-(const SpinQuantum& q) { 

         return SpinQuantum(-q.gSz()); 

      }

      /**
       * overload output stream operator
       */
      friend std::ostream& operator<< (std::ostream& ost, const SpinQuantum& q) {

         ost << "(" << std::setw(2) << q.Sz << ")";

         return ost;

      }

      inline bool parity() const { return false; }

      /**
       * @return a SpinQuantum object initialized on zero
       */
      const static SpinQuantum zero() {

         return SpinQuantum(0); 

      }

   private:
      friend class boost::serialization::access;
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version) { ar & Sz; }
      //! the z component of the spin
      int Sz;

};

#endif // SPINQUANTUM_H
