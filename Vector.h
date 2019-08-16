//============================================================
// Vector.h
// 
// Vector represents an element of Euclidean 3-space.
//

#ifndef TCODS_VECTOR_H
#define TCODS_VECTOR_H

#include <iostream>

namespace tcods
{
   class Vector
   {
      public:
         Vector();                               // default constructor -- all coordinates are zero
         Vector( double x, double y, double z ); // construct w/ specified coordinates

         double& operator[] ( const int& index );             // returns reference to specified coordinate (0-based indexing: x, y, z)
         const double& operator[] ( const int& index ) const; // returns const reference to specified coordinate (0-based indexing: x, y, z)
         Vector operator+  ( const Vector& v ) const;         // vector addition
         Vector operator-  ( const Vector& v ) const;         // vector subtraction
         Vector operator-  ( void ) const;                    // negation
         double operator*  ( const Vector& v ) const;         // dot product
         Vector operator^  ( const Vector& v ) const;         // cross product
         Vector operator*  ( const double& c ) const;         // scalar product
         Vector operator/  ( const double& c ) const;         // scalar division
         void   operator+= ( const Vector& v );               // vector addition / assignment
         void   operator-= ( const Vector& v );               // vector subtraction / assignment
         void   operator*= ( const double& c );               // scalar product / assignment
         void   operator/= ( const double& c );               // scalar division / assignment
         double  norm( void ) const;                          // returns Euclidean length
         double norm2( void ) const;                          // returns Euclidean length squared
         Vector  unit( void ) const;                          // returns vector divided by norm
         void normalize( void );                              // divides by norm

         double x, y, z; // coordinates
   };

   Vector operator* ( const double& c, const Vector& v );         // scalar product
   std::ostream& operator << (std::ostream& os, const Vector& o); // prints coordinates
}

#endif

