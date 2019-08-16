#include "Vector.h"
#include <math.h>

namespace tcods
{
   Vector :: Vector( void )
   : x( 0. ),
     y( 0. ),
     z( 0. )
   {}
   
   Vector :: Vector( double x0,
                     double y0,
                     double z0 )
   : x( x0 ),
     y( y0 ),
     z( z0 )
   {}
   
   double& Vector :: operator[]( const int& index )
   {
      return ( &x )[ index ];
   }
   
   const double& Vector :: operator[]( const int& index ) const
   {
      return ( &x )[ index ];
   }
   
   Vector Vector :: operator+( const Vector& v ) const
   {
      return Vector( x + v.x,
                     y + v.y,
                     z + v.z );
   }
   
   Vector Vector :: operator-( const Vector& v ) const
   {
      return Vector( x - v.x,
                     y - v.y,
                     z - v.z );
   }
   
   Vector Vector :: operator-( void ) const
   {
      return Vector( -x,
                     -y,
                     -z );
   }
   
   double Vector :: operator*( const Vector& v ) const
   {
      return x*v.x +
             y*v.y +
             z*v.z ;
   }
   
   Vector Vector :: operator^( const Vector& v ) const
   {
      return Vector( y*v.z - z*v.y,
                     z*v.x - x*v.z,
                     x*v.y - y*v.x );
   }
   
   Vector Vector :: operator*( const double& c ) const
   {
      return Vector( x*c,
                     y*c,
                     z*c );
   }
   
   Vector operator*( const double& c, const Vector& v )
   {
      return v*c;
   }
   
   Vector Vector :: operator/( const double& c ) const
   {
      return (*this) * ( 1./c );
   }
   
   void Vector :: operator+=( const Vector& v )
   {
      x += v.x;
      y += v.y;
      z += v.z;
   }
   
   void Vector :: operator-=( const Vector& v )
   {
      x -= v.x;
      y -= v.y;
      z -= v.z;
   }
   
   void Vector :: operator*=( const double& c )
   {
      x *= c;
      y *= c;
      z *= c;
   }
   
   void Vector :: operator/=( const double& c )
   {
      (*this) *= ( 1./c );
   }
   
   double Vector :: norm( void ) const
   {
      return sqrt( norm2());
   }
   
   double Vector :: norm2( void ) const
   {
      return (*this) * (*this);
   }
   
   void Vector :: normalize( void )
   {
      (*this) /= norm();
   }
   
   Vector Vector :: unit( void ) const
   {
      return (*this) / norm();
   }
   
   std::ostream& operator << (std::ostream& os, const Vector& o)
   {
      os << "[ "
         << o.x << " "
         << o.y << " "
         << o.z
         << " ]";
   
      return os;
   }
}

