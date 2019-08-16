#include "Connection.h"
#include <cmath>
#include <iostream>

using namespace std;
using namespace cm;

namespace tcods
{
   Connection :: Connection( Mesh& _mesh )
   : mesh( _mesh ),
     K( common ),
     b( common ),
     QR( NULL )
   {
      build();
   }
   
   Connection :: ~Connection( void )
   {
      if( QR != NULL ) SuiteSparseQR_free<double>( &QR, common );
   }
   
   double sgn( double x )
   {
      if( x > 0. ) return  1.;
      if( x < 0. ) return -1.;
      return 0.;
   }
   
   void Connection :: build( void )
   // build the system of constraint equations
   {
      // construct a list of basis cycles, represented by a vector
      // of oriented edges (really: halfedges)
      std::vector<Mesh::Cycle> basisCycles;
   
      // append contractible basis cycles
      append1RingBases( basisCycles );
   
      unsigned int nContractibleCycles = static_cast<unsigned int>(basisCycles.size());
   
      // append noncontractible basis cycles
      mesh.appendDualGenerators( basisCycles );
   
      nBasisCycles = static_cast<unsigned int>(basisCycles.size());
   
      // append derivative constraints used to specify directional constraints in faces
      vector<double> directionalHolonomies;
      mesh.appendDirectionalConstraints( basisCycles, directionalHolonomies );
   
      // build constraint matrix
      Sparse A( common, static_cast<int>(mesh.edges.size()), static_cast<int>(basisCycles.size()) );
      buildCycleMatrix( A, basisCycles );
      applyCotanWeights( A );
   
      // prefactor [Q,R,E] = qr( A' )
      if( QR != NULL ) { SuiteSparseQR_free<double>( &QR, common ); QR = NULL; }
      QR = SuiteSparseQR_factorize <double> (7, -2., *A, common);
      
      // compute Riemannian holonomy of basis cycles
      K = Dense( common, static_cast<int>(basisCycles.size()), 1 );
      b = Dense( common, static_cast<int>(basisCycles.size()), 1 );
      generatorOnBoundary.resize( mesh.nGenerators() );
      for( unsigned int i = 0; i < nContractibleCycles; i++ )
      {
         K( i ) = -mesh.vertex( i )->defect();
      }
      for( unsigned int i = nContractibleCycles;
                        i < nContractibleCycles + mesh.nGenerators();
                        i++ )
      {
         if( Mesh::isDualBoundaryLoop( basisCycles[i] ))
         {
            K( i ) = -mesh.boundaryLoopCurvature( basisCycles[i] );
            generatorOnBoundary[ i-nContractibleCycles ] = true;
         }
         else
         {
            K( i ) = -mesh.defect( basisCycles[i] );
            generatorOnBoundary[ i-nContractibleCycles ] = false;
         }
      }
      
      // specify change in angle for each directional constraint
      for( unsigned int i = 0; i < directionalHolonomies.size(); i++ )
      {
         K( i+nBasisCycles ) = -directionalHolonomies[i];
      }
      
      // setup the right hand side using the Riemannian holonomy, ignoring
      // singularities for now; also make sure rhsChanged() is initially true
      for( int i = 0; i < b.length(); i++ )
      {
         b( i ) = K( i );
      }
   }
   
   void Connection :: setupRHS( void )
   // add 2*pi*k to the right hand side, where k is the vector of singularity/generator indices
   {
      double indexSum = 0;
   
      // iterate over vertices
      for( VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++ )
      {
         int i = vertex2row[ v->index ]; // get the row index of the current vertex
   
         b(i) = K(i) + 2*M_PI*v->k;
   
         indexSum += v->k;
      }
   
      // iterate over generators
      int nContractibleCycles = nBasisCycles - static_cast<int>(mesh.generatorIndices.size());
      for( int i = 0; i < (int)mesh.generatorIndices.size(); i++ )
      {
         int j = nContractibleCycles + i;
         b(j) = K(j) + 2.*M_PI*mesh.generatorIndices[i];
   
         if( generatorOnBoundary[i] )
         {
            indexSum += mesh.generatorIndices[i];
         }
      }
   
      // display a warning if the sum of singular indices does not
      // add up to the Euler characteristic
      if( abs( indexSum-mesh.eulerCharacteristic() ) > 1e-7 )
      {
         cerr << endl;
         cerr << "  *************************************************************" << endl;
         cerr << "    Warning: indices do not add up to Euler characteristic!" << endl;
         cerr << "             (solution may have unwanted singularities)" << endl;
         cerr << endl;
         cerr << "             Euler characteristic: " << mesh.eulerCharacteristic() << endl;
         cerr << "                   sum of indices: " << indexSum << endl;
         cerr << "  *************************************************************" << endl;
         cerr << endl;
      }
   }
   
   void Connection :: resetRHS( void )
   {
      // make a copy of the right hand side used in the most recent solve
      for( int i = 0; i < b.length(); i++ )
      {
         b(i) = K(i);
      }
   }
   
   bool Connection :: update( void )
   {
      Dense x( common ), y( common );
   
      // specify singular values in right hand side b
      setupRHS();
   
      // solve y = R'\(E'*b)
      y = SuiteSparseQR_solve (SPQR_RTX_EQUALS_ETB, QR, *b, common) ;
   
      // compute w = Q*y
      x = SuiteSparseQR_qmult (SPQR_QX, QR, *y, common) ;
   
      applyCotanWeights( x );
   
      for( EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); e++ )
      {
         e->theta = x( e->index );
      }
   
      // restore original right hand side
      resetRHS();
   
      return true;
   }
   
   void Connection :: applyCotanWeights( Sparse& A )
   {
      for( Sparse::iterator e = A.begin(); e != A.end(); e++ )
      {
         int row = e->first.second;
         double& val( e->second.first );
         double s = mesh.edge( row )->star;
   
         val *= sqrt( s );
      }
   }
   
   void Connection :: applyCotanWeights( cm::Dense& x )
   {
      for( EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); e++ )
      {
         x( e->index ) *= sqrt( e->star );
      }
   }
   
   void Connection :: buildCycleMatrix( Sparse& A, vector<Mesh::Cycle>& cycles ) const
   {
      for( unsigned int l = 0; l < cycles.size(); l++ )
      {
         for( Mesh::Cycle::iterator h  = cycles[l].begin();
                                    h != cycles[l].end();
                                    h ++ )
         {
            int k = (*h)->edge->index;
            int i = (*h)->from->index;
            int j = (*h)->flip->from->index;
   
            if( i > j ) A( k, l ) = -1.;
            else        A( k, l ) =  1.;
         }
      }
   }
   
   void Connection :: append1RingBases( vector<Mesh::Cycle>& cycles )
   {
      // contractible bases
      vertex2row.resize( mesh.vertices.size() );
      for( VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++ )
      {
         if( v->onBoundary() )
         {
            vertex2row[ v->index ] = -1;
            continue;
         }
   
         Mesh::Cycle c;
         HalfEdgeIter he = v->out;
         do
         {
            c.push_back( he );
            he = he->flip->next;
         }
         while( he != v->out );
   
         vertex2row[ v->index ] = static_cast<int>(cycles.size());
         cycles.push_back( c );
      }
   }
}

