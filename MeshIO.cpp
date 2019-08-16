//============================================================
// MeshIO.cpp
// Keenan Crane
// 

#include "MeshIO.h"
#include "HalfEdge.h"
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <cmath>
#include <set>

using namespace std;

namespace tcods
{
   class Index
   {
      public:
         Index( void )
         {}

         Index( int p, int t, int n )
         : position( p ), texcoord( t ), normal( n )
         {}

         int position;
         int texcoord;
         int normal;
   };

   class MeshData
   {
      public:
         vector<Vector> positions;
         vector<Vector> texcoords;
         vector<Vector> normals;
         vector< vector< Index > > indices;
   };

   void MeshIO :: readOBJ( istream& in, Mesh& mesh )
   {
      MeshData data;

      readMeshData( in, data );
      buildMesh( data, mesh );
   }

   void MeshIO :: writeOBJ( ostream& out, const Mesh& mesh )
   {
      int currentIndex = 1;
      map< VertexCIter, int > vertexIndex;
      const vector<Vertex>& vertices( mesh.vertices );
      const vector<Face>& faces( mesh.faces );

      for( VertexCIter i = vertices.begin(); i != vertices.end(); i++ )
      {
         out << "v " << i->position[0] << " "
                     << i->position[1] << " "
                     << i->position[2] << endl;

         vertexIndex[ i ] = currentIndex;
         currentIndex++;
      }

      for( FaceCIter f = faces.begin(); f != faces.end(); f++ )
      {
         HalfEdgeIter he = f->he;

         for( int j = 0; j < 3; j++ )
         {
            out << "vt " << he->texcoord.x << " " << he->texcoord.y << endl;
            he = he->next;
         }
      }

      for( FaceCIter i = faces.begin(); i != faces.end(); i++ )
      {
         HalfEdgeIter he = i->he;

         // don't write boundary faces
         if( he->onBoundary )
         {
            continue;
         }

         out << "f ";

         int j = 0;
         do
         {
            out << vertexIndex[ he->from ] << "/" << 1+(i->index*3+j) << " ";
            he = he->next;
            j++;
         }
         while( he != i->he );

         out << endl;
      }
   }

   void MeshIO :: writeOBJX( ostream& out, const Mesh& mesh )
   {
      out.precision( 10 );

      int currentIndex = 1;
      map< VertexCIter, int > vertexIndex;
      const vector<Vertex>& vertices( mesh.vertices );
      const vector<Face>& faces( mesh.faces );

      for( VertexCIter i = vertices.begin(); i != vertices.end(); i++ )
      {
         out << "v " << i->position[0] << " "
                     << i->position[1] << " "
                     << i->position[2] << endl;

         vertexIndex[ i ] = currentIndex;
         currentIndex++;
      }

      for( VertexCIter i = vertices.begin(); i != vertices.end(); i++ )
      {
         Vector u( 0., 0., 0. );
         HalfEdgeIter he = i->out;
         do
         {
            FaceCIter f = he->face;
            if( !f->isBoundary() )
            {
               double alpha = f->alpha;
               Vector w( cos(alpha), sin(alpha), 0. );
               u += f->toGlobal( w );
            }
            he = he->flip->next;
         }
         while( he != i->out );
         u.normalize();
         
         out << "vf " << u.x << " " << u.y << " " << u.z << endl;
      }

      for( FaceCIter i = faces.begin(); i != faces.end(); i++ )
      {
         HalfEdgeIter he = i->he;

         // don't write boundary faces
         if( he->onBoundary )
         {
            continue;
         }

         out << "f ";

         do
         {
            out << vertexIndex[ he->from ] << " ";
            he = he->next;
         }
         while( he != i->he );

         out << endl;
      }
   }

   void MeshIO :: writeEOBJ( ostream& out, const Mesh& mesh )
   {
      out.precision( 10 );

      int currentIndex = 1;
      map< VertexCIter, int > vertexIndex;
      map< FaceCIter, int > faceIndex;
      const vector<Vertex>& vertices( mesh.vertices );
      const vector<Face>& faces( mesh.faces );

      for( VertexCIter i = vertices.begin(); i != vertices.end(); i++ )
      {
         out << "v " << i->position[0] << " "
                     << i->position[1] << " "
                     << i->position[2] << endl;

         vertexIndex[ i ] = currentIndex;
         currentIndex++;
      }

      currentIndex = 1;
      for( FaceCIter i = faces.begin(); i != faces.end(); i++ )
      {
         HalfEdgeIter he = i->he;

         // don't write boundary faces
         if( he->onBoundary )
         {
            continue;
         }

         faceIndex[ i ] = currentIndex;
         currentIndex++;

         out << "f ";

         do
         {
            out << vertexIndex[ he->from ] << " ";
            he = he->next;
         }
         while( he != i->he );

         out << endl;
      }

      for( FaceCIter i = faces.begin(); i != faces.end(); i++ )
      {
         HalfEdgeIter he = i->he;

         // don't write vectors on boundary faces
         if( he->onBoundary )
         {
            continue;
         }

         out << "# attrs f " << faceIndex[ i ] << " ";

         double alpha = i->alpha;
         Vector w( cos(alpha), sin(alpha), 0. );
         Vector u = i->toGlobal( w );

         out << u.x << " " << u.y << " " << u.z;

         out << endl;
      }
   }

   int modulo( int a, int b )
   {
      a += b*(1+abs(a/b));
      return a % b;
   }

   void MeshIO :: writeJVX( ostream& out, const Mesh& mesh, const vector<int>& singularityIndices )
   {
      out.precision( 10 );

      out << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\" standalone=\"no\"?>\n";
      out << "<!DOCTYPE jvx-model SYSTEM \"http://www.javaview.de/rsrc/jvx.dtd\">\n";
      out << "<jvx-model>\n";
      out << "\t<geometries>\n";
      out << "\t\t<geometry name=\"trivialconnections_output\">\n";
      out << "\t\t\t<pointSet dim=\"3\">\n";
      out << "\t\t\t\t<points num=\"" << (int)mesh.vertices.size() << "\">\n";
      for( vector<VertexIter>::const_iterator i  = mesh.index2vertex.begin();
                                              i != mesh.index2vertex.end();
                                              i ++ )
      {
         Vector p = (*i)->position;
         out << "\t\t\t\t\t<p>" << p.x << " " << p.y << " " << p.z << "</p>\n";
      }
      out << "\t\t\t\t</points>\n";
      out << "\t\t\t</pointSet>\n";
      out << "\t\t\t<faceSet>\n";
      out << "\t\t\t\t<faces num=\"" << (int)mesh.faces.size() << "\">\n";
      for( vector<FaceIter>::const_iterator i  = mesh.index2face.begin();
                                            i != mesh.index2face.end();
                                            i ++ )
      {
         FaceCIter f = *i;
         int ii = f->he->from->index;
         int j = f->he->next->from->index;
         int k = f->he->next->next->from->index;
         out << "\t\t\t\t\t<f>" << ii << " " << j << " " << k << "</f>\n";
      }
      out << "\t\t\t\t</faces>\n";
      out << "\t\t\t</faceSet>\n";
      out << "\t\t\t<vectorField name=\"First field\" base=\"element\">\n";
      out << "\t\t\t\t<vectors num=\"" << (int)mesh.faces.size() << "\">\n";
      for( vector<FaceIter>::const_iterator i  = mesh.index2face.begin();
                                            i != mesh.index2face.end();
                                            i ++ )
      {
         FaceCIter f = *i;
         double alpha = f->alpha;
         Vector e1, e2; f->frame( e1, e2 );
         Vector v = cos(alpha)*e1 + sin(alpha)*e2;
         out << "\t\t\t\t\t<v>" << v.x << " " << v.y << " " << v.z << "</v>\n";
      }
      out << "\t\t\t\t</vectors>\n";
      out << "\t\t\t</vectorField>\n";
      out << "\t\t\t<vectorField name=\"Second field\" base=\"element\">\n";
      out << "\t\t\t\t<vectors num=\"" << (int)mesh.faces.size() << "\">\n";
      for( vector<FaceIter>::const_iterator i  = mesh.index2face.begin();
                                            i != mesh.index2face.end();
                                            i ++ )
      {
         FaceCIter f = *i;
         double alpha = f->alpha + M_PI/2.;
         Vector e1, e2; f->frame( e1, e2 );
         Vector v = cos(alpha)*e1 + sin(alpha)*e2;
         out << "\t\t\t\t\t<v>" << v.x << " " << v.y << " " << v.z << "</v>\n";
      }
      out << "\t\t\t\t</vectors>\n";
      out << "\t\t\t</vectorField>\n";
      out << "\t\t\t<matchingSet>\n";
      out << "\t\t\t\t<matchings num=\"" << (int)mesh.edges.size() << "\">\n";
      for( vector<EdgeIter>::const_iterator i  = mesh.index2edge.begin();
                                            i != mesh.index2edge.end();
                                            i ++ )
      {
         EdgeCIter e = *i;
         HalfEdgeCIter he = e->he;
         if( he->from->index > he->flip->from->index ) he = he->flip;
         double theta = e->theta;
         double alphaI = he->face->alpha;
         double alphaJ = he->flip->face->alpha;
         double alphaIp = Mesh::parallelTransport( alphaI, he ) - theta;
         int m = (int) round(( alphaJ - alphaIp ) / ( M_PI / 2. ));
         m = modulo( m, 4 );
         int ii = he->face->index;
         int j = he->flip->face->index;
         out << "\t\t\t\t\t<m>" << ii << " " << j << " " << m << "</m>\n";
      }
      out << "\t\t\t\t</matchings>\n";
      out << "\t\t\t</matchingSet>\n";
      out << "\t\t</geometry>\n";
      out << "\t</geometries>\n";

      out << "<!-- singular vertices -->\n";
      for( vector<int>::const_iterator i = singularityIndices.begin(); i != singularityIndices.end(); i++ )
      {
         out << "<!-- " << *i << " -->\n";
      }

      out << "</jvx-model>\n";
   }

   void MeshIO :: readMeshData( istream& in, MeshData& data )
   {
      string line;

      while( getline( in, line ))
      {
         stringstream ss( line );
         string token;

         ss >> token;

         if( token == "v"  ) readPosition( ss, data );
         if( token == "vt" ) readTexCoord( ss, data );
         if( token == "vn" ) readNormal  ( ss, data );
         if( token == "f"  ) readFace    ( ss, data );
      }
   }

   void MeshIO :: preallocateMeshElements( const MeshData& data, Mesh& mesh )
   {
      // count the number of edges
      set< pair<int,int> > edges;
      for( vector< vector< Index > >::const_iterator f  = data.indices.begin();
                                                     f != data.indices.end();
                                                     f ++ )
      {
         for( unsigned int I = 0; I < f->size(); I++ )
         {
            int J = (I+1) % f->size();
            int i = (*f)[I].position;
            int j = (*f)[J].position;

            if( i > j ) swap( i, j );

            edges.insert( pair<int,int>( i, j ));
         }
      }

      int nV = static_cast<int>(data.positions.size());
      int nE = static_cast<int>(edges.size());
      int nF = static_cast<int>(data.indices.size());
      int nHE = 2*nE;
      int chi = nV - nE + nF;
      int nB = 2 - chi; // (conservative approximation of number of boundary cycles)

      mesh.vertices.reserve( nV );
      mesh.edges.reserve( nE );
      mesh.faces.reserve( nF + nB );
      mesh.halfedges.reserve( nHE );
   }

   void MeshIO :: buildMesh( const MeshData& data, Mesh& mesh )
   {
      map< pair< int, int >, HalfEdgeIter > existingHalfEdges;
      map< int, VertexIter > indexToVertex;
      map< HalfEdgeIter, bool > hasFlipEdge;

      preallocateMeshElements( data, mesh );

      // allocate a vertex for each position in the data and construct
      // a map from vertex indices to vertex pointers
      for( unsigned int i = 0; i < data.positions.size(); i++ )
      {
         VertexIter newVertex = mesh.vertices.insert( mesh.vertices.end(), Vertex() );
         newVertex->position = data.positions[ i ];
         indexToVertex[ i ] = newVertex;
      }

      // insert each face into the mesh
      for( vector< vector< Index > >::const_iterator f  = data.indices.begin();
                                                     f != data.indices.end();
                                                     f ++ )
      {
         // create a new face
         FaceIter newFace = mesh.faces.insert( mesh.faces.end(), Face());

         // create a new half edge for each edge of the current face
         int N = static_cast<int>(f->size());
         vector< HalfEdgeIter > hes( N );
         for( int i = 0; i < N; i++ )
         {
            hes[ i ] = mesh.halfedges.insert( mesh.halfedges.end(), HalfEdge());
         }

         // initialize these new halfedges
         for( int i = 0; i < N; i++ )
         {
            // the current halfedge goes from vertex a to vertex b
            int a = (*f)[     i     ].position;
            int b = (*f)[ (i+1) % N ].position;

            // set current halfedge's attributes
            hes[ i ]->next = hes[ (i+1) % N ];
            hes[ i ]->from = indexToVertex[ a ];
            int t = (*f)[i].texcoord;
            int n = (*f)[i].normal;
            if( t >= 0 ) hes[ i ]->texcoord = data.texcoords[ t ];
            else         hes[ i ]->texcoord = Vector( 0., 0., 0. );
            if( n >= 0 ) hes[ i ]->normal   = data.normals  [ n   ];
            else         hes[ i ]->normal   = Vector( 0., 0., 0. );
            hes[ i ]->onBoundary = false;

            // keep track of which halfedges have flip edges defined (for detecting boundaries)
            hasFlipEdge[ hes[ i ]] = false;

            // point vertex a at the current halfedge
            indexToVertex[ a ]->out = hes[ i ];

            // point the new face and this half edge to each-other
            hes[ i ]->face = newFace;
            newFace->he = hes[ i ];

            // if we've created an edge between a and b in the past, it is the
            // flip edge of the current halfedge
            if( existingHalfEdges.find( pair<int,int>( a, b )) != existingHalfEdges.end())
            {
               hes[ i ]->flip = existingHalfEdges[ pair<int,int>( a, b ) ];
               hes[ i ]->flip->flip = hes[ i ];
               hes[ i ]->edge = hes[ i ]->flip->edge;
               hasFlipEdge[ hes[ i ]] = true;
               hasFlipEdge[ hes[ i ]->flip ] = true;
            }
            else // otherwise, create an edge connected to the current halfedge
            {
               hes[ i ]->edge = mesh.edges.insert( mesh.edges.end(), Edge());
               hes[ i ]->edge->he = hes[i];
            }

            // record the fact that we've created a halfedge from a to b
            existingHalfEdges[ pair<int,int>( a, b ) ] = hes[ i ];
            existingHalfEdges[ pair<int,int>( b, a ) ] = hes[ i ];
         }
      }

      // insert extra faces for each boundary cycle
      for( HalfEdgeIter currentHE  = mesh.halfedges.begin();
                        currentHE != mesh.halfedges.end();
                        currentHE ++ )
      {
         // if we find a halfedge with no flip edge defined, create
         // a new face and link it to the corresponding boundary cycle

         if( !hasFlipEdge[ currentHE ] )
         {
            // create a new face
            FaceIter newFace = mesh.faces.insert( mesh.faces.end(), Face());

            // walk along this boundary cycle
            vector<HalfEdgeIter> boundaryCycle;
            HalfEdgeIter he = currentHE;
            do
            {
               // create a new halfedge on the boundary face
               HalfEdgeIter newHE = mesh.halfedges.insert( mesh.halfedges.end(), HalfEdge());

               // mark only the halfedge on the boundary face as being on the boundary
               newHE->onBoundary = true;

               // link the current halfedge in the cycle to its new flip edge
               he->flip = newHE;

               // grab the next halfedge along the boundary by finding
               // the next halfedge around the current vertex that doesn't
               // have a flip edge defined
               HalfEdgeIter nextHE = he->next;
               while( hasFlipEdge[ nextHE ] )
               {
                  nextHE = nextHE->flip->next;
               }

               // set attributes for the flip edge (we'll set ->next below)
               newHE->flip = he;
               newHE->from = nextHE->from;
               newHE->edge = he->edge;
               newHE->face = newFace;
               newHE->normal = nextHE->normal;
               newHE->texcoord = nextHE->texcoord;

               // point the new face to this half edge
               newFace->he = newHE;

               // keep track of all the new halfedges in the boundary cycle
               boundaryCycle.push_back( newHE );

               // continue to walk along the cycle
               he = nextHE;

            } while( he != currentHE );

            // link together the cycle of boundary halfedges
            unsigned int N = static_cast<unsigned int>(boundaryCycle.size());
            for( unsigned int i = 0; i < N; i++ )
            {
               boundaryCycle[ i ]->next = boundaryCycle[ (i+N-1)%N ];
               hasFlipEdge[ boundaryCycle[i] ] = true;
               hasFlipEdge[ boundaryCycle[i]->flip ] = true;
            }
         }
      }
   }

   void MeshIO :: readPosition( stringstream& ss, MeshData& data )
   {
      double x, y, z;

      ss >> x >> y >> z;

      data.positions.push_back( Vector( x, y, z ));
   }

   void MeshIO :: readTexCoord( stringstream& ss, MeshData& data )
   {
      double u, v;

      ss >> u >> v;

      data.texcoords.push_back( Vector( u, v, 0. ));
   }

   void MeshIO :: readNormal( stringstream& ss, MeshData& data )
   {
      double x, y, z;

      ss >> x >> y >> z;

      data.normals.push_back( Vector( x, y, z ));
   }

   void MeshIO :: readFace( stringstream& ss, MeshData &data )
   {
      vector<Index> faceIndices;
      string token;

      while( ss >> token )
      {
         faceIndices.push_back( parseFaceIndex( token ));
      }

      data.indices.push_back( faceIndices );
   }

   Index MeshIO :: parseFaceIndex( const string& token )
   {
      // parse indices of the form
      //
      // p/[t]/[n]
      //
      // where p is an index into positions, t is an index into
      // texcoords, n is an index into normals, and [.] indicates
      // that an index is optional
      
      stringstream in( token );
      string indexstring;
      int indices[3] = { -1, -1, -1 };
      int i = 0;
   
      while( getline( in, indexstring, '/' ))
      {
         stringstream ss( indexstring );
         ss >> indices[i++];
      }
   
      // decrement since indices in OBJ files are 1-based
      return Index( indices[0]-1,
                    indices[1]-1,
                    indices[2]-1 );
   }
}

