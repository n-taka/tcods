//============================================================
// MeshIO.h
// 
// MeshIO contains a collection of static members for polygon
// mesh input/output via the Mesh class.
//

#ifndef TCODS_MESHIO_H
#define TCODS_MESHIO_H

#include <iosfwd>
#include <string>
#include <sstream>
#include <string>
#include <vector>

namespace tcods
{
   class Index;
   class MeshData;
   class Mesh;

   class MeshIO
   {
      public:
         // routines for reading and writing various mesh formats
         // assumes valid, open streams that point to the target file

         // read
         static void readOBJ( std::istream& in, Mesh& mesh ); // reads WavefrontOBJ

         // write
         static void  writeOBJ( std::ostream& out, const Mesh& mesh ); // writes WavefrontOBJ
         static void writeOBJX( std::ostream& out, const Mesh& mesh ); // writes WavefrontOBJ plus tangent vector per vertex
         static void writeEOBJ( std::ostream& out, const Mesh& mesh ); // writes WavefrontOBJ plus tangent vector per face

         static void  writeJVX( std::ostream& out, const Mesh& mesh, const std::vector<int>& singularityIndices );
         // writes XML containing mesh data, plus pair of orthogonal tangent vector fields, singularity
         // locations, and "matchings" between neighboring triangles which define a quad-covering of the surface
         // (can be used as input for QuadCover)

      protected:
         static void readMeshData( std::istream& in, MeshData& data );
         static void buildMesh( const MeshData& data, Mesh& mesh );
         static void readPosition( std::stringstream& ss, MeshData& data );
         static void readTexCoord( std::stringstream& ss, MeshData& data );
         static void readNormal  ( std::stringstream& ss, MeshData& data );
         static void readFace    ( std::stringstream& ss, MeshData& data );
         static Index parseFaceIndex( const std::string& token );
         static void preallocateMeshElements( const MeshData& data, Mesh& mesh );
   };
}

#endif

