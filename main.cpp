#include "Problem.h"
#include <iostream>
#include <fstream>

using namespace std;
using namespace tcods;

int main( int argc, char** argv )
{
   if( argc != 2 )
   {
      cerr << "usage: " << argv[0] << " problem.txt" << endl;
      return 1;
   }

   ifstream in( argv[1] );
   if( !in.is_open() )
   {
      cerr << "Error: could not open problem file " << argv[1] << endl;
      return 1;
   }

   Problem problem( in );
   problem.solve();

   return 0;
}

