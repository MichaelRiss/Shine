/*  Shine - NMR NOESY simulation program
    Copyright (C) 2015  Michael Riss

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. */

#include <string.h>
#include <iostream>
#include <fstream>
#include "../include/amino.h"

using namespace std;

int main( int argc, char** argv ){
  
  if( argc != 3 ){
    cout << "Usage: iupac2xplor_pdb <input.pdb> <output.pdb>" << endl;
    exit(1);
  }

  ifstream input( argv[1] );
  if( !input ){
    cout << "Can't open inputfile argv[1]." << endl;
    exit(1);
  }
  ofstream output( argv[2] );
  if( !output ){
    cout << "Can't open outputfile argv[2]." << endl;
    exit(1);
  }

  char line[1024];
  char buffer[5];
  char threeletter[4];
  while( true ){
    input.getline( line, 1024 );
    if( input.eof() )
      break;

    // search for ATOM
    if( line[0] == 'A' &&
	line[1] == 'T' &&
	line[2] == 'O' &&
	line[3] == 'M' ){
      buffer[0] = line[12];
      buffer[1] = line[13];
      buffer[2] = line[14];
      buffer[3] = line[15];
      buffer[4] = '\0';
      sscanf( line+12, "%s", buffer );
      sscanf( line+17, "%s", threeletter );
      cout << "invoke ac_iupac( " << threeletter << ", " 
	   << buffer << " );" << endl;
      const char* conv = ac_threeletter_iupac2xplor( threeletter, buffer);
      if( conv != NULL ){
	line[12] = ' ';
	line[13] = ' ';
	line[14] = ' ';
	line[15] = ' ';
	strncpy( line+12, conv , strlen( conv ) );
      }
    }
    output << line << endl;
  }

  input.close();
  output.close();
  
  return 0;
}
