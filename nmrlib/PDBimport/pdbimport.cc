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

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <iostream>
#include "../include/pdbimport.h"

using namespace std;

vector<pdbatom>* readatoms( const char* pdbfilename ){
  FILE* pdbfile;
  char* scanline = NULL;
  size_t linelength;
  vector<pdbatom>* atoms;
  
  pdbfile = fopen( pdbfilename, "r" );
  if( pdbfile == NULL ){
    cout << "Can't open file: " << pdbfilename << endl;
    exit(1);
  }

  atoms = new vector<pdbatom>();

  while( getline( &scanline, &linelength, pdbfile ) != -1 ){
    if( strncmp( scanline, "END", 3 ) == 0 )
      break;
    if( strncmp( scanline, "REMARK", 6 ) == 0 )
      continue;
    if( strncmp( scanline, "ATOM  ", 6 ) == 0 ){
      pdbatom tmp;

      // serial
      char buffer[9];
      char* chopped;
      buffer[0] = scanline[6];
      buffer[1] = scanline[7];
      buffer[2] = scanline[8];
      buffer[3] = scanline[9];
      buffer[4] = scanline[10];
      buffer[5] = (char) NULL;
      tmp.serial = atoi( buffer );

      // atomname
      buffer[0] = scanline[12];
      buffer[1] = scanline[13];
      buffer[2] = scanline[14];
      buffer[3] = scanline[15];
      buffer[4] = (char) NULL;
      sscanf( buffer, "%as", &chopped );
      tmp.atomname = chopped;
      free( chopped );
      
      // residname
      buffer[0] = scanline[17];
      buffer[1] = scanline[18];
      buffer[2] = scanline[19];
      buffer[3] = (char) NULL;
      sscanf( buffer, "%as", &chopped );
      tmp.residname = chopped;
      free( chopped );

      // chainID
      tmp.chainID = scanline[21];

      // residnumber
      buffer[0] = scanline[22];
      buffer[1] = scanline[23];
      buffer[2] = scanline[24];
      buffer[3] = scanline[25];
      buffer[4] = (char) NULL;
      tmp.residnumber = atoi( buffer );

      // xcoord
      buffer[0] = scanline[30];
      buffer[1] = scanline[31];
      buffer[2] = scanline[32];
      buffer[3] = scanline[33];
      buffer[4] = scanline[34];
      buffer[5] = scanline[35];
      buffer[6] = scanline[36];
      buffer[7] = scanline[37];
      buffer[8] = (char) NULL;
      tmp.xcoord = (float) atof( buffer );

      // ycoord
      buffer[0] = scanline[38];
      buffer[1] = scanline[39];
      buffer[2] = scanline[40];
      buffer[3] = scanline[41];
      buffer[4] = scanline[42];
      buffer[5] = scanline[43];
      buffer[6] = scanline[44];
      buffer[7] = scanline[45];
      buffer[8] = (char) NULL;
      tmp.ycoord = (float) atof( buffer );

      // zcoord
      buffer[0] = scanline[46];
      buffer[1] = scanline[47];
      buffer[2] = scanline[48];
      buffer[3] = scanline[49];
      buffer[4] = scanline[50];
      buffer[5] = scanline[51];
      buffer[6] = scanline[52];
      buffer[7] = scanline[53];
      buffer[8] = (char)NULL;
      tmp.zcoord = (float) atof( buffer );

      // append to vector
      atoms->push_back( tmp );
    }
  }

  fclose( pdbfile );

  return atoms;
}
