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

#include <stdio.h>
#include <iostream>
#include <fstream>
#include "../include/nmrFileIO.h"

using namespace std;

vector<betasheet> readBetaSheets( const char* filename ){
  vector<betasheet> result;
  
  ifstream input( filename );
  if( !input ){
    cout << "ReadBetaSheets: cannot open file " << filename << endl;
    exit(1);
  }

  char buffer[1024];
  while( true ){
    input.getline( buffer, 1024 );
    if( input.eof() )
      break;
    int begin = -1;
    int end = -1;
    char chainID = '\0';

    sscanf( buffer, "%i %i \"%c\"", &begin, &end, &chainID );

    if( begin == -1 || end == -1 || chainID == '\0' ){
      cout << "Parse error: begin: " << begin << ", end: " 
	   << end << ", chainID: " << chainID << endl;
      exit(1);
    }

    betasheet tmp;
    tmp.begin = begin;
    tmp.end = end;
    tmp.chainID = chainID;
    result.push_back( tmp );
  }

  input.close();
  return result;
}
