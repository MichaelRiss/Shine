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
#include <iomanip>
#include <fstream>
#include <sstream>
#include "../include/nmrFileIO.h"
#include "../include/amino.h"

vector<aminoshift> readShifts( const char* filename ){
  ifstream shiftfile( filename );
  vector<aminoshift> list;
  if( !shiftfile ){
    cerr << "Can't open Shiftfile " << filename << endl;
    shiftfile.close();
    throw FileNotFoundException( filename );
  }
  
  //vector<aminoshift>* list = new vector<aminoshift>();

  int curseq = -1;
  aminoshift tmp;
  bool tmpvalid = false;

  char line[1024];
  while( true ){
    shiftfile.getline( line, 1024 );
    if( shiftfile.eof() ) break;
    
    stringstream linestr( line );
    string amino_nr;
    string atom;
    string HshiftS;
    string CNshiftS;
    char amino;
    int seqnr;
    float Hshift = 0.0;
    float CNshift = 0.0;
    // save parse and check for invalid lines
    linestr >> amino_nr >> atom >> HshiftS >> CNshiftS;

    if( amino_nr == "END" || 
	atom == "-" || 
	HshiftS == "-" || 
	CNshiftS == "-" ){
      continue;
    }
    
    // real parse
    linestr.seekg( 0, ios_base::beg );

    linestr.clear();

    linestr >> amino_nr >> atom >> Hshift >> CNshift;

    stringstream amino_nrS( amino_nr );
    amino_nrS >> amino >> seqnr;
    
    if( seqnr != curseq ){
      if( tmpvalid ){
	list.push_back( tmp );
	tmpvalid = false;
      }
      curseq = seqnr;
      tmp.seqnr = seqnr;
      tmp.amino = amino;
      tmp.CProt.clear();
      tmp.NProt.clear();
      tmp.Prot.clear();
      tmpvalid = true;
    }
    
    shiftpair tmp2;
    tmp2.hetero = CNshift;
    tmp2.proton = Hshift;
    tmp2.pname = atom;

    const char* heavy = ac_oneletter_longlist2heavy( amino, atom.c_str() );
    if( heavy[0] == 'C' ){
      tmp2.htype = 'C';
      tmp.CProt[atom] = tmp2;
    }
    else{
      tmp2.htype = 'N';
      tmp.NProt[atom] = tmp2;
    }
    tmp.Prot[atom] = tmp2;
    
  }
  if( tmpvalid ){
    list.push_back( tmp );
    tmpvalid = false;
  }

  shiftfile.close();
  return list;
}
