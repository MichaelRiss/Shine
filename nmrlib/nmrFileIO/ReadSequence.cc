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

#include <fstream>
#include <iostream>
#include <vector>
#include "../include/nmrFileIO.h"

using namespace std;

vector<seqitem>* readSeq( const char* filename ){
  ifstream seqfile( filename );
  if( !seqfile.good() ){
    cerr << "Can't open sequence file " << filename << endl;
    seqfile.close();
    return NULL;
  }

  vector<seqitem>* list = new vector<seqitem>();

  int seqnr = -1;
  seqitem tmp;

  seqfile >> seqnr;
  if( !seqfile.good() ){
    cerr << "Error reading sequence file " << filename << endl;
    seqfile.close();
    return NULL;
  }
  
  for( ;; ){
    seqfile >> tmp.amino;
    if( seqfile.eof() ) break;
    tmp.seqnr = seqnr++;
    list->push_back( tmp );
  }

  seqfile.close();
  return list;
}
