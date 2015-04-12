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

#include <iostream>
#include <fstream>
#include "../include/nmrFileIO.h"

using namespace std;

vector<peakentry>* readPeaklist( const char* filename ){
  ifstream input( filename );
  if( !input.good() ){
    cerr << "Can't open Peaklist " << filename << endl;
    input.close();
    return NULL;
  }

  vector<peakentry>* list = new vector<peakentry>();
  
  while( !input.eof() ){
    peakentry tmp;
    input >> tmp.x >> tmp.y >> tmp.z >> tmp.intensity
	  >> tmp.residshort1 >> tmp.sequencenumber1 >> tmp.atom1
	  >> tmp.residshort2 >> tmp.sequencenumber2 >> tmp.atom2;
    if( !input.fail() ){
      list->push_back( tmp );
    }
  }

  input.close();

  return list;
}

vector<peakentry_light>* readPeaklist_light( const char* filename ){
  ifstream input( filename );
  if( !input.good() ){
    cerr << "Can't open Peaklist " << filename << endl;
    input.close();
    return NULL;
  }

  vector<peakentry_light>* list = new vector<peakentry_light>();

  while( !input.eof() ){
    peakentry_light tmp;
    input >> tmp.x >> tmp.y >> tmp.z >> tmp.intensity;
    if( !input.fail() ){
      list->push_back( tmp );
    }
  }

  input.close();

  return list;
}
