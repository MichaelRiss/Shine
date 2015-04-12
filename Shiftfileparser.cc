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
#include "Shiftfileparser.h"

using namespace std;

shiftfileitem::shiftfileitem(){
}

vector<shiftfileitem> readshiftfile( const char* filename ){
  vector<shiftfileitem> atoms;

  ifstream input( filename );
  if( !input ){
    cout << "Can't open " << filename << "."<< endl;
    exit(1);
  }
  
  while( true ){
    shiftfileitem tmp;
    input >> tmp.resid >> tmp.resname >> tmp.atomname >> tmp.shift;
    if( !input.good() ){
      break;
    }
    atoms.push_back( tmp );
  }
  input.close();

  return atoms;
}
