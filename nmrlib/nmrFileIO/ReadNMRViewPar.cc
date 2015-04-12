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
#include "../include/nmrFileIO.h"

using namespace std;

nmrviewpar* readNMRViewPar( const char* filename, const char* name ){
  nmrviewpar tmp; // temporaray structure with reversed vectors for parsing
  strncpy( tmp.name, name, 4 );
  ifstream infile( filename );
  if( !infile ){
    cout << "Can't open " << filename << endl;
    exit(1);
  }
  string dummy;
  string keyword;
  // consume header
  infile >> keyword >> dummy >> dummy;
  if( keyword != "header" )
    cout << "header not found!" << endl;
  // read #dimensions
  infile >> keyword >> tmp.dim;
  if( keyword != "dim" )
    cout << "dim not found!" << endl;
  // read in resolution for each dimension
  for( int i = 0; i < tmp.dim; i++ ){
    int tmpres;
    infile >> tmpres;
    tmp.res.push_back( tmpres );
  }
  // consume submatrix information for each dimension
  for( int i = 0; i < tmp.dim; i++ ){
    infile >> dummy;
  }
  // consume posneg parameter
  infile >> keyword >> dummy;
  if( keyword != "posneg" )
    cout << "posneg not found!" << endl;
  // consume lvl parameter
  infile >> keyword >> dummy;
  if( keyword != "lvl" )
    cout << "lvl not found!" << endl;
  // consume scale parameter
  infile >> keyword >> dummy;
  if( keyword != "scale" )
    cout << "scale not found!" << endl;
  // consume rdims parameter
  infile >> keyword >> dummy;
  if( keyword != "rdims" )
    cout << "rdims not found!" << endl;
  // consume datatype parameter
  infile >> keyword >> dummy;
  if( keyword != "datatype" )
    cout << "datatype not found!" << endl;
  // consume byteorder parameter;
  infile >> keyword >> dummy;
  if( keyword != "byteorder" )
    cout << "byteorder not found!" << endl;
  for( int i = 0; i < tmp.dim; i++ ){
    float tmpfloat;
    string tmpstring;
    // read in sw
    infile >> keyword >> dummy >> tmpfloat;
    if( keyword != "sw" )
      cout << "sw not found!" << endl;
    tmp.sw.push_back( tmpfloat );
    // read in sf
    infile >> keyword >> dummy >> tmpfloat;
    if( keyword != "sf" )
      cout << "sf not found!" << endl;
    tmp.sf.push_back( tmpfloat );
    // read in ref
    infile >> keyword >> dummy >> tmpfloat >> dummy;
    if( keyword != "ref" )
      cout << "ref not found!" << endl;
    tmp.ref.push_back( tmpfloat );
    // read in label
    infile >> keyword >> dummy >> tmpstring;
    if( keyword != "label" )
      cout << "label not found!" << endl;
    tmp.label.push_back( tmpstring );
  }
  
  // invert vectors
  nmrviewpar* result = new nmrviewpar();
  strncpy( result->name, tmp.name, 4 );
  result->dim = tmp.dim;
  for( int i = 0; i < tmp.dim; i++ ){
    result->res.push_back( tmp.res.back() );
    tmp.res.pop_back();
    result->sw.push_back( tmp.sw.back() );
    tmp.sw.pop_back();
    result->sf.push_back( tmp.sf.back() );
    tmp.sf.pop_back();
    result->ref.push_back( tmp.ref.back() );
    tmp.ref.pop_back();
    result->label.push_back( tmp.label.back() );
    tmp.label.pop_back();
  }
  infile.close();

  return result;
}
