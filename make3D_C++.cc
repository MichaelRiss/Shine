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
#include <sstream>
#include <vector>
#include <math.h>

using namespace std;

// This is a port of make3d of the fortran spirit package.

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

#define LW2RHO 2.35482

class entry{
public:
  int residx;
  int residy;
  int residz;
  string atmx;
  string atmy;
  string atmz;
  float shx;
  float shy;
  float shz;
  float lwx;
  float lwy;
  float lwz;
  float vol;
};

int main( int argc, char** argv ){
  if( argc != 3 ){
    cout << "Usage: make3D_C++ <peaklistfile> <spectrafile>" << endl;
    exit(1);
  }

  ifstream in( argv[1] );
  if( !in.good() ){
    cout << "Can't open " << argv[1] << ". Exiting." << endl;
    exit(1);
  }

  int nrmix;
  float mixt;
  unsigned int resx, resy, resz;
  float swx, swy, swz;

  in >> nrmix >> mixt;
  in >> resx >> resy >> resz;
  in >> swx >> swy >> swz;

  float rsx, rsy, rsz;

  rsx = (float)resx / swx;
  rsy = (float)resy / swy;
  rsz = (float)resz / swz;

  string sshx, sshy, sshz;

  vector<entry> list;

  char buffer[1024];

  in.getline( buffer, 1024 ); // finish current line

  while( true ){
    entry e;
    in.getline( buffer, 1024 );
    if( in.eof() ){
      break;
    }
    stringstream pbuf( buffer );
    pbuf >> e.residx >> e.atmx >> e.residy >> e.atmy >> e.residz >> e.atmz 
	 >> sshx >> sshy >> sshz >> e.lwx >> e.lwy >> e.lwz;
    if( sshx == "nan" ){
      e.shx = NAN;
    }
    else{
      e.shx = atof( sshx.c_str() ) - 1.0;
    }
    if( sshy == "nan" ){
      e.shy = NAN;
    }
    else{
      e.shy = atof( sshy.c_str() ) - 1.0;
    }
    if( sshz == "nan" ){
      e.shz = NAN;
    }
    else{
      e.shz = atof( sshz.c_str() ) - 1.0;
    }
    e.lwx = e.lwx * rsx;
    e.lwy = e.lwy * rsy;
    e.lwz = e.lwz * rsz;
    in.getline( buffer, 1024 );
    if( in.eof() ){
      break;
    }
    stringstream pbuf2( buffer );
    pbuf2 >> e.vol;
    e.vol = e.vol * ( 1000000.0 * 1.0 / ( e.lwx * e.lwy * e.lwz ) ) * 10000.0;
    list.push_back( e );
  }

  float* matrix = new float[resz * resy * resx];

  for( unsigned int z = 0; z < resz; z++ ){
    for( unsigned int y = 0; y < resy; y++ ){
      for( unsigned int x = 0; x < resx; x++ ){
	matrix[x + y*resx + z*resy*resx] = 0.0;
      }
    }
  }

  for( vector<entry>::iterator it = list.begin();
       it != list.end();
       it++ ){
    if( isnan( it->shx ) || isnan( it->shy ) || isnan( it->shz ) ){
      continue;
    }
    float lwzmax = 3.6 * it->lwz + 1.0;
    int z1 = (int) ( it->shz - lwzmax );
    int z2 = (int) ( it->shz + lwzmax );
    float lwymax = 3.6 * it->lwy + 1.0;
    int y1 = (int) ( it->shy - lwymax );
    int y2 = (int) ( it->shy + lwymax );
    float lwxmax = 3.6 * it->lwx + 1.0;
    int x1 = (int) ( it->shx - lwxmax );
    int x2 = (int) ( it->shx + lwxmax );
    float rhoz = it->lwz / LW2RHO;
    float rhoy = it->lwy / LW2RHO;
    float rhox = it->lwx / LW2RHO;
    for( unsigned int z = (z1+resz) % resz; 
	 z != (z2+1) % resz; z = (z+1) % resz ){
      for( unsigned int y = (y1+resy) % resy; 
	   y != (y2+1) % resy; y = (y+1) % resy ){
	for( unsigned int x = (x1+resx) % resx; 
	     x != (x2+1) %resx; x = (x+1) % resx ){
	  float dz = fabsf( (float) z - it->shz );
	  float dzf = fabsf( (float) resz - dz );
	  float dy = fabsf( (float) y - it->shy );
	  float dyf = fabsf( (float) resy - dy );
	  float dx = fabsf( (float) x - it->shx );
	  float dxf = fabsf( (float) resx - dx );
	  float add = it->vol *
	    ( expf( -( dz * dz ) / ( 2.0 * rhoz*rhoz ) ) +
	      expf( -( dzf * dzf ) / ( 2.0 * rhoz*rhoz ) ) ) *
	    ( expf( -( dy * dy ) / ( 2.0 * rhoy*rhoy ) ) +
	      expf( -( dyf * dyf ) / ( 2.0 * rhoy*rhoy ) ) ) *
	    ( expf( -( dx * dx ) / ( 2.0 * rhox*rhox ) ) +
	      expf( -( dxf * dxf ) / ( 2.0 * rhox*rhox ) ) );
	  matrix[x + y*resx + z*resy*resx] += add;
	}
      }
    }
  }

  ofstream out( argv[2] );
  if( !out.good() ){
    cout << "Cannot open " << argv[2] << ". Exiting." << endl;
    exit(1);
  }
  out.write( (char*) matrix, resz * resy * resx * sizeof(float) );
  out.close();

  delete matrix;

  return 0;
}
