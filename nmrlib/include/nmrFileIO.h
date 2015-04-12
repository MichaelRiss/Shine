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

#include <math.h>
#include <vector>
#include <unordered_map>

using namespace std;
using namespace __gnu_cxx;

typedef struct{
  int seqnr;
  char amino;
} seqitem;

extern vector<seqitem>* readSeq( const char* filename );

class shiftpair{
 public:
  float proton;
  float hetero;
  string pname;
  char htype;
  shiftpair();
};

class aminoshift{
 public:
  int seqnr;
  char amino;
  unordered_map<string, shiftpair> CProt;
  unordered_map<string, shiftpair> NProt;
  unordered_map<string, shiftpair> Prot;
  aminoshift();
};

extern vector<aminoshift> readShifts( const char* filename );

typedef struct{
  float x;
  float y;
  float z;
  float intensity;
  char residshort1;
  char residshort2;
  int sequencenumber1;
  int sequencenumber2;
  char atom1[4];
  char atom2[4];
} peakentry;

typedef struct{
  float x;
  float y;
  float z;
  float intensity;
} peakentry_light;

vector<peakentry>* readPeaklist( const char* filename );

vector<peakentry_light>* readPeaklist_light( const char* filename );

typedef struct{
  char name[4];
  int dim;
  vector<int> res;
  vector<float> sw;
  vector<float> sf;
  vector<float> ref;
  vector<string> label;
} nmrviewpar;

nmrviewpar* readNMRViewPar( const char* filename, const char* name );

typedef struct{
  int begin;
  int end;
  char chainID;
} betasheet;

vector<betasheet> readBetaSheets( const char* filename );

typedef struct{
  float H2O_line;
  float H2O_width;
  float diag_width;
} nogo;

nogo readNogo( const char* filename );

class FileNotFoundException{
 public:
  string filename;
  FileNotFoundException( string );
};
