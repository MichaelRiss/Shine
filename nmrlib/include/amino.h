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

#include <vector>
#include <unordered_map>

using namespace std;
using namespace __gnu_cxx;

enum amintype { ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, 
		LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL };

#define ACIDIC 1
#define ACYCLIC 2
#define ALIPHATIC 4
#define AROMATIC 8
#define BASIC 16
#define BURIED 32
#define CHARGED 64
#define CYCLIC 128
#define HYDROPHOBIC 256
#define LARGE 512
#define MEDIUM 1024
#define NEGATIVE 2048
#define NEUTRAL 4096
#define POLAR 8192
#define POSITIVE 16384
#define SMALL 32768
#define SURFACE 65536


typedef unsigned int aminoflags;

extern unordered_map<string, aminoflags> aflags;

extern void InitAminoFlags();

extern unordered_map<string, float> centroidsize;

extern void InitCentroidSize();

// Conversion functions

// 1   one-letter code
// 2   three-letter code
// 3   xplor format
// 4   iupac format
// 5   longlist format
// 6   general format = old-PDB u.A
// 7   heavy atom
// 8   xplor pseudo-atom
// 9   Aqua pseudo-atom

extern const char* ac_oneletter_longlist2heavy( const char oneletter, 
						const char* longlist );

extern const char* ac_oneletter_xplor2heavy( const char oneletter, 
					     const char* xplor );

extern const char* ac_oneletter2threeletter( const char oneletter );

extern const char ac_threeletter2oneletter( const char* threeletter );

extern const vector<const char*> 
ac_threeletter_xplor2longlist( const char* threeletter,
			       const char* xplor );

extern const char* ac_threeletter_xplor2heavy( const char* threeletter,
					       const char* xplor );
extern const char* ac_threeletter_iupac2xplor( const char* threeletter,
					       const char* iupac );

