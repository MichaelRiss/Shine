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

#ifndef __NMR_DATA
#define __NMR_DATA

#include <vector>
#include <string>

using namespace std;

class hetero;

class proton{
 public:
  string Name;
  double shift;
  double spec_shift1;
  double spec_shift2;
  hetero* hetatom;
  double x;
  double y;
  double z;

  double leak;
  double fold2;

  /* lw */
  double linewidth;
  double lw1;
  double lw2;
  double valvef;
  proton( string Name = "", char amino = ' ' );
};

class amino;

class hetero{
 public:
  string Name;
  double shift;
  double spec_shift2;
  double spec_shift3;
  vector<proton*> prot;
  amino* residue;
  double x;
  double y;
  double z;

  double fold2;
  double fold3;
  /* lw */
  double linewidth;
  double lw2;
  double lw3;
  double rj;
  double trans;
  double valveg;
  hetero( string Name = "", char aminoCode = ' ', amino* residue = NULL );
  void addproton( proton* h );
};

class amino{
 public:
  int residueNumber;
  char residueName;
  vector<hetero*> C13;
  vector<hetero*> N15;
  vector<proton*> H;
  amino( char residueName = '0', int residueNumber = -1 );
  void addC13( hetero* het );
  void addN15( hetero* het );
};

class protein{
 public:
  string Name;
  vector<amino*> chain;
  vector<hetero*> C13;
  vector<hetero*> N15;
  vector<proton*> H;
  vector<proton*> HN;
  vector<proton*> HC;
  protein();
};

#endif /* __NMR_DATA */
