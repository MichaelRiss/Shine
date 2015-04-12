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

#include <cmath>
#include "NMRdata.h"

#define LW_C 150.0
#define LW_N 30.0
#define LW_H 35.0
#define LW_HN 16.0

#define LEAK0 0.5
#define LEAK1 0.65

/* Originally the LEAK2 parameter has been set to 0.75 in the 
   spirit fortran program, but as the LEAK2 parameter has not been 
   declared it defaults to integer an the value of 0.75 mapped to 0.0
   Bug or feature?
   For compatibility we also set the value to 0.0, but s.o. with chemical
   knowledge should find out what the correct value for LEAK2 is.
   #define LEAK2 0.75
*/
#define LEAK2 0.0

proton::proton( string Name, char aminoCode ){
  this->Name = Name;
  shift = NAN;
  spec_shift1 = NAN;
  spec_shift2 = NAN;
  hetatom = NULL;
  x = NAN;
  y = NAN;
  z = NAN;
  leak = NAN;
  if( Name == "HA" || 
      ( ( aminoCode == 'R' ) && 
	( Name == "HB1" || Name == "HB2" || 
	  Name == "HG1" || Name == "HG2" ||
	  Name == "HD1" || Name == "HD2" ) ) ||
      ( ( aminoCode == 'N' ) && 
	( Name == "HB1" || Name == "HB2" ) ) ||
      ( ( aminoCode == 'D' ) &&
	( Name == "HB1" || Name == "HB2" ) ) ||
      ( ( aminoCode == 'C' ) &&
	( Name == "HB1" || Name == "HB2" ) ) ||
	( ( aminoCode == 'Q' ) &&
	( Name == "HB1" || Name == "HB2" ||
	  Name == "HG1" || Name == "HG2" ) ) ||
      ( ( aminoCode == 'E' ) &&
	( Name == "HB1" || Name == "HB2" ||
	  Name == "HG1" || Name == "HG2" ) ) ||
      ( ( aminoCode == 'G' ) &&
	( Name == "HA1" || Name == "HA2" ) ) ||
      ( ( aminoCode == 'H' ) &&
	( Name == "HB1" || Name == "HB2" ||
	  Name == "HD1" || Name == "HD2" ||
	  Name == "HE1" ) ) ||
      ( ( aminoCode == 'I' ) &&
	( Name == "HB" || 
	  Name == "HG11" || Name == "HG12" ) ) ||
      ( ( aminoCode == 'L' ) &&
	( Name == "HB1" || Name == "HB2" ||
	  Name == "HG" ) ) ||
      ( ( aminoCode == 'K' ) &&
	( Name == "HB1" || Name == "HB2" ||
	  Name == "HG1" || Name == "HG2" ||
	  Name == "HD1" || Name == "HD2" ||
	  Name == "HE1" || Name == "HE2" ) ) ||
      ( ( aminoCode == 'M' ) &&
	( Name == "HB1" || Name == "HB2" ||
	  Name == "HG1" || Name == "HG2" ) ) ||
      ( ( aminoCode == 'F' ) &&
	( Name == "HB1" || Name == "HB2" ||
	  Name == "HD1" || Name == "HD2" ||
	  Name == "HE1" || Name == "HE2" ||
	  Name == "HZ" ) ) ||
      ( ( aminoCode == 'P' ) &&
	( Name == "HB1" || Name == "HB2" ||
	  Name == "HG1" || Name == "HG2" ||
	  Name == "HD1" || Name == "HD2" ) ) ||
      ( ( aminoCode == 'S' ) &&
	( Name == "HB1" || Name == "HB2" ) ) ||
      ( ( aminoCode == 'T' ) &&
	( Name == "HB" ) ) ||
      ( ( aminoCode == 'W' ) &&
	( Name == "HB1" || Name == "HB2" ||
	  Name == "HD1" || 
	  Name == "HE1" || Name == "HE3" ||
	  Name == "HZ2" || Name == "HZ3" ||
	  Name == "HH2" ) ) ||
      ( ( aminoCode == 'Y' ) &&
	( Name == "HB1" || Name == "HB2" ||
	  Name == "HD1" || Name == "HD2" ||
	  Name == "HE1" || Name == "HE2" ) ) ||
      ( ( aminoCode == 'V' ) &&
	( Name == "HB" ) ) ){
    leak = LEAK0;
  }
  if( ( aminoCode == 'A' &&
        ( Name == "HB1" || Name == "HB2" || Name == "HB3" ) ) ||
      ( ( aminoCode == 'I' ) &&
	( Name == "HD11" || Name == "HD12" || Name == "HD13" ||
	  Name == "HG21" || Name == "HG22" || Name == "HG23" ) ) ||
      ( ( aminoCode == 'L' ) &&
	( Name == "HD11" || Name == "HD12" || Name == "HD13" ||
	  Name == "HD21" || Name == "HD22" || Name == "HD23" ) ) ||
      ( ( aminoCode == 'M' ) && 
	( Name == "HE1" || Name == "HE2" || Name == "HE3" ) ) ||
      ( ( aminoCode == 'T' ) &&
	( Name == "HG21" || Name == "HG22" || Name == "HG23" ) ) ||
      ( ( aminoCode == 'V' ) &&
	( Name == "HG11" || Name == "HG12" || Name == "HG13" ||
	  Name == "HG21" || Name == "HG22" || Name == "HG23" ) )
      ){
    leak = LEAK1;
  }
  if( Name == "H" || Name == "HE" ){
    leak = LEAK2;
  }
  fold2 = 1.0;
  if( Name == "H" ||
      (aminoCode == 'R' && Name == "HE") ){
    linewidth = LW_HN;
  }
  else{
    linewidth = LW_H;
  }
  lw1 = NAN;
  lw2 = NAN;
  valvef = NAN;
}

#define ALL_N 94.0

#define ALA_A 142.0
#define ALA_B 125.0

#define ARG_A 142.0
#define ARG_B 120.0
#define ARG_G 120.0
#define ARG_D 130.0

#define ASN_A 142.0
#define ASN_B 130.0

#define ASP_A 142.0
#define ASP_B 130.0

#define CYS_A 142.0
#define CYS_B 135.0

#define GLN_A 142.0
#define GLN_B 120.0
#define GLN_G 130.0

#define GLU_A 142.0
#define GLU_B 120.0
#define GLU_G 130.0
 
#define GLY_A 142.0

#define HIS_A 142.0
#define HIS_B 125.0
#define HIS_D1 100.0
#define HIS_D2 198.0
#define HIS_E1 218.0
#define HIS_E2 100.0

#define ILE_A 142.0
#define ILE_B 115.0
#define ILE_G1 120.0
#define ILE_D1 125.0
#define ILE_G2 125.0

#define LEU_A 142.0
#define LEU_B 120.0
#define LEU_G 115.0
#define LEU_D 125.0

#define LYS_A 142.0
#define LYS_B 120.0
#define LYS_G 120.0
#define LYS_D 120.0
#define LYS_E 130.0

#define MET_A 142.0
#define MET_B 120.0
#define MET_G 125.0
#define MET_E 130.0

#define PHE_A 142.0
#define PHE_B 125.0
#define PHE_D 158.0
#define PHE_E 158.0
#define PHE_Z 158.0

#define PRO_A 148.0
#define PRO_B 120.0
#define PRO_G 120.0
#define PRO_D 130.0

#define SER_A 142.0
#define SER_B 130.0

#define THR_A 142.0
#define THR_B 130.0
#define THR_G 125.0

#define TRP_A 142.0
#define TRP_B 130.0
#define TRP_D 180.0
#define TRP_E1 100.0
#define TRP_Z2 158.0
#define TRP_H 158.0
#define TRP_E3 158.0
#define TRP_Z3 158.0

#define TYR_A 142.0
#define TYR_B 125.0
#define TYR_D 158.0
#define TYR_E 158.0
#define VAL_A 142.0
#define VAL_B 115.0
#define VAL_G 125.0

hetero::hetero( string Name, char aminoCode, amino* residue ){
  this->Name = Name;
  shift = NAN;
  spec_shift2 = NAN;
  spec_shift3 = NAN;
  this->residue = residue;
  x = NAN;
  y = NAN;
  z = NAN;
  fold2 = NAN;
  fold3 = NAN;
  linewidth = NAN;
  lw2 = NAN;
  lw3 = NAN;
  valveg = NAN;

  if( Name[0] == 'N' &&
      !( aminoCode == 'H' && Name == "ND1" ) ){
    linewidth = LW_N;
  }
  // Is this exception correct?
  if( aminoCode == 'H' && Name == "ND1" ){
    linewidth = LW_C;
  }

  if( Name[0] == 'C' && 
      !( aminoCode == 'P' && Name == "CA" ) ){
    linewidth = LW_C;
  }
  // Is this exception correct?
  if( aminoCode == 'P' && Name == "CA" ){
    linewidth = LW_N;
  }

  rj = NAN;
  switch( aminoCode ){
  case 'A':
    if( Name == "N" ){
      rj = ALL_N;
    }
    if( Name == "CA" ){
      rj = ALA_A;
    }
    if( Name == "CB" ){
      rj = ALA_B;
    }
    break;
  case 'R':
    if( Name == "N" ){
      rj = ALL_N;
    }
    if( Name == "CA" ){
      rj = ARG_A;
    }
    if( Name == "CB" ){
      rj = ARG_B;
    }
    if( Name == "CG" ){
      rj = ARG_G;
    }
    if( Name == "CD" ){
      rj = ARG_D;
    }
    if( Name == "NE" ){
      rj = ALL_N;
    }
    break;
  case 'N':
    if( Name == "N" ){
      rj = ALL_N;
    }
    if( Name == "CA" ){
      rj = ASN_A;
    }
    if( Name == "CB" ){
      rj = ASN_B;
    }
    break;
  case 'D':
    if( Name == "N" ){
      rj = ALL_N;
    }
    if( Name == "CA" ){
      rj = ASP_A;
    }
    if( Name == "CB" ){
      rj = ASP_B;
    }
    break;
  case 'C':
    if( Name == "N" ){
      rj = ALL_N;
    }
    if( Name == "CA" ){
      rj = CYS_A;
    }
    if( Name == "CB" ){
      rj = CYS_B;
    }
    break;
  case 'Q':
    if( Name == "N" ){
      rj = ALL_N;
    }
    if( Name == "CA" ){
      rj = GLN_A;
    }
    if( Name == "CB" ){
      rj = GLN_B;
    }
    if( Name == "CG" ){
      rj = GLN_G;
    }
    break;
  case 'E':
    if( Name == "N" ){
      rj = ALL_N;
    }
    if( Name == "CA" ){
      rj = GLU_A;
    }
    if( Name == "CB" ){
      rj = GLU_B;
    }
    if( Name == "CG" ){
      rj = GLU_G;
    }
    break;
  case 'G':
    if( Name == "N" ){
      rj = ALL_N;
    }
    if( Name == "CA" ){
      rj = GLY_A;
    }
    break;
  case 'H':
    if( Name == "N" ){
      rj = ALL_N;
    }
    if( Name == "CA" ){
      rj = HIS_A;
    }
    if( Name == "CB" ){
      rj = HIS_B;
    }
    if( Name == "ND1" ){
      rj = HIS_D1;
    }
    if( Name == "CD2" ){
      rj = HIS_D2;
    }
    if( Name == "CE1" ){
      rj = HIS_E1;
    }
    break;
  case 'I':
    if( Name == "N" ){
      rj = ALL_N;
    }
    if( Name == "CA" ){
      rj = ILE_A;
    }
    if( Name == "CB" ){
      rj = ILE_B;
    }
    if( Name == "CG1" ){
      rj = ILE_G1;
    }
    if( Name == "CD1" ){
      rj = ILE_D1;
    }
    if( Name == "CG2" ){
      rj = ILE_G2;
    }
    break;
  case 'L':
    if( Name == "N" ){
      rj = ALL_N;
    }
    if( Name == "CA" ){
      rj = LEU_A;
    }
    if( Name == "CB" ){
      rj = LEU_B;
    }
    if( Name == "CG" ){
      rj = LEU_G;
    }
    if( Name == "CD1" ){
      rj = LEU_D;
    }
    if( Name == "CD2" ){
      rj = LEU_D;
    }
    break;
  case 'K':
    if( Name == "N" ){
      rj = ALL_N;
    }
    if( Name == "CA" ){
      rj = LYS_A;
    }
    if( Name == "CB" ){
      rj = LYS_B;
    }
    if( Name == "CG" ){
      rj = LYS_G;
    }
    if( Name == "CD" ){
      rj = LYS_D;
    }
    if( Name == "CE" ){
      rj = LYS_E;
    }
    break;
  case 'M':
    if( Name == "N" ){
      rj = ALL_N;
    }
    if( Name == "CA" ){
      rj = MET_A;
    }
    if( Name == "CB" ){
      rj = MET_B;
    }
    if( Name == "CG" ){
      rj = MET_G;
    }
    if( Name == "CE" ){
      rj = MET_E;
    }
    break;
  case 'F':
    if( Name == "N" ){
      rj = ALL_N;
    }
    if( Name == "CA" ){
      rj = PHE_A;
    }
    if( Name == "CB" ){
      rj = PHE_B;
    }
    if( Name == "CD1" ){
      rj = PHE_D;
    }
    if( Name == "CD2" ){
      rj = PHE_D;
    }
    if( Name == "CE1" ){
      rj = PHE_E;
    }
    if( Name == "CE2" ){
      rj = PHE_E;
    }
    if( Name == "CZ" ){
      rj = PHE_Z;
    }
    break;
  case 'P':
    if( Name == "CA" ){
      // FIXME: Is this really correct?
      rj = ALL_N;
    }
    if( Name == "CB" ){
      rj = PRO_B;
    }
    if( Name == "CG" ){
      rj = PRO_G;
    }
    if( Name == "CD" ){
      rj = PRO_D;
    }
    break;
  case 'S':
    if( Name == "N" ){
      rj = ALL_N;
    }
    if( Name == "CA" ){
      rj = SER_A;
    }
    if( Name == "CB" ){
      rj = SER_B;
    }
    break;
  case 'T':
    if( Name == "N" ){
      rj = ALL_N;
    }
    if( Name == "CA" ){
      rj = THR_A;
    }
    if( Name == "CB" ){
      rj = THR_B;
    }
    if( Name == "CG2" ){
      rj = THR_G;
    }
    break;
  case 'W':
    if( Name == "N" ){
      rj = ALL_N;
    }
    if( Name == "CA" ){
      rj = TRP_A;
    }
    if( Name == "CB" ){
      rj = TRP_B;
    }
    if( Name == "CD1" ){
      rj = TRP_D;
    }
    if( Name == "NE1" ){
      rj = TRP_E1;
    }
    if( Name == "CZ2" ){
      rj = TRP_Z2;
    }
    if( Name == "CH2" ){
      rj = TRP_H;
    }
    if( Name == "CE3" ){
      rj = TRP_E3;
    }
    if( Name == "CZ3" ){
      rj = TRP_Z3;
    }
    break;
  case 'Y':
    if( Name == "N" ){
      rj = ALL_N;
    }
    if( Name == "CA" ){
      rj = TYR_A;
    }
    if( Name == "CB" ){
      rj = TYR_B;
    }
    if( Name == "CD1" ){
      rj = TYR_D;
    }
    if( Name == "CD2" ){
      rj = TYR_D;
    }
    if( Name == "CE1" ){
      rj = TYR_E;
    }
    if( Name == "CE2" ){
      rj = TYR_E;
    }
    break;
  case 'V':
    if( Name == "N" ){
      rj = ALL_N;
    }
    if( Name == "CA" ){
      rj = VAL_A;
    }
    if( Name == "CB" ){
      rj = VAL_B;
    }
    if( Name == "CG1" ){
      rj = VAL_G;
    }
    if( Name == "CG2" ){
      rj = VAL_G;
    }
    break;
  default:
    break;
  }
}

void hetero::addproton( proton* h ){
  prot.push_back( h );
  h->hetatom = this;
}

amino::amino( char residueName, int residueNumber ){
  this->residueName = residueName;
  this->residueNumber = residueNumber;
}

void amino::addC13( hetero* het ){
  C13.push_back( het );
}

void amino::addN15( hetero* het ){
  N15.push_back( het );
}

protein::protein(){
  Name = "undefined";
}
