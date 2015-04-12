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
#include "NMRdata.h"
#include "BuildSequence.h"

using namespace std;

protein readSequence( const char* filename, int seq_start ){
  ifstream input( filename );
  if( !input ){
    cout << "Can't open " << filename << "." << endl;
    exit(1);
  }
  int seq_index = seq_start;

  protein prot;

  char aminochar;

  while( true ){
    input >> aminochar;

    if( !input.good() )
      break;

    amino* buildamino = new amino( aminochar, seq_index++ );
    prot.chain.push_back( buildamino );

    hetero* buildhetero;

    switch( aminochar ){
      
    case 'A': // ALA
      // add N, H(N), CA, HA, CB, HB1, HB2, HB3
      buildhetero = new hetero( "N", aminochar, buildamino );
      buildhetero->addproton( new proton( "H", 'A' ) );
      buildamino->addN15( buildhetero );

      buildhetero = new hetero( "CA", aminochar, buildamino );
      buildhetero->addproton( new proton( "HA", 'A' ) );
      buildamino->addC13( buildhetero );
      
      buildhetero = new hetero ( "CB", aminochar, buildamino );
      buildhetero->addproton( new proton( "HB1", 'A' ) );
      buildhetero->addproton( new proton( "HB2", 'A' ) );
      buildhetero->addproton( new proton( "HB3", 'A' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "C", aminochar, buildamino );
      buildamino->addC13( buildhetero );
      break;

    case 'C': // CYS
      // add N, H(N), CA, HA, CB, HB1, HB2
      buildhetero = new hetero( "N", aminochar, buildamino );
      buildhetero->addproton( new proton( "H", 'C' ) );
      buildamino->addN15( buildhetero );
      
      buildhetero = new hetero( "CA", aminochar, buildamino );
      buildhetero->addproton( new proton( "HA", 'C' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CB", aminochar, buildamino );
      buildhetero->addproton( new proton( "HB1", 'C' ) );
      buildhetero->addproton( new proton( "HB2", 'C' ) );
      buildamino->addC13( buildhetero );

      // TODO: not tested, some atoms might still be missing
      buildhetero = new hetero( "C", aminochar, buildamino );
      buildamino->addC13( buildhetero );
      break;

    case 'D': // ASP
      // add N, H(N), CA, HA, CB, HB1, HB2
      buildhetero = new hetero( "N", aminochar, buildamino );
      buildhetero->addproton( new proton( "H", 'D' ) );
      buildamino->addN15( buildhetero );

      buildhetero = new hetero( "CA", aminochar, buildamino );
      buildhetero->addproton( new proton( "HA", 'D' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CB", aminochar, buildamino );
      buildhetero->addproton( new proton( "HB1", 'D' ) );
      buildhetero->addproton( new proton( "HB2", 'D' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "CG", aminochar, buildamino );
      buildamino->addC13( buildhetero );
      
      buildhetero = new hetero( "C", aminochar, buildamino );
      buildamino->addC13( buildhetero );
      break;

    case 'E': // GLU
      // add N, H(N), CA, HA, CB, HB1, HB2, CG, HG1, HG2
      buildhetero = new hetero( "N", aminochar, buildamino );
      buildhetero->addproton( new proton( "H", 'E' ) );
      buildamino->addN15( buildhetero );

      buildhetero = new hetero( "CA", aminochar, buildamino );
      buildhetero->addproton( new proton( "HA", 'E' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CB", aminochar, buildamino );
      buildhetero->addproton( new proton( "HB1", 'E' ) );
      buildhetero->addproton( new proton( "HB2", 'E' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CG", aminochar, buildamino );
      buildhetero->addproton( new proton( "HG1", 'E' ) );
      buildhetero->addproton( new proton( "HG2", 'E' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "CD", aminochar, buildamino );
      buildamino->addC13( buildhetero );
      
      buildhetero = new hetero( "C", aminochar, buildamino );
      buildamino->addC13( buildhetero );
      break;

    case 'F': // PHE
      // add N, H(N), CA, HA, CB, HB1, HB2, CD1, HD1, CD2, HD2, 
      // CE1, HE1, CE2, HE2
      buildhetero = new hetero( "N", aminochar, buildamino );
      buildhetero->addproton( new proton( "H", 'F' ) );
      buildamino->addN15( buildhetero );

      buildhetero = new hetero( "CA", aminochar, buildamino );
      buildhetero->addproton( new proton( "HA", 'F' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CB", aminochar, buildamino );
      buildhetero->addproton( new proton( "HB1", 'F' ) );
      buildhetero->addproton( new proton( "HB2", 'F' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "CG", aminochar, buildamino );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "CD1", aminochar, buildamino );
      buildhetero->addproton( new proton( "HD1", 'F' ) );
      buildamino->addC13( buildhetero );
      
      buildhetero = new hetero( "CD2", aminochar, buildamino );
      buildhetero->addproton( new proton( "HD2", 'F' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "CE1", aminochar, buildamino );
      buildhetero->addproton( new proton( "HE1", 'F' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "CE2", aminochar, buildamino );
      buildhetero->addproton( new proton( "HE2", 'F' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "CZ", aminochar, buildamino );
      buildhetero->addproton( new proton( "HZ", 'F' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "C", aminochar, buildamino );
      buildamino->addC13( buildhetero );
      break;

    case 'G': // GLY
      // add N, H(N), CA, HA, HA1, HA2
      buildhetero = new hetero( "N", aminochar, buildamino );
      buildhetero->addproton( new proton( "H", 'G' ) );
      buildamino->addN15( buildhetero );

      buildhetero = new hetero( "CA", aminochar, buildamino );
      buildhetero->addproton( new proton( "HA1", 'G' ) );
      buildhetero->addproton( new proton( "HA2", 'G' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "C", aminochar, buildamino );
      buildamino->addC13( buildhetero );
      break;

    case 'H': // HIS
      // add N, H(N), CA, HA, CB, HB1, HB2, CD2, HD2, ND1, HD1, 
      // NE2, HE2, CE1, HE1
      buildhetero = new hetero( "N", aminochar, buildamino );
      buildhetero->addproton( new proton( "H", aminochar ) );
      buildamino->addN15( buildhetero );

      buildhetero = new hetero( "CA", aminochar, buildamino );
      buildhetero->addproton( new proton( "HA", aminochar ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "CB", aminochar, buildamino );
      buildhetero->addproton( new proton( "HB1", aminochar ) );
      buildhetero->addproton( new proton( "HB2", aminochar ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "CG", aminochar, buildamino );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "ND1", aminochar, buildamino );
      buildhetero->addproton( new proton( "HD1", aminochar ) );
      buildamino->addN15( buildhetero );

      buildhetero = new hetero( "CD2", aminochar, buildamino );
      buildhetero->addproton( new proton( "HD2", aminochar ) );
      buildamino->addC13( buildhetero );
      
      buildhetero = new hetero( "CE1" , aminochar, buildamino);
      buildhetero->addproton( new proton( "HE1", aminochar ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "NE2", aminochar, buildamino );
      buildamino->addN15( buildhetero );

      buildhetero = new hetero( "C", aminochar, buildamino );
      buildamino->addC13( buildhetero );
      break;

    case 'I': // HIS
      // add N, H(N), CA, HA, CB, HB, CG1, HG11, HG12, CG2, HG21, HG22, HG23
      // CD1, HD11, HD12, HD13
      buildhetero = new hetero( "N", aminochar, buildamino );
      buildhetero->addproton( new proton( "H", 'I' ) );
      buildamino->addN15( buildhetero );

      buildhetero = new hetero( "CA", aminochar, buildamino );
      buildhetero->addproton( new proton( "HA", 'I' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CB", aminochar, buildamino );
      buildhetero->addproton( new proton( "HB", 'I' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CG1", aminochar, buildamino );
      buildhetero->addproton( new proton( "HG11", 'I' ) );
      buildhetero->addproton( new proton( "HG12", 'I' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CG2", aminochar, buildamino );
      buildhetero->addproton( new proton( "HG21", 'I' ) );
      buildhetero->addproton( new proton( "HG22", 'I' ) );
      buildhetero->addproton( new proton( "HG23", 'I' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CD1", aminochar, buildamino );
      buildhetero->addproton( new proton( "HD11", 'I' ) );
      buildhetero->addproton( new proton( "HD12", 'I' ) );
      buildhetero->addproton( new proton( "HD13", 'I' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "C", aminochar, buildamino );
      buildamino->addC13( buildhetero );
      break;

    case 'K': // LYS
      // add N, H(N), CA, HA, CB, HB1, HB2, CG, HG1, HG2, CD, HD1, HD2
      // CE, HE1, HE2, NZ, HZ1, HZ2, HZ3
      buildhetero = new hetero( "N", aminochar, buildamino );
      buildhetero->addproton( new proton( "H", 'K' ) );
      buildamino->addN15( buildhetero );

      buildhetero = new hetero( "CA", aminochar, buildamino );
      buildhetero->addproton( new proton( "HA", 'K' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CB", aminochar, buildamino );
      buildhetero->addproton( new proton( "HB1", 'K' ) );
      buildhetero->addproton( new proton( "HB2", 'K' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CG", aminochar, buildamino );
      buildhetero->addproton( new proton( "HG1", 'K' ) );
      buildhetero->addproton( new proton( "HG2", 'K' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CD", aminochar, buildamino );
      buildhetero->addproton( new proton( "HD1", 'K' ) );
      buildhetero->addproton( new proton( "HD2", 'K' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CE", aminochar, buildamino );
      buildhetero->addproton( new proton( "HE1", 'K' ) );
      buildhetero->addproton( new proton( "HE2", 'K' ) );
      buildamino->addC13( buildhetero );

      // Commented out because of missing leak parameters
      // TODO: get leak parameters and include it again

//       buildhetero = new hetero( "NZ", aminochar, buildamino );
//       buildhetero->addproton( new proton( "HZ1", 'K' ) );
//       buildhetero->addproton( new proton( "HZ2", 'K' ) );
//       buildhetero->addproton( new proton( "HZ3", 'K' ) );
//       buildamino->addN15( buildhetero );

      buildhetero = new hetero( "C", aminochar, buildamino );
      buildamino->addC13( buildhetero );
      break;

    case 'L': // LEU
      // add N, H(N), CA, HA, CB, HB1, HB2, CG, HG, CD1, HD11, HD12, HD13,
      // CD2, HD21, HD22, HD23
      buildhetero = new hetero( "N", aminochar, buildamino );
      buildhetero->addproton( new proton( "H", 'L' ) );
      buildamino->addN15( buildhetero );

      buildhetero = new hetero( "CA", aminochar, buildamino );
      buildhetero->addproton( new proton( "HA", 'L' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CB", aminochar, buildamino );
      buildhetero->addproton( new proton( "HB1", 'L' ) );
      buildhetero->addproton( new proton( "HB2", 'L' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CG", aminochar, buildamino );
      buildhetero->addproton( new proton( "HG", 'L' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CD1", aminochar, buildamino );
      buildhetero->addproton( new proton( "HD11", 'L' ) );
      buildhetero->addproton( new proton( "HD12", 'L' ) );
      buildhetero->addproton( new proton( "HD13", 'L' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CD2", aminochar, buildamino );
      buildhetero->addproton( new proton( "HD21", 'L' ) );
      buildhetero->addproton( new proton( "HD22", 'L' ) );
      buildhetero->addproton( new proton( "HD23", 'L' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "C", aminochar, buildamino );
      buildamino->addC13( buildhetero );
      break;

    case 'M': // MET
      // add N, H(N), CA, HA, CB, HB1, HB2, CG, HG1, HG2, CE, HE1, HE2, HE3
      buildhetero = new hetero( "N", aminochar, buildamino );
      buildhetero->addproton( new proton( "H", 'M' ) );
      buildamino->addN15( buildhetero );

      buildhetero = new hetero( "CA", aminochar, buildamino );
      buildhetero->addproton( new proton( "HA", 'M' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CB", aminochar, buildamino );
      buildhetero->addproton( new proton( "HB1", 'M' ) );
      buildhetero->addproton( new proton( "HB2", 'M' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CG", aminochar, buildamino );
      buildhetero->addproton( new proton( "HG1", 'M' ) );
      buildhetero->addproton( new proton( "HG2", 'M' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CE", aminochar, buildamino );
      buildhetero->addproton( new proton( "HE1", 'M' ) );
      buildhetero->addproton( new proton( "HE2", 'M' ) );
      buildhetero->addproton( new proton( "HE3", 'M' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "C", aminochar, buildamino );
      buildamino->addC13( buildhetero );
      break;

    case 'N': // ASN
      // add N, H(N), CA, HA, CB, HB1, HB2, ND2, HD21, HD22
      buildhetero = new hetero( "N", aminochar, buildamino );
      buildhetero->addproton( new proton( "H", 'N' ) );
      buildamino->addN15( buildhetero );

      buildhetero = new hetero( "CA", aminochar, buildamino );
      buildhetero->addproton( new proton( "HA", 'N' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "CB", aminochar, buildamino );
      buildhetero->addproton( new proton( "HB1", 'N' ) );
      buildhetero->addproton( new proton( "HB2", 'N' ) );
      buildamino->addC13( buildhetero );
      
      buildhetero = new hetero( "CG", aminochar, buildamino );
      buildamino->addC13( buildhetero );

      // Commented out because of missing leak parameters
      // TODO: get leak parameters and include it again

//       buildhetero = new hetero( "ND2", aminochar, buildamino );
//       buildhetero->addproton( new proton( "HD21", 'N' ) );
//       buildhetero->addproton( new proton( "HD22", 'N' ) );
//       buildamino->addN15( buildhetero );

      buildhetero = new hetero( "C", aminochar, buildamino );
      buildamino->addC13( buildhetero );
      break;

    case 'P': // PRO
      // add N, H(N), CA, HA, CB, HB1, HB2, CG, HG1, HG2, CD, HD1, HD2
      buildhetero = new hetero( "N", aminochar, buildamino );
      buildamino->addN15( buildhetero );

      buildhetero = new hetero( "CA", aminochar, buildamino );
      buildhetero->addproton( new proton( "HA", 'P' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CB", aminochar, buildamino );
      buildhetero->addproton( new proton( "HB1", 'P' ) );
      buildhetero->addproton( new proton( "HB2", 'P' ) );
      buildamino->addC13( buildhetero );
      
      buildhetero = new hetero ( "CG", aminochar, buildamino );
      buildhetero->addproton( new proton( "HG1", 'P' ) );
      buildhetero->addproton( new proton( "HG2", 'P' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CD", aminochar, buildamino );
      buildhetero->addproton( new proton( "HD1", 'P' ) );
      buildhetero->addproton( new proton( "HD2", 'P' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "C", aminochar, buildamino );
      buildamino->addC13( buildhetero );
      break;

    case 'Q': // GLN
      // add N, H(N), CA, HA, CB, HB1, HB2, CG, HG1, HG2, NE2, HE21, HE22
      buildhetero = new hetero( "N", aminochar, buildamino );
      buildhetero->addproton( new proton( "H", 'Q' ) );
      buildamino->addN15( buildhetero );

      buildhetero = new hetero( "CA", aminochar, buildamino );
      buildhetero->addproton( new proton( "HA", 'Q' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CB", aminochar, buildamino );
      buildhetero->addproton( new proton( "HB1", 'Q' ) );
      buildhetero->addproton( new proton( "HB2", 'Q' ) );
      buildamino->addC13( buildhetero );
      
      buildhetero = new hetero ( "CG", aminochar, buildamino );
      buildhetero->addproton( new proton( "HG1", 'Q' ) );
      buildhetero->addproton( new proton( "HG2", 'Q' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "CD", aminochar, buildamino );
      buildamino->addC13( buildhetero );

      // Commented out because of missing leak parameters
      // TODO: get leak parameters and include it again

//       buildhetero = new hetero ( "NE2", aminochar, buildamino );
//       buildhetero->addproton( new proton( "HE21", 'Q' ) );
//       buildhetero->addproton( new proton( "HE22", 'Q' ) );
//       buildamino->addN15( buildhetero );

      buildhetero = new hetero( "C", aminochar, buildamino );
      buildamino->addC13( buildhetero );
      break;

    case 'R': // ARG
      // add N, H(N), CA, HA, CB, HB1, HB2, CG, HG1, HG2, CD, HD1, HD2, NE, HE,
      // NH1, HH11, HH12, NH2, HH21, HH22
      buildhetero = new hetero( "N", aminochar, buildamino );
      buildhetero->addproton( new proton( "H", 'R' ) );
      buildamino->addN15( buildhetero );

      buildhetero = new hetero( "CA", aminochar, buildamino );
      buildhetero->addproton( new proton( "HA", 'R' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CB", aminochar, buildamino );
      buildhetero->addproton( new proton( "HB1", 'R' ) );
      buildhetero->addproton( new proton( "HB2", 'R' ) );
      buildamino->addC13( buildhetero );
      
      buildhetero = new hetero ( "CG", aminochar, buildamino );
      buildhetero->addproton( new proton( "HG1", 'R' ) );
      buildhetero->addproton( new proton( "HG2", 'R' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CD", aminochar, buildamino );
      buildhetero->addproton( new proton( "HD1", 'R' ) );
      buildhetero->addproton( new proton( "HD2", 'R' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "NE", aminochar, buildamino );
      buildhetero->addproton( new proton( "HE", 'R' ) );
      buildamino->addN15( buildhetero );

      buildhetero = new hetero( "CZ", aminochar, buildamino );
      buildamino->addC13( buildhetero );

      // Commented out because of missing leak parameters
      // TODO: get leak parameters and include it again

//       buildhetero = new hetero ( "NH1", aminochar, buildamino );
//       buildhetero->addproton( new proton( "HH11", 'R' ) );
//       buildhetero->addproton( new proton( "HH12", 'R' ) );
//       buildamino->addN15( buildhetero );

//       buildhetero = new hetero ( "NH2", aminochar, buildamino );
//       buildhetero->addproton( new proton( "HH21", 'R' ) );
//       buildhetero->addproton( new proton( "HH22", 'R' ) );
//       buildamino->addN15( buildhetero );

      buildhetero = new hetero( "C", aminochar, buildamino );
      buildamino->addC13( buildhetero );
      break;

    case 'S': // SER
      // add N, H(N), CA, HA, CB, HB1, HB2
      buildhetero = new hetero( "N", aminochar, buildamino );
      buildhetero->addproton( new proton( "H", 'S' ) );
      buildamino->addN15( buildhetero );

      buildhetero = new hetero( "CA", aminochar, buildamino );
      buildhetero->addproton( new proton( "HA", 'S' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CB", aminochar, buildamino );
      buildhetero->addproton( new proton( "HB1", 'S' ) );
      buildhetero->addproton( new proton( "HB2", 'S' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "C", aminochar, buildamino );
      buildamino->addC13( buildhetero );
      break;

    case 'T': // THR
      // add N, H(N), CA, HA, CB, HB, CG2, HG21, HG22, HG23
      buildhetero = new hetero( "N", aminochar, buildamino );
      buildhetero->addproton( new proton( "H", 'T' ) );
      buildamino->addN15( buildhetero );

      buildhetero = new hetero( "CA", aminochar, buildamino );
      buildhetero->addproton( new proton( "HA", 'T' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CB", aminochar, buildamino );
      buildhetero->addproton( new proton( "HB", 'T' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CG2", aminochar, buildamino );
      buildhetero->addproton( new proton( "HG21", 'T' ) );
      buildhetero->addproton( new proton( "HG22", 'T' ) );
      buildhetero->addproton( new proton( "HG23", 'T' ) );
      buildamino->addC13( buildhetero );
      
      buildhetero = new hetero( "C", aminochar, buildamino );
      buildamino->addC13( buildhetero );
      break;

    case 'V': // VAL
      // add N, H(N), CA, HA, CB, HB, CG1, HG11, HG12, HG13, 
      // CG2, HG21, HG22, HG23
      buildhetero = new hetero( "N", aminochar, buildamino );
      buildhetero->addproton( new proton( "H", 'V' ) );
      buildamino->addN15( buildhetero );

      buildhetero = new hetero( "CA", aminochar, buildamino );
      buildhetero->addproton( new proton( "HA", 'V' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CB", aminochar, buildamino );
      buildhetero->addproton( new proton( "HB", 'V' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CG1", aminochar, buildamino );
      buildhetero->addproton( new proton( "HG11", 'V' ) );
      buildhetero->addproton( new proton( "HG12", 'V' ) );
      buildhetero->addproton( new proton( "HG13", 'V' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CG2", aminochar, buildamino );
      buildhetero->addproton( new proton( "HG21", 'V' ) );
      buildhetero->addproton( new proton( "HG22", 'V' ) );
      buildhetero->addproton( new proton( "HG23", 'V' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "C", aminochar, buildamino );
      buildamino->addC13( buildhetero );
      break;

    case 'W': // TRP
      // add N, H(N), CA, HA, CB, HB1, HB2, CD1, HD1,
      // NE1, HE1, CE3, HE3, CZ2, HZ2, CZ3, HZ3, CH2, HH2
      buildhetero = new hetero( "N", aminochar, buildamino );
      buildhetero->addproton( new proton( "H", 'W' ) );
      buildamino->addN15( buildhetero );

      buildhetero = new hetero( "CA", aminochar, buildamino );
      buildhetero->addproton( new proton( "HA", 'W' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CB", aminochar, buildamino );
      buildhetero->addproton( new proton( "HB1", 'W' ) );
      buildhetero->addproton( new proton( "HB2", 'W' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "CD1", aminochar, buildamino );
      buildhetero->addproton( new proton( "HD1", 'W' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "NE1", aminochar, buildamino );
      buildhetero->addproton( new proton( "HE1", 'W' ) );
      buildamino->addN15( buildhetero );

      buildhetero = new hetero( "CE3", aminochar, buildamino );
      buildhetero->addproton( new proton( "HE3", 'W' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "CZ2", aminochar, buildamino );
      buildhetero->addproton( new proton( "HZ2", 'W' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "CZ3", aminochar, buildamino );
      buildhetero->addproton( new proton( "HZ3", 'W' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "CH2", aminochar, buildamino );
      buildhetero->addproton( new proton( "HH2", 'W' ) );
      buildamino->addC13( buildhetero );
      break;

    case 'Y': // TYR
      // add N, H(N), CA, HA, CB, HB1, HB2, CD1, HD1, 
      // CD2, HD2, CE1, HE1, CE2, HE2
      buildhetero = new hetero( "N", aminochar, buildamino );
      buildhetero->addproton( new proton( "H", 'Y' ) );
      buildamino->addN15( buildhetero );

      buildhetero = new hetero( "CA", aminochar, buildamino );
      buildhetero->addproton( new proton( "HA", 'Y' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero ( "CB", aminochar, buildamino );
      buildhetero->addproton( new proton( "HB1", 'Y' ) );
      buildhetero->addproton( new proton( "HB2", 'Y' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "CD1", aminochar, buildamino );
      buildhetero->addproton( new proton( "HD1", 'Y' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "CD2", aminochar, buildamino );
      buildhetero->addproton( new proton( "HD2", 'Y' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "CE1", aminochar, buildamino );
      buildhetero->addproton( new proton( "HE1", 'Y' ) );
      buildamino->addC13( buildhetero );

      buildhetero = new hetero( "CE2", aminochar, buildamino );
      buildhetero->addproton( new proton( "HE2", 'Y' ) );
      buildamino->addC13( buildhetero );
      break;
    }
  }

  // Create indices
  vector<amino*>::iterator amin_it = prot.chain.begin();
  while( amin_it != prot.chain.end() ){

    vector<hetero*>::iterator het_it = (*amin_it)->C13.begin();
    while( het_it != (*amin_it)->C13.end() ){
      vector<proton*>::iterator prot_it = (*het_it)->prot.begin();
      while( prot_it != (*het_it)->prot.end() ){
	(*amin_it)->H.push_back( *prot_it );
	prot.H.push_back( *prot_it );
	prot.HC.push_back( *prot_it );
	prot_it++;
      }
      prot.C13.push_back( *het_it );
      het_it++;
    }

    het_it = (*amin_it)->N15.begin();
    while( het_it != (*amin_it)->N15.end() ){
      vector<proton*>::iterator prot_it = (*het_it)->prot.begin();
      while( prot_it != (*het_it)->prot.end() ){
	(*amin_it)->H.push_back( *prot_it );
	prot.H.push_back( *prot_it );
	prot.HN.push_back( *prot_it );
	prot_it++;
      }
      prot.N15.push_back( *het_it );
      het_it++;
    }

    amin_it++;
  }

  // Set valveg/valvef
  for( vector<proton*>::iterator prot_it = prot.H.begin();
       prot_it != prot.H.end(); prot_it++ ){
    if( (*prot_it)->hetatom->prot.size() == 3 ){
      (*prot_it)->valvef = 3.0;
      (*prot_it)->hetatom->valveg = 3.0;
    }
    else{
      (*prot_it)->valvef = 1.0;
      (*prot_it)->hetatom->valveg = 1.0;
    }
  }

  return prot;
}
