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
#include <vector>
#include <complex>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <pdbimport.h>
#include <amino.h>
#include <nmrFileIO.h>
#include "NMRparseparameters.h"
#include "Shiftfileparser.h"
#include "BuildSequence.h"
#include "NMRdata.h"

#define ZEROP 0.0000001
#define ONEM 0.9999999
#define CON_HH 56963.0
#define CON_HC 3603.21
#define CON_HN 585.414
#define HCMAX 6.0
#define HNMAX 3.0

using namespace std;

double fnJ( double x, double t ){
  return t / ( 1.0 + pow( 6.2832 * x * t, 2.0 ) );
}

int main( int argc, char** argv ){
  /* Parse parameters */
  if( argc != 2 ){
    cout << "Usage: Shine <parameterfile>" << endl;
    exit(1);
  }
  char* parameterfile = argv[1];


  /* Parse input files */
  
  /* First the parameter file */
  struct parameters param = NMRparseparameters( parameterfile );
  
  cout << "Parameters parsed." << endl;
  
  /* Read and build the sequence */
  protein prot = readSequence ( param.sequence.c_str(), param.seq_start );

  cout << "Sequence read." << endl;
  /* Then the PDB file */
  vector<pdbatom>* pdbdata = readatoms( param.pdbfile.c_str() );

  cout << "PDB file read." << endl;
  
  /* The Shift file (longlist format for now, ok changed to a
     special inputfile format with an explizit assignment to HB1/HB2 etc. */
  vector<shiftfileitem> shifts = readshiftfile( param.shiftfile.c_str() );

  cout << "Shifts read." << endl;

  /* reorganize data, preprocessing */

  // coordinate and axis configuration
  if( param.heteromode == "single" ){
    if( param.delay_proton ){
      param.cors2 = -1.0;
    }
    else{
      param.cors2 = 1.0;
    }
    if( param.delay_hetbound ){
      param.cors3 = -1.0;
    }
    else{
      param.cors3 = 1.0;
    }
  }
  if( param.heteromode == "double" ){
    if( param.delay_hetfree ){
      param.cors2 = -1.0;
    }
    else{
      param.cors2 = 1.0;
    }
    if( param.delay_hetbound ){
      param.cors3 = -1.0;
    }
    else{
      param.cors3 = 1.0;
    }
  }

  // Insert PDB coordinates into protein

  vector<amino*>::iterator resit;
  vector<pdbatom>::iterator pdbit = pdbdata->begin();
  vector<proton*>::iterator protit;
  vector<hetero*>::iterator hetit;

  while( pdbit != pdbdata->end() ){
    resit = prot.chain.begin();
    while( resit != prot.chain.end() ){
      if( (*resit)->residueNumber == pdbit->residnumber ){
	switch( pdbit->atomname.c_str()[0]){
	case 'H':
	  protit = (*resit)->H.begin();
	  while( protit != (*resit)->H.end() ){
	    if( (*protit)->Name == pdbit->atomname ){
	      (*protit)->x = pdbit->xcoord;
	      (*protit)->y = pdbit->ycoord;
	      (*protit)->z = pdbit->zcoord;
	    }
	    protit++;
	  }
	  break;
	case 'C':
	  hetit = (*resit)->C13.begin();
	  while( hetit != (*resit)->C13.end() ){
	    if( (*hetit)->Name == pdbit->atomname ){
	      (*hetit)->x = pdbit->xcoord;
	      (*hetit)->y = pdbit->ycoord;
	      (*hetit)->z = pdbit->zcoord;
	    }
	    hetit++;
	  }
	  break;
	case 'N':
	  hetit = (*resit)->N15.begin();
	  while( hetit != (*resit)->N15.end() ){
	    if( (*hetit)->Name == pdbit->atomname ){
	      (*hetit)->x = pdbit->xcoord;
	      (*hetit)->y = pdbit->ycoord;
	      (*hetit)->z = pdbit->zcoord;
	    }
	    hetit++;
	  }
	  break;
	case 'O':
	  break;
	case 'S':
	  break;
	default:
	  cout << "Problem while inserting Coordinates." << endl;
	  cout << "Atom Name: " << pdbit->atomname << endl;
	  break;
	}
      }
      resit++;
    }
    pdbit++;
  }

  cout << "PDB coordinates inserted into the sequence." << endl;

  // Insert Shifts into protein
  
  vector<shiftfileitem>::iterator sit = shifts.begin();
  while( sit != shifts.end() ){
    // Convert -999.900 values to NaN
    if( sit->shift < -990.0 ){
      sit->shift = NAN;
    }

    resit = prot.chain.begin();
    while( resit != prot.chain.end() ){
      if( (*resit)->residueNumber == sit->resid ){
	switch( sit->atomname[0] ){
	case 'H':
	  protit = (*resit)->H.begin();
	  while( protit != (*resit)->H.end() ){
	    if( (*protit)->Name == sit->atomname ){
	      (*protit)->shift = sit->shift;
	    }
	    protit++;
	  }
	  break;
	case 'C':
	  hetit = (*resit)->C13.begin();
	  while( hetit != (*resit)->C13.end() ){
	    if( (*hetit)->Name == sit->atomname ){
	      (*hetit)->shift = sit->shift;
	    }
	    hetit++;
	  }
	  break;
	case 'N':
	  hetit = (*resit)->N15.begin();
	  while( hetit != (*resit)->N15.end() ){
	    if( (*hetit)->Name == sit->atomname ){
	      (*hetit)->shift = sit->shift;
	    }
	    hetit++;
	  }
	  break;
	default:
	  cout << "Problem while inserting Shift: " << sit->atomname 
	       << " " << sit->resid 
	       << endl;
	  break;
	}
      }
      resit++;
    }
    sit++;
  }

  cout << "Shifts inserted into sequence." << endl;

  // Fix incomplete methyl group shifts

  for( protit = prot.H.begin(); protit != prot.H.end(); protit++ ){
    if( std::isnan( (*protit)->shift ) && (*protit)->hetatom->prot.size() == 3 ){
      for( unsigned int i = 0; i < 3; i++ ){
	if( !std::isnan( (*protit)->hetatom->prot[i]->shift ) ){
	  (*protit)->shift = (*protit)->hetatom->prot[i]->shift;
	}
      }
    }
  }

  // Calculate trans variables
  if( param.heterobound == "C13" ){
    hetit = prot.C13.begin();
    while( hetit != prot.C13.end() ){
      if( param.hetbound_sensitivity_enhanced ){
	(*hetit)->trans = 
	  sin( M_PI * (*hetit)->rj * param.inept_hetbound / 1000.0 ) * 
	  ( sin( M_PI * (*hetit)->rj * param.inept_hetbound_rev / 1000.0 ) +
	    sin( M_PI * (*hetit)->rj * param.inept_hetbound_MQ / 1000.0 ) );
      }
      else{
	(*hetit)->trans = 
	  sin( M_PI * (*hetit)->rj * param.inept_hetbound / 1000.0 ) * 
	  sin( M_PI * (*hetit)->rj * param.inept_hetbound_rev / 1000.0 );
      }
      hetit++;
    }
  }
  if( param.heterobound == "N15" ){
    hetit = prot.N15.begin();
    while( hetit != prot.N15.end() ){
      if( param.hetbound_sensitivity_enhanced ){
	(*hetit)->trans = 
	  sin( M_PI * (*hetit)->rj * param.inept_hetbound / 1000.0 ) * 
	  ( sin( M_PI * (*hetit)->rj * param.inept_hetbound_rev / 1000.0 ) +
	    sin( M_PI * (*hetit)->rj * param.inept_hetbound_MQ / 1000.0 ) );
      }
      else{
	(*hetit)->trans = 
	  sin( M_PI * (*hetit)->rj * param.inept_hetbound / 1000.0 ) * 
	  sin( M_PI * (*hetit)->rj * param.inept_hetbound_rev / 1000.0 );
      }
      hetit++;
    }
  }

  if( param.heterofree == "C13" ){
    hetit = prot.C13.begin();
    while( hetit != prot.C13.end() ){
      if( param.hetfree_sensitivity_enhanced ){
	(*hetit)->trans = 
	  sin( M_PI * (*hetit)->rj * param.inept_hetfree / 1000.0 ) * 
	  ( sin( M_PI * (*hetit)->rj * param.inept_hetfree_rev / 1000.0 ) +
	    sin( M_PI * (*hetit)->rj * param.inept_hetfree_MQ / 1000.0 ) );
      }
      else{
	(*hetit)->trans = 
	  sin( M_PI * (*hetit)->rj * param.inept_hetfree / 1000.0 ) * 
	  sin( M_PI * (*hetit)->rj * param.inept_hetfree_rev / 1000.0 );
      }
      hetit++;
    }
  }
  if( param.heterofree == "N15" ){
    hetit = prot.N15.begin();
    while( hetit != prot.N15.end() ){
      if( param.hetfree_sensitivity_enhanced ){
	(*hetit)->trans = 
	  sin( M_PI * (*hetit)->rj * param.inept_hetfree / 1000.0 ) * 
	  ( sin( M_PI * (*hetit)->rj * param.inept_hetfree_rev / 1000.0 ) +
	    sin( M_PI * (*hetit)->rj * param.inept_hetfree_MQ / 1000.0 ) );
      }
      else{
	(*hetit)->trans = 
	  sin( M_PI * (*hetit)->rj * param.inept_hetfree / 1000.0 ) * 
	  sin( M_PI * (*hetit)->rj * param.inept_hetfree_rev / 1000.0 );
      }
      hetit++;
    }
  }
  
  cout << "trans variables calculated." << endl;

  // Calculate folding and inversion of hetero shifts
  hetit = prot.C13.begin();
  while( hetit != prot.C13.end() ){
    (*hetit)->fold3 = 1.0;
    if( !std::isnan( (*hetit)->shift ) ){
      (*hetit)->spec_shift2 = param.refpoint2 - 
	( (*hetit)->shift - param.refppm2 ) * 
	double( param.size2 ) / param.sweepwidth2 * param.freq2;
      // folded?
      if( (*hetit)->spec_shift2 < 0.0 ){
	(*hetit)->spec_shift2 += double( param.size2 );
	(*hetit)->fold2 = param.cors2;
      }
      // double folded?
      if( (*hetit)->spec_shift2 < 0.0 ){
	(*hetit)->spec_shift2 += double( param.size2 );
	(*hetit)->fold2 = param.cors2 * param.cors2;
      }

      // folded?
      if( (*hetit)->spec_shift2 > double( param.size2 ) ){
	(*hetit)->spec_shift2 -= double( param.size2 );
	(*hetit)->fold2 = param.cors2;
      }
      // double folded?
      if( (*hetit)->spec_shift2 > double( param.size2 ) ){
	(*hetit)->spec_shift2 -= double( param.size2 );
	(*hetit)->fold2 = param.cors2 * param.cors2;
      }
      

      (*hetit)->spec_shift3 = param.refpoint3 - 
	( (*hetit)->shift - param.refppm3 ) * 
	double( param.size3 ) / param.sweepwidth3 * param.freq3;
      // folded?
      if( (*hetit)->spec_shift3 < 0.0 ){
	(*hetit)->spec_shift3 += double( param.size3 );
	(*hetit)->fold3 = param.cors3;
      }
      // double folded?
      if( (*hetit)->spec_shift3 < 0.0 ){
	(*hetit)->spec_shift3 += double( param.size3 );
	(*hetit)->fold3 = param.cors3 * param.cors3;
      }

      // folded?
      if( (*hetit)->spec_shift3 > double( param.size3 ) ){
	(*hetit)->spec_shift3 -= double( param.size3 );
	(*hetit)->fold3 = param.cors3;
      }
      // double folded?
      if( (*hetit)->spec_shift3 > double( param.size3 ) ){
	(*hetit)->spec_shift3 -= double( param.size3 );
	(*hetit)->fold3 = param.cors3 * param.cors3;
      }
    }
    hetit++;
  }
  hetit = prot.N15.begin();
  while( hetit != prot.N15.end() ){
    (*hetit)->fold3 = 1.0;
    if( !std::isnan( (*hetit)->shift ) ){
      (*hetit)->spec_shift2 = param.refpoint2 - 
	( (*hetit)->shift - param.refppm2 ) * 
	double( param.size2 ) / param.sweepwidth2 * param.freq2;
      // folded?
      if( (*hetit)->spec_shift2 < 0.0 ){
	(*hetit)->spec_shift2 += double( param.size2 );
	(*hetit)->fold2 = param.cors2;
      }
      // double folded?
      if( (*hetit)->spec_shift2 < 0.0 ){
	(*hetit)->spec_shift2 += double( param.size2 );
	(*hetit)->fold2 = param.cors2 * param.cors2;
      }

      // folded?
      if( (*hetit)->spec_shift2 > double( param.size2 ) ){
	(*hetit)->spec_shift2 -= double( param.size2 );
	(*hetit)->fold2 = param.cors2;
      }
      // double folded?
      if( (*hetit)->spec_shift2 > double( param.size2 ) ){
	(*hetit)->spec_shift2 -= double( param.size2 );
	(*hetit)->fold2 = param.cors2 * param.cors2;
      }


      (*hetit)->spec_shift3 = param.refpoint3 - 
	( (*hetit)->shift - param.refppm3 ) * 
	double( param.size3 ) / param.sweepwidth3 * param.freq3;
      // folded?
      if( (*hetit)->spec_shift3 < 0.0 ){
	(*hetit)->spec_shift3 += double( param.size3 );
	(*hetit)->fold3 = param.cors3;
      }
      // double folded?
      if( (*hetit)->spec_shift3 < 0.0 ){
	(*hetit)->spec_shift3 += double( param.size3 );
	(*hetit)->fold3 = param.cors3 * param.cors3;
      }

      // folded?
      if( (*hetit)->spec_shift3 > double( param.size3 ) ){
	(*hetit)->spec_shift3 -= double( param.size3 );
	(*hetit)->fold3 = param.cors3;
      }
      // double folded?
      if( (*hetit)->spec_shift3 > double( param.size3 ) ){
	(*hetit)->spec_shift3 -= double( param.size3 );
	(*hetit)->fold3 = param.cors3 * param.cors3;
      }
    }
    hetit++;
  }

  // Calculate folding and inversion of proton shifts
  protit = prot.H.begin();
  while( protit != prot.H.end() ){
    if( !std::isnan( (*protit)->shift ) ){
      (*protit)->spec_shift1 = param.refpoint1 - 
	( (*protit)->shift - param.refppm1 ) * double( param.size1 ) / 
	param.sweepwidth1 * param.freq1;
      if( (*protit)->spec_shift1 < 0.0 ){
	(*protit)->spec_shift1 += double( param.size1 );
      }
      if( (*protit)->spec_shift1 > double( param.size1 ) ){
	(*protit)->spec_shift1 -= double( param.size1 );
      }

      (*protit)->spec_shift2 = param.refpoint2 - 
	( (*protit)->shift - param.refppm2 ) * double( param.size2 ) / 
	param.sweepwidth2 * param.freq2;
      if( (*protit)->spec_shift2 < 0.0 ){
	(*protit)->spec_shift2 += double( param.size2 );
	(*protit)->fold2 = param.cors2;
      }
      if( (*protit)->spec_shift2 > double( param.size2 ) ){
	(*protit)->spec_shift2 -= double( param.size2 );
	(*protit)->fold2 = param.cors2;
      }
    }
    protit++;
  }

  cout << "Folding and inversion of shifts calculated." << endl;

  // Calculate linewidths
  hetit = prot.C13.begin();
  while( hetit != prot.C13.end() ){
    (*hetit)->lw2 = (*hetit)->linewidth * param.linewidthfactor2;
    (*hetit)->lw3 = (*hetit)->linewidth * param.linewidthfactor3;
    hetit++;
  }
  hetit = prot.N15.begin();
  while( hetit != prot.N15.end() ){
    (*hetit)->lw2 = (*hetit)->linewidth * param.linewidthfactor2;
    (*hetit)->lw3 = (*hetit)->linewidth * param.linewidthfactor3;
    hetit++;
  }
  protit = prot.H.begin();
  while( protit != prot.H.end() ){
    (*protit)->lw1 = (*protit)->linewidth * param.linewidthfactor1;
    (*protit)->lw2 = (*protit)->linewidth * param.linewidthfactor2;
    protit++;
  }

  cout << "Linewidths calculated." << endl;

  // Only for the isotropic case !!!
  double tau1 = param.isocorrel;
  tau1 = tau1 / 1000.0;
  double D1 = 1.0 / ( 6.0 * tau1 );
  double D2 = D1;

  // Calculation of ta, tb, tc, whatever they are good for ...
  double ta = 1.0 / ( 6.0 * D2 );
  double tb = 1.0 / ( 5.0 * D2 + D1 );
  double tc = 1.0 / ( 2.0 * D2 + 4.0 * D1 );

  // Handling of symmetry axis
  // Predefining the atoms defining the symmetry axis
  double temp = -0.288675135;
  double xs1 =  temp + 1.0;
  double xs2 = -temp + 1.0;
  double ys1 =  temp + 2.0;
  double ys2 = -temp + 2.0;
  double zs1 =  temp + 3.0;
  double zs2 = -temp + 3.0;
  
  hetit = prot.C13.begin();
  while( hetit != prot.C13.end() ){
    if( (*hetit)->residue->residueNumber == param.symatom1_resID &&
	(*hetit)->Name == param.symatom1_Name ){
      xs1 = (*hetit)->x;
      ys1 = (*hetit)->y;
      zs1 = (*hetit)->z;
    }
    if( (*hetit)->residue->residueNumber == param.symatom2_resID &&
	(*hetit)->Name == param.symatom2_Name ){
      xs2 = (*hetit)->x;
      ys2 = (*hetit)->y;
      zs2 = (*hetit)->z;
    }
    hetit++;
  }
  
  hetit = prot.N15.begin();
  while( hetit != prot.N15.end() ){
    if( (*hetit)->residue->residueNumber == param.symatom1_resID &&
	(*hetit)->Name == param.symatom1_Name ){
      xs1 = (*hetit)->x;
      ys1 = (*hetit)->y;
      zs1 = (*hetit)->z;
    }
    if( (*hetit)->residue->residueNumber == param.symatom2_resID &&
	(*hetit)->Name == param.symatom2_Name ){
      xs2 = (*hetit)->x;
      ys2 = (*hetit)->y;
      zs2 = (*hetit)->z;
    }
    hetit++;
  }

  protit = prot.H.begin();
  while( protit != prot.H.end() ){
    if( (*protit)->hetatom->residue->residueNumber == param.symatom1_resID &&
	(*protit)->Name == param.symatom1_Name ){
      xs1 = (*protit)->x;
      ys1 = (*protit)->y;
      zs1 = (*protit)->z;
    }
    if( (*protit)->hetatom->residue->residueNumber == param.symatom2_resID &&
	(*protit)->Name == param.symatom2_Name ){
      xs2 = (*protit)->x;
      ys2 = (*protit)->y;
      zs2 = (*protit)->z;
    }
    protit++;
  }

  // original of the axis system
  double x0 = 0.5 * (xs2 + xs1);
  double y0 = 0.5 * (ys2 + ys1);
  double z0 = 0.5 * (zs2 + zs1);
  
  // S (z) axis
  double dxs = xs2 - xs1;
  double dys = ys2 - ys1;
  double dzs = zs2 - zs1;
  double ds = sqrt( dxs*dxs + dys*dys + dzs*dzs );
  dxs = dxs / ds;
  dys = dys / ds;
  dzs = dzs / ds;


  // P (x) axis
  double dxp = -x0;
  double dyp = -y0;
  double dzp = ( dxs*x0 + dys*y0 ) / dzs;
  double dp  = sqrt( dxp * dxp + dyp * dyp + dzp * dzp );
  dxp = dxp / dp;
  dyp = dyp / dp;
  dzp = dzp / dp;

  // Q (y) axis
  double dxq = z0 * ( dzs * dyp - dzp * dys ) / ( dxs * dyp - dxp * dys );
  double dyq = z0 * ( dzs * dxp - dzp * dxs ) / ( dys * dxp - dyp * dxs );
  double dzq = -z0;
  double dq = sqrt( dxq * dxq + dyq * dyq + dzq * dzq );
  dxq = dxq / dq;
  dyq = dyq / dq;
  dzq = dzq / dq;

  cout << "Calculation of various variables done." << endl;

  /* calculate Relax Matrix */
  
  /* allocate some memory for matrices */
  double* XR = (double*) calloc( prot.H.size() * prot.H.size(), 
				 sizeof(double) );
  double* XRh = (double*) calloc( prot.H.size() * prot.H.size(), 
				  sizeof(double) );
  double* XRc = (double*) calloc( prot.H.size() * prot.C13.size(), 
				  sizeof(double) );
  double* XRn = (double*) calloc( prot.H.size() * prot.N15.size(), 
				  sizeof(double) );
  double* RD = (double*) calloc( prot.H.size(), sizeof(double) );
  double* RG = (double*) calloc( prot.H.size() * prot.H.size(),
				 sizeof(double) );
  double* P = (double*) calloc( prot.H.size(), sizeof(double) );
  double* G = (double*) calloc( prot.H.size() * prot.H.size(),
				sizeof(double) );
  double* DD = (double*) calloc( prot.H.size(), sizeof(double) );

  for( unsigned int i = 0; i < prot.H.size(); i++ ){
    for( unsigned int j = i+1; j  < prot.H.size(); j++ ){

      // non-methyl or methyl on same hetero atom
      if( ( prot.H[i]->hetatom->prot.size() != 3 && 
	    prot.H[j]->hetatom->prot.size() != 3 ) || 
	  ( ( prot.H[i]->hetatom->prot.size() == 3 && 
	      prot.H[j]->hetatom->prot.size() == 3 ) &&
	    prot.H[i]->hetatom == prot.H[j]->hetatom ) ){
	// Calculate distance
	double dx = prot.H[i]->x - prot.H[j]->x;
	double dy = prot.H[i]->y - prot.H[j]->y;
	double dz = prot.H[i]->z - prot.H[j]->z;
	double r = sqrt( dx*dx + dy*dy + dz*dz );
	if( r < ZEROP || std::isnan(r) ){
	  cout << "Problems with the distance between two protons: " 
	       << r << endl;
	  cout << "Proton i: " 
	       << prot.H[i]->hetatom->residue->residueNumber << " "
	       << prot.H[i]->hetatom->residue->residueName << " "
	       << prot.H[i]->Name << " "
	       << prot.H[i]->x << ", "
	       << prot.H[i]->y << ", "
	       << prot.H[i]->z 
	       << ", Proton j: " 
	       << prot.H[j]->hetatom->residue->residueNumber << " "
	       << prot.H[j]->hetatom->residue->residueName << " "
	       << prot.H[j]->Name << " "
	       << prot.H[j]->x << ", "
	       << prot.H[j]->y << ", "
	       << prot.H[j]->z << endl; 
	  exit(EXIT_FAILURE);
	}
	double an = ( dxs*dx + dys*dx + dzs*dz ) / r;
	if( fabs(an) > ONEM ){
	  an = ONEM;
	}
	an = acos(an);
	
	double A = 0.25 * pow( 3.0 * cos(an) * cos(an) - 1, 2.0 );
	double B = 3.0 * pow( sin(an) * cos(an), 2.0 );
	double C = 0.75 * pow( sin(an), 4.0 );
	double J0 = 
	  A * fnJ( 0.0, ta ) + 
	  B * fnJ( 0.0, tb ) + 
	  C * fnJ( 0.0, tc );
	double J1 = 
	  A * fnJ( 1.0 * param.freq1, ta ) +
	  B * fnJ( 1.0 * param.freq1, tb ) +
	  C * fnJ( 1.0 * param.freq1, tc );
	double J2 = 
	  A * fnJ( 2.0 * param.freq1, ta ) + 
	  B * fnJ( 2.0 * param.freq1, tb ) +
	  C * fnJ( 2.0 * param.freq1, tc );
	
	XR[i + j * prot.H.size()] = 
	  CON_HH / pow( r, 6 ) * ( 6.0 * J2 - J0 );
	XRh[i + j * prot.H.size()] = 
	  CON_HH / pow( r, 6 ) * ( J0 + 3.0 * J1 + 6.0 * J2 );

	if( std::isnan( XR[i + j * prot.H.size()] ) ){
	  cout 
	    << "NAN in calculation: non-methyl or methyl on same hetero atom: "
	    << endl;
	}

	if( std::isnan( XRh[i + j * prot.H.size()] ) ){
	  cout 
	    << "NAN in calculation: non-methyl or methyl on same hetero atom: "
	    << endl;
	}
      }
      
      // i ^= methyl proton, j ^= normal proton
      if( ( prot.H[i]->hetatom->prot.size() == 3 && 
	    prot.H[j]->hetatom->prot.size() != 3 ) ){
	complex<double> cJ0( 0.0, 0.0 );
	complex<double> cJ1p( 0.0, 0.0 );
	complex<double> cJ1m( 0.0, 0.0 );
	complex<double> cJ2p( 0.0, 0.0 );
	complex<double> cJ2m( 0.0, 0.0 );
	
	for( int k = 0; k < 3; k++ ){
	  double dx = prot.H[j]->x - prot.H[i]->hetatom->prot[k]->x;
	  double dy = prot.H[j]->y - prot.H[i]->hetatom->prot[k]->y;
	  double dz = prot.H[j]->z - prot.H[i]->hetatom->prot[k]->z;
	  double r = sqrt( dx*dx + dy*dy + dz*dz );
	  
	  double ceta = ( dxs*dx + dys*dy + dzs*dz ) / r;
	  if( fabs( ceta ) > ONEM ){
	    ceta = ONEM;
	  }
	  ceta = acos( ceta );

	  double vp = dx * dxp + dy * dyp + dz * dzp;
	  double vq = dx * dxq + dy * dyq + dz * dzq;
	  double phi = vp / sqrt( vp*vp + vq*vq );
	  if( fabs( phi ) > ONEM ){
	    phi = ONEM;
	  }
	  phi = acos( phi );

	  cJ0 = cJ0 + ( 3.0 * pow( cos( ceta ), 2.0 ) - 1.0 ) / pow( r, 3.0 ); 
	  cJ1p = cJ1p + sin( ceta ) * cos( ceta ) * 
	    complex<double>( cos( phi ), sin( phi ) ) / pow( r, 3.0 );
	  cJ1m = cJ1m + sin( ceta ) * cos( ceta ) *
	    complex<double>( cos( phi ), -sin( phi ) ) / pow( r, 3.0 );
	  cJ2p = cJ2p + pow( sin( ceta ), 2.0 ) *
	    pow( complex<double>( cos( phi ), sin( phi ) ), 2.0 ) / 
	    pow( r, 3.0 );
	  cJ2m = cJ2m + pow( sin( ceta ), 2.0 ) *
	    pow( complex<double>( cos( phi ), -sin( phi )), 2.0 ) /
	    pow( r, 3.0 );
	}

	double A = 0.25 * abs( cJ0 ) * abs( cJ0 ) / pow( 3.0, 2.0 );
	double B = 1.5 * ( pow( abs( cJ1p ), 2.0 ) + 
			   pow( abs( cJ1m ), 2.0 ) ) / pow( 3.0, 2.0 );
	double C = 0.375 * ( pow( abs( cJ2p ), 2.0 ) + 
			     pow( abs( cJ2m ), 2.0 ) ) / pow( 3.0, 2.0 );
	double J0 = 
	  A * fnJ( 0.0, ta ) + 
	  B * fnJ( 0.0, tb ) + 
	  C * fnJ( 0.0, tc );
	double J1 = 
	  A * fnJ( 1.0 * param.freq1, ta ) + 
	  B * fnJ( 1.0 * param.freq1, tb ) + 
	  C * fnJ( 1.0 * param.freq1, tc );
	double J2 = 
	  A * fnJ( 2.0 * param.freq1, ta ) + 
	  B * fnJ( 2.0 * param.freq1, tb ) + 
	  C * fnJ( 2.0 * param.freq1, tc );
	
	XR[i + j * prot.H.size()] = CON_HH * ( 6.0 * J2 - J0 );
	XRh[i + j * prot.H.size()] = CON_HH * ( J0 + 3.0 * J1 + 6.0 * J2 );

	if( std::isnan( XR[i + j * prot.H.size()] ) ){
	  cout 
	    << "NAN in calculation: i ^= methyl proton, j ^= normal proton: "
	    << endl;
	}

	if( std::isnan( XRh[i + j * prot.H.size()] ) ){
	  cout 
	    << "NAN in calculation: i ^= methyl proton, j ^= normal proton: "
	    << endl;
	}
      }

      // i ^= normal proton, j ^= methyl proton
      if( ( prot.H[i]->hetatom->prot.size() != 3 && 
	    prot.H[j]->hetatom->prot.size() == 3 ) ){
	complex<double> cJ0( 0.0, 0.0 );
	complex<double> cJ1p( 0.0, 0.0 );
	complex<double> cJ1m( 0.0, 0.0 );
	complex<double> cJ2p( 0.0, 0.0 );
	complex<double> cJ2m( 0.0, 0.0 );
	
	for( int k = 0; k < 3; k++ ){
	  double dx = prot.H[i]->x - prot.H[j]->hetatom->prot[k]->x;
	  double dy = prot.H[i]->y - prot.H[j]->hetatom->prot[k]->y;
	  double dz = prot.H[i]->z - prot.H[j]->hetatom->prot[k]->z;
	  double r = sqrt( dx*dx + dy*dy + dz*dz );
	  
	  double ceta = ( dxs*dx + dys*dy + dzs*dz ) / r;
	  if( fabs( ceta ) > ONEM ){
	    ceta = ONEM;
	  }
	  ceta = acos( ceta );

	  double vp = dx * dxp + dy * dyp + dz * dzp;
	  double vq = dx * dxq + dy * dyq + dz * dzq;
	  double phi = vp / sqrt( vp*vp + vq*vq );
	  if( fabs( phi ) > ONEM ){
	    phi = ONEM;
	  }
	  phi = acos( phi );

	  cJ0 = cJ0 + ( 3.0 * pow( cos( ceta ), 2.0 ) - 1.0 ) / pow( r, 3.0 ); 
	  cJ1p = cJ1p + sin( ceta ) * cos( ceta ) * 
	    complex<double>( cos( phi ), sin( phi ) ) / pow( r, 3.0 );
	  cJ1m = cJ1m + sin( ceta ) * cos( ceta ) *
	    complex<double>( cos( phi ), -sin( phi ) ) / pow( r, 3.0 );
	  cJ2p = cJ2p + pow( sin( ceta ), 2.0 ) *
	    pow( complex<double>( cos( phi ), sin( phi ) ), 2.0 ) / 
	    pow( r, 3.0 );
	  cJ2m = cJ2m + pow( sin( ceta ), 2.0 ) *
	    pow( complex<double>( cos( phi ), -sin( phi )), 2.0 ) /
	    pow( r, 3.0 );
	}

	double A = 0.25 * abs( cJ0 ) * abs( cJ0 ) / pow( 3.0, 2.0 );
	double B = 1.5 * ( pow( abs( cJ1p ), 2.0 ) + 
			   pow( abs( cJ1m ), 2.0 ) ) / pow( 3.0, 2.0 );
	double C = 0.375 * ( pow( abs( cJ2p ), 2.0 ) + 
			     pow( abs( cJ2m ), 2.0 ) ) / pow( 3.0, 2.0 );
	double J0 = 
	  A * fnJ( 0.0, ta ) + 
	  B * fnJ( 0.0, tb ) + 
	  C * fnJ( 0.0, tc );
	double J1 = 
	  A * fnJ( 1.0 * param.freq1, ta ) + 
	  B * fnJ( 1.0 * param.freq1, tb ) + 
	  C * fnJ( 1.0 * param.freq1, tc );
	double J2 = 
	  A * fnJ( 2.0 * param.freq1, ta ) + 
	  B * fnJ( 2.0 * param.freq1, tb ) + 
	  C * fnJ( 2.0 * param.freq1, tc );
	
	XR[i + j * prot.H.size()] = CON_HH * ( 6.0 * J2 - J0 );
	XRh[i + j * prot.H.size()] = CON_HH * ( J0 + 3.0 * J1 + 6.0 * J2 );

	if( std::isnan( XR[i + j * prot.H.size()] ) ){
	  cout 
	    << "NAN in calculation: i ^= normal proton, j ^= methyl proton: "
	    << endl;
	}

	if( std::isnan( XRh[i + j * prot.H.size()] ) ){
	  cout 
	    << "NAN in calculation: i ^= normal proton, j ^= methyl proton: "
	    << endl;
	}
      }

      // methyl protons on different hetero atoms
      if( prot.H[i]->hetatom->prot.size() == 3 &&
	  prot.H[j]->hetatom->prot.size() == 3 &&
	  prot.H[i]->hetatom != prot.H[j]->hetatom ){
	complex<double> cJ0( 0.0, 0.0 );
	complex<double> cJ1p( 0.0, 0.0 );
	complex<double> cJ1m( 0.0, 0.0 );
	complex<double> cJ2p( 0.0, 0.0 );
	complex<double> cJ2m( 0.0, 0.0 );

	hetero* het1 = prot.H[i]->hetatom;
	hetero* het2 = prot.H[j]->hetatom;

	for( int k = 0; k < 3; k++ ){
	  for( int l = 0; l < 3; l++ ){
	    
	    double dx = het1->prot[k]->x - het2->prot[l]->x;
	    double dy = het1->prot[k]->y - het2->prot[l]->y;
	    double dz = het1->prot[k]->z - het2->prot[l]->z;
	    double r = sqrt( dx*dx + dy*dy + dz*dz );

	    double ceta = ( dxs*dx + dys*dy + dzs*dz ) / r;
	    if( fabs( ceta ) > ONEM ){
	      ceta = ONEM;
	    }
	    ceta = acos( ceta );

	    double vp = dx * dxp + dy * dyp + dz * dzp;
	    double vq = dx * dxq + dy * dyq + dz * dzq;
	    double phi = vp / sqrt( vp*vp + vq*vq );
	    if( fabs( phi ) > ONEM ){
	      phi = ONEM;
	    }
	    phi = acos( phi );

	    cJ0 = cJ0 + ( 3.0 * pow( cos( ceta ), 2.0 ) - 1.0 ) / 
	      pow( r, 3.0 ); 
	    cJ1p = cJ1p + sin( ceta ) * cos( ceta ) * 
	      complex<double>( cos( phi ), sin( phi ) ) / pow( r, 3.0 );
	    cJ1m = cJ1m + sin( ceta ) * cos( ceta ) *
	      complex<double>( cos( phi ), -sin( phi ) ) / pow( r, 3.0 );
	    cJ2p = cJ2p + pow( sin( ceta ), 2.0 ) *
	      pow( complex<double>( cos( phi ), sin( phi ) ), 2.0 ) / 
	      pow( r, 3.0 );
	    cJ2m = cJ2m + pow( sin( ceta ), 2.0 ) *
	      pow( complex<double>( cos( phi ), -sin( phi )), 2.0 ) /
	      pow( r, 3.0 );
	  }
	}

	double A = 0.25 * abs( cJ0 ) * abs( cJ0 ) / 
	  ( pow( 3.0, 2.0 ) * pow( 3.0, 2.0 ) );
	double B = 1.5 * ( pow( abs( cJ1p ), 2.0 ) + 
			   pow( abs( cJ1m ), 2.0 ) ) / 
	  ( pow( 3.0, 2.0 ) * pow( 3.0, 2.0 ) );
	double C = 0.375 * ( pow( abs( cJ2p ), 2.0 ) + 
			     pow( abs( cJ2m ), 2.0 ) ) / 
	  ( pow( 3.0, 2.0 ) * pow( 3.0, 2.0 ) );

	if( std::isnan( A ) ){
	  cout << "std::isnan( A ) == true" << endl;
	}

	if( std::isnan( B ) ){
	  cout << "std::isnan( B ) == true" << endl;
	}
	
	if( std::isnan( C ) ){
	  cout << "std::isnan( C ) == true" << endl;
	}

	double J0 = 
	  A * fnJ( 0.0, ta ) + 
	  B * fnJ( 0.0, tb ) + 
	  C * fnJ( 0.0, tc );
	double J1 = 
	  A * fnJ( 1.0 * param.freq1, ta ) + 
	  B * fnJ( 1.0 * param.freq1, tb ) + 
	  C * fnJ( 1.0 * param.freq1, tc );
	double J2 = 
	  A * fnJ( 2.0 * param.freq1, ta ) + 
	  B * fnJ( 2.0 * param.freq1, tb ) + 
	  C * fnJ( 2.0 * param.freq1, tc );
	
	XR[i + j * prot.H.size()] = CON_HH * ( 6.0 * J2 - J0 );
	XRh[i + j * prot.H.size()] = CON_HH * ( J0 + 3.0 * J1 + 6.0 * J2 );

	if( std::isnan( XR[i + j * prot.H.size()] ) ){
	  cout 
	    << "NAN in calculation: methyl protons on different hetero atoms: "
	    << "XR" << endl;
	}

	if( std::isnan( XRh[i + j * prot.H.size()] ) ){
	  cout 
	    << "NAN in calculation: methyl protons on different hetero atoms: "
	    << "XRh" << endl;
	}
      }
    }
  }

  cout << "Relax matric set up." << endl;
  
  // Mirror Matrix to the other side
  for( unsigned int i = 0; i < prot.H.size(); i++ ){
    for( unsigned int j = i+1; j < prot.H.size(); j++ ){
      XR[j + i * prot.H.size()] = XR[i + j * prot.H.size()];
      XRh[j + i * prot.H.size()] = XRh[i + j * prot.H.size()];
    }
  }
  // set diagonal to 0.0
  for( unsigned int i = 0; i < prot.H.size(); i++ ){
    XR[i + i * prot.H.size()] = 0.0;
  }
  // add the ith row of XRh but without the diagonal element to
  // the diagonal element of XR
  for( unsigned int i = 0; i < prot.H.size(); i++ ){
    for( unsigned int k = 0; k < prot.H.size(); k++ ){
      if( i != k ){
	XR[i + i * prot.H.size()] = 
	  XR[i + i * prot.H.size()] + XRh[i + k * prot.H.size()];
      }
    }
  }

  double fc = param.freq1 * 0.251506;
  double fz = param.freq1 - fc;
  double fd = param.freq1 + fc;

  // store proton-carbon interactions in extra matrix (XRc)
  
  for( unsigned int i = 0; i < prot.H.size(); i++ ){
    for( unsigned int j = 0; j < prot.C13.size(); j++ ){
      double dxc = prot.H[i]->x - prot.C13[j]->x;
      double dyc = prot.H[i]->y - prot.C13[j]->y;
      double dzc = prot.H[i]->z - prot.C13[j]->z;
      double rc = sqrt( dxc*dxc + dyc*dyc + dzc*dzc );
      
      if( rc < HCMAX ){
	double an = ( dxs*dxc + dys*dyc + dzs*dzc ) / ( ds*rc );
	if( fabs( an ) > ONEM ){
	  an = ONEM;
	}
	an = acos( an );
	
	double A = 0.25 * pow( 3.0 * cos(an) * cos(an) - 1, 2.0 );
	double B = 3.0 * pow( sin(an) * cos(an), 2.0 );
	double C = 0.75 * pow( sin(an), 4.0 );
	double J0c = A * fnJ( fz, ta ) + B * fnJ( fz, tb ) + C * fnJ( fz, tc );
	double J1c = 
	  A * fnJ( param.freq1, ta ) + 
	  B * fnJ( param.freq1, tb ) + 
	  C * fnJ( param.freq1, tc );
	double J2c = A * fnJ( fd, ta ) + B * fnJ( fd, tb ) + C * fnJ( fd, tc );
	
	XRc[i + j * prot.H.size()] = 
	  CON_HC / pow( rc, 6.0 ) * ( J0c + 3.0 * J1c + 6.0 * J2c );
      }
      else{
	XRc[i + j * prot.H.size()] = 0.0;
      }
    }
  }

  double fn = param.freq1 * 0.101376;
  fz = param.freq1 - fn;
  fd = param.freq1 + fn;

  // store proton-nitrogen interactions in extra matrix (XRn)
  for( unsigned int i = 0; i < prot.H.size(); i++ ){
    for( unsigned int j = 0; j < prot.N15.size(); j++ ){
      double dxn = prot.H[i]->x - prot.N15[j]->x;
      double dyn = prot.H[i]->y - prot.N15[j]->y;
      double dzn = prot.H[i]->z - prot.N15[j]->z;
      double rn = sqrt( dxn*dxn + dyn*dyn + dzn*dzn );
      
      if( rn < HNMAX ){
	double an = ( dxs*dxn + dys*dyn + dzs*dzn ) / ( ds*rn );
	if( fabs( an ) > ONEM ){
	  an = ONEM;
	}
	an = acos( an );

	double A = 0.25 * pow( 3.0 * cos( an ) * cos( an ) - 1, 2.0 );
	double B = 3.0 * pow( sin( an ) * cos( an ), 2.0 );
	double C = 0.75 * pow( sin( an ), 4.0 );
	double J0n = A * fnJ( fz, ta ) + B * fnJ( fz, tb ) + C * fnJ( fz, tc );
	double J1n = 
	  A * fnJ( param.freq1, ta ) + 
	  B * fnJ( param.freq1, tb ) + 
	  C * fnJ( param.freq1, tc );
	double J2n = A * fnJ( fd, ta ) + B * fnJ( fd, tb ) + C * fnJ( fd, tc );

	XRn[i + j * prot.H.size()] = 
	  CON_HN / pow( rn, 6.0 ) * ( J0n + 3.0 * J1n + 6.0 * J2n );
      }
      else{
	XRn[i + j * prot.H.size()] = 0.0;
      }
    }
  }

  // merge the extra matrices with the XR matrix
  // case: C13
  if( param.label == "C13" ){
    for( unsigned int i = 0; i < prot.H.size(); i++ ){
      for( unsigned int k = 0; k < prot.C13.size(); k++ ){
	XR[i + i * prot.H.size()] = 
	  XR[i + i * prot.H.size()] + XRc[i + k * prot.H.size()];
      }
    }
  }

  // case: N15
  if( param.label == "N15" ){
    for( unsigned int i = 0; i < prot.H.size(); i++ ){
      for( unsigned int k = 0; k < prot.N15.size(); k++ ){
	XR[i + i * prot.H.size()] = 
	  XR[i + i * prot.H.size()] + XRn[i + k * prot.H.size()];
      }
    }
  }

  // case: double
  if( param.label == "double" ){
    for( unsigned int i = 0; i < prot.H.size(); i++ ){
      for( unsigned int k = 0; k < prot.C13.size(); k++ ){
	XR[i + i * prot.H.size()] =
	  XR[i + i * prot.H.size()] + XRc[i + k * prot.H.size()];
      }
      for( unsigned int k = 0; k < prot.N15.size(); k++ ){
	XR[i + i * prot.H.size()] =
	  XR[i + i * prot.H.size()] + XRn[i + k * prot.H.size()];
      }
    }
  }

  // add the Rext elements to the diagonal of XR, Rext ^= Leak parameter of 
  // the ppm file

  for( unsigned int i = 0; i < prot.H.size(); i++ ){
    XR[i + i * prot.H.size()] = 
      XR[i + i * prot.H.size()] + prot.H[i]->leak;
  }
  
  // TODO, WARNING: here we leave out the averaging step to do the proper
  // weighting of PDB populations, currently we use only one structure

  // Eigenvalue computation
  gsl_matrix* A = gsl_matrix_calloc( prot.H.size(), prot.H.size() );

  for( unsigned int i = 0; i < prot.H.size(); i++ ){
    for( unsigned int j = 0; j < prot.H.size(); j++ ){
      gsl_matrix_set( A, i, j, XR[i + j * prot.H.size()] );
    }
  }

  gsl_vector* eval = gsl_vector_calloc( prot.H.size() );
  gsl_matrix* evec = gsl_matrix_calloc( prot.H.size(), prot.H.size() );

  gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc( prot.H.size() );

  cout << "Starting Eigenvalue computation." << endl;
  gsl_eigen_symmv( A, eval, evec, w );
  cout << "Eigenvalue computation done." << endl;

  gsl_eigen_symmv_free( w );
  
  /*
    result = gsl_eigen_symmv_sort( eval, evec, GSL_EIGEN_SORT_ABS_ASC );
  */

  double vector_length;
  for( unsigned int i = 0; i < prot.H.size(); i++ ){
    vector_length = 0.0;
    for( unsigned int j = 0; j < prot.H.size(); j++ ){
      vector_length += gsl_matrix_get( evec, j, i ) * 
	gsl_matrix_get( evec, j, i );
    }
  }

  for( unsigned int i = 0; i < prot.H.size(); i++ ){
    for( unsigned int j = 0; j < prot.H.size(); j++ ){
      XR[i + j * prot.H.size()] = gsl_matrix_get( evec, i, j );
    }
  }

  for( unsigned int i = 0; i < prot.H.size(); i++ ){
    RD[i] = exp( -gsl_vector_get( eval, i) * param.Trex );
  }

  // initialize RG with 0.0
  for( unsigned int i = 0; i < prot.H.size(); i++ ){
    for( unsigned int j = 0; j < prot.H.size(); j++ ){
      RG[i] = 0.0;
    }
  }

  // calculate a triangular matrix RG from XR and RD
  unsigned int maxHarray = prot.H.size();
  for( unsigned int i = 0; i < maxHarray; i++ ){
    for( unsigned int j = i+1; j < maxHarray; j++ ){
      for( unsigned int k = 0; k < maxHarray; k++ ){
	RG[i + j * maxHarray] += 
	  XR[i + k * maxHarray] * XR[j + k * maxHarray] *
	  RD[k];
      }
    }
  }
  
  // mirror the triangular matrix RG
  for( unsigned int i = 0; i < prot.H.size(); i++ ){
    for( unsigned int j = i+1; j < prot.H.size(); j++ ){
      RG[j + i * prot.H.size()] = RG[i + j * prot.H.size()];
    }
  }

  // calculate the diagonal elements of RG, but why isn't this step integrated
  // into the calculation of the triangular matrix above?
  // MEMO: Check for integration possibility
  for( unsigned int i = 0; i < prot.H.size(); i++ ){
    for( unsigned int k = 0; k < prot.H.size(); k++ ){
      RG[i + i * prot.H.size()] += 
	XR[i + k * prot.H.size()] * XR[i + k * prot.H.size()] * 
	RD[k];
    }
  }

  for( unsigned int i = 0; i < prot.H.size(); i++ ){
    P[i] = 0.0;
  }

  // sum up columns of RG and store them in vector P, although
  // modified: P -> 1-P, P corresponds to recovery percentage
  for( unsigned int i = 0; i < prot.H.size(); i++ ){
    for( unsigned int k = 0; k < prot.H.size(); k++ ){
      P[i] += RG[i + k * prot.H.size()];
    }
    P[i] = 1.0 - P[i];
  }

  // And here the differences between single heteroatom mode and
  // double heteroatom mode start

  // Write out transfer.dat and recover.dat
  vector<proton*>* editprot = NULL;

  if( param.heteromode == "single" ){
    if( param.heterobound == "C13" ){
      editprot = &prot.HC;
    }
    if( param.heterobound == "N15" ){
      editprot = &prot.HN;
    }
  }
  vector<proton*>* hetbound = NULL;
  vector<proton*>* hetfree = NULL;
  if( param.heteromode == "double" ){
    if( param.heterobound == "C13" ){
      hetbound = &prot.HC;
    }
    if( param.heterobound == "N15" ){
      hetbound = &prot.HN;
    }
    if( param.heterofree == "C13" ){
      hetfree = &prot.HC;
    }
    if( param.heterofree == "N15" ){
      hetfree = &prot.HN;
    }
  }

  ofstream output;

  // The output of recover.dat and transfer.dat is missing here
  // TODO: reimplement

  // initialize G with 0.0
  for( unsigned int i = 0; i < prot.H.size(); i++ ){
    for( unsigned int j = i; j < prot.H.size(); j++ ){
      G[i + j * prot.H.size()] = 0.0;
    }
  }

  double* noe = NULL;
  double* doublenoe = NULL;
  if( param.heteromode == "single" ){
    noe = (double*) calloc( editprot->size() * prot.H.size(),
			    sizeof(double) );
  }
  if( param.heteromode == "double" ){
    doublenoe = (double*) calloc( hetbound->size() * hetfree->size(),
				  sizeof(double) );
  }

  for( unsigned int i = 0; i < prot.H.size(); i++ ){
    DD[i] = exp( -gsl_vector_get( eval,i ) * param.mixt[0] );
  }

  // calculate triangular matrix G without diagonal
  for( unsigned int i = 0; i < prot.H.size(); i++ ){
    for( unsigned int j = i+1; j < prot.H.size(); j++ ){
      for( unsigned int k = 0; k < prot.H.size(); k++ ){
	G[i + j * prot.H.size()] += 
	  XR[i + k * prot.H.size()] * XR[j + k * prot.H.size()] * DD[k];
      }
    }
  }

  // mirror triangular matrix G 
  for( unsigned int i = 0; i < prot.H.size(); i++ ){
    for( unsigned int j = i+1; j < prot.H.size(); j++ ){
      G[j + i * prot.H.size()] = G[i + j * prot.H.size()];
    }
  }

  // calculate diagonal of matrix G, again: Why isn't this integrated in
  // the calculation step of the triangular matrix. 
  // MEMO: Check for the possibility of integration
  for( unsigned int i = 0; i < prot.H.size(); i++ ){
    for( unsigned int k = 0; k < prot.H.size(); k++ ){
      G[i + i * prot.H.size()] += 
	XR[i + k * prot.H.size()] * XR[i + k * prot.H.size()] * DD[k];
    }
  }

  // HSQC step?

  if( param.heteromode == "single" ){
    for( unsigned int i = 0; i < editprot->size(); i++ ){
      for( unsigned int j = 0; j < prot.H.size(); j++ ){
	for( unsigned int k = 0; k < prot.H.size(); k++ ){
	  if( (*editprot)[i] == prot.H[k] ){
	    noe[i + j*editprot->size()] = P[j] * G[j + k * prot.H.size()] *
	      (*editprot)[i]->hetatom->trans * 
	      (*editprot)[i]->hetatom->fold3 * 
	      prot.H[j]->fold2;
	  }
	}
      }
    }
  }

  cout << "Starting formerly very inefficient ~O(n^4) calculation." << endl;

  if( param.heteromode == "double" ){
    for( unsigned int i = 0; i < hetfree->size(); i++ ){
      for( unsigned int k = 0; k < prot.H.size(); k++ ){
	if( (*hetfree)[i] == prot.H[k] ){
	  for( unsigned int j = 0; j < hetbound->size(); j++ ){
	    for( unsigned int l = 0; l < prot.H.size(); l++ ){
	      if( (*hetbound)[j] == prot.H[l] ) {
		doublenoe[i * hetbound->size() + j] = 
		  P[l] * G[l + k * prot.H.size()] * 
		  (*hetfree)[i]->hetatom->trans *
		  (*hetbound)[j]->hetatom->trans *
		  (*hetfree)[i]->hetatom->fold3 *
		  (*hetbound)[j]->hetatom->fold3 *
		  (*hetbound)[j]->fold2; // TODO: Is this necessary?
	      }
	    }
	  }
	}
      }
    }
  }

  cout << "Formerly inefficient calculation done." << endl;
    
  gsl_matrix_free( evec );
  gsl_vector_free( eval );
  gsl_matrix_free( A );

  /* Write out results */

  output.open( "output.txt" );
  output.precision(5);
  if( !output.good() ){
    cout << "Can't open output.txt, exiting." << endl;
    exit(1);
  }

  output << param.Nmix << " ";
  for( unsigned int i = 0; i < param.Nmix; i++ ){
    output << param.mixt[i] << " ";
  }
  output << endl;
  output << param.size1 << " " << param.size2 << " " << param.size3 << endl;
  output << param.sweepwidth1 << " " << param.sweepwidth2 << " " 
	 << param.sweepwidth3 << endl;

  if( param.heteromode == "single" ){
    for( unsigned int i = 0; i < editprot->size(); i++ ){
      for( unsigned int j = 0; j < prot.H.size(); j++ ){
	double value = param.noemin / ( (*editprot)[i]->hetatom->valveg *
					prot.H[j]->valvef );
	if( fabs( noe[i + j * editprot->size()] ) > value ){
	  output << (*editprot)[i]->hetatom->residue->residueNumber << " "
		 << (*editprot)[i]->Name << " "
		 << prot.H[j]->hetatom->residue->residueNumber << " "
		 << prot.H[j]->Name << " "
		 << (*editprot)[i]->hetatom->residue->residueNumber << " "
		 << (*editprot)[i]->hetatom->Name << " "
		 << (*editprot)[i]->spec_shift1 << " "
		 << prot.H[j]->spec_shift2 << " "
		 << (*editprot)[i]->hetatom->spec_shift3 << " "
		 << (*editprot)[i]->lw1 << " "
		 << prot.H[j]->lw2 << " "
		 << (*editprot)[i]->hetatom->lw3 << " "
		 << endl;
	  output << noe[i + j * editprot->size()] << endl;
	}
      }
    }
  }

  ofstream shiftoutput( "shiftoutput.txt" );
  if( !shiftoutput.good() ){
    cout << "Can't open shiftoutput.txt, exiting." << endl;
    exit(1);
  }

  if( param.heteromode == "double" ){
    for( unsigned int i = 0; i < hetfree->size(); i++ ){
      for( unsigned int j = 0; j < hetbound->size(); j++ ){
	double value = param.noemin / ( (*hetfree)[i]->hetatom->valveg *
					(*hetbound)[j]->valvef );
	if( fabs( doublenoe[i * hetbound->size() + j] ) > value ){
	  output << (*hetbound)[j]->hetatom->residue->residueNumber << " "
		 << (*hetbound)[j]->Name << " "
		 << (*hetfree)[i]->hetatom->residue->residueNumber << " "
		 << (*hetfree)[i]->hetatom->Name << " "
		 << (*hetbound)[j]->hetatom->residue->residueNumber << " "
		 << (*hetbound)[j]->hetatom->Name << " "
		 << (*hetbound)[j]->spec_shift1 << " "
		 << (*hetfree)[i]->hetatom->spec_shift2 << " "
		 << (*hetbound)[j]->hetatom->spec_shift3 << " "
		 << (*hetbound)[j]->lw1 << " "
		 << (*hetfree)[i]->hetatom->lw2 << " "
		 << (*hetbound)[j]->hetatom->lw3 
		 << endl;
	  output << doublenoe[i * hetbound->size() + j] << endl;
	  shiftoutput 
	    << (*hetbound)[j]->hetatom->residue->residueNumber << " "
	    << (*hetbound)[j]->Name << " "
	    << (*hetfree)[i]->hetatom->residue->residueNumber << " "
	    << (*hetfree)[i]->hetatom->Name << " "
	    << (*hetbound)[j]->hetatom->residue->residueNumber << " "
	    << (*hetbound)[j]->hetatom->Name << " "
	    << (*hetbound)[j]->shift << " "
	    << (*hetfree)[i]->hetatom->shift << " "
	    << (*hetbound)[j]->hetatom->shift << " "
	    << (*hetbound)[j]->lw1 << " "
	    << (*hetfree)[i]->hetatom->lw2 << " "
	    << (*hetbound)[j]->hetatom->lw3 
	    << endl;
	  shiftoutput << doublenoe[i * hetbound->size() + j] << endl;
	}
	else{
	  shiftoutput << "Filtered: "
	    << (*hetbound)[j]->hetatom->residue->residueNumber << " "
	    << (*hetbound)[j]->Name << " "
	    << (*hetfree)[i]->hetatom->residue->residueNumber << " "
	    << (*hetfree)[i]->hetatom->Name << " "
	    << (*hetbound)[j]->hetatom->residue->residueNumber << " "
	    << (*hetbound)[j]->hetatom->Name << " "
	    << (*hetbound)[j]->shift << " "
	    << (*hetfree)[i]->hetatom->shift << " "
	    << (*hetbound)[j]->hetatom->shift << " "
	    << (*hetbound)[j]->lw1 << " "
	    << (*hetfree)[i]->hetatom->lw2 << " "
	    << (*hetbound)[j]->hetatom->lw3 
	    << endl;
	  shiftoutput << doublenoe[i * hetbound->size() + j] << endl;
	}
      }
    }
  }
  
  shiftoutput.close();
  output.close();

  /* free memory */
  free( doublenoe );
  free( noe );
  free( DD );
  free( G );
  free( P );
  free( RG );
  free( RD );
  free( XRn );
  free( XRc );
  free( XRh );
  free( XR );
  return 0;
}
