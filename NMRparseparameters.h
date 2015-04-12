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

#ifndef __parseparameter
#define __parseparameter

#include <string>

using namespace std;

class parameters{
 public:
  double freq1;
  double freq2;
  double freq3;
  double sweepwidth1;
  double sweepwidth2;
  double sweepwidth3;
  int size1;
  int size2;
  int size3;
  double refpoint1;
  double refpoint2;
  double refpoint3;
  double refppm1;
  double refppm2;
  double refppm3;
  double linewidthfactor1;
  double linewidthfactor2;
  double linewidthfactor3;
  string pdbfile;
  string shiftfile;
  string sequence;
  string tumbling;
  double isocorrel;
  string label;
  double Trex;
  unsigned int Nmix;
  double* mixt;
  string heteromode;
  string heterobound;
  string heterofree;
  double inept_hetbound;
  double inept_hetbound_rev;
  double inept_hetbound_MQ;
  double inept_hetfree;
  double inept_hetfree_rev;
  double inept_hetfree_MQ;
  bool delay_proton;
  bool delay_proton_set;
  bool delay_hetbound;
  bool delay_hetbound_set;
  bool delay_hetfree;
  bool delay_hetfree_set;
  bool hetbound_sensitivity_enhanced;
  bool hetbound_sensitivity_enhanced_set;
  bool hetfree_sensitivity_enhanced;
  bool hetfree_sensitivity_enhanced_set;
  double cors2;
  double cors3;
  double noemin;
  int symatom1_resID;
  int symatom2_resID;
  string symatom1_Name;
  string symatom2_Name;
  int seq_start;
  parameters();
};

extern struct parameters NMRparseparameters( const char* filename );

#endif /* __parseparameter */
