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
#include <string.h>
#include <math.h>
#include "NMRparseparameters.h"

parameters::parameters(){
  freq1 = NAN;
  freq2 = NAN;
  freq3 = NAN;
  sweepwidth1 = NAN;
  sweepwidth2 = NAN;
  sweepwidth3 = NAN;
  size1 = -1;
  size2 = -1;
  size3 = -1;
  refpoint1 = NAN;
  refpoint2 = NAN;
  refpoint3 = NAN;
  refppm1 = NAN;
  refppm2 = NAN;
  refppm3 = NAN;
  linewidthfactor1 = NAN;
  linewidthfactor2 = NAN;
  linewidthfactor3 = NAN;
  pdbfile = "";
  shiftfile = "";
  sequence = "";
  tumbling = "";
  isocorrel = NAN;
  label = "";
  Trex = NAN;
  Nmix = 0;
  mixt = NULL;
  heteromode = "";
  heterobound = "";
  heterofree = "";
  inept_hetbound = NAN;
  inept_hetbound_rev = NAN;
  inept_hetbound_MQ = NAN;
  inept_hetfree = NAN;
  inept_hetfree_rev = NAN;
  inept_hetfree_MQ = NAN;
  delay_proton_set = false;
  delay_hetbound_set = false;
  delay_hetfree_set = false;
  hetbound_sensitivity_enhanced_set = false;
  hetfree_sensitivity_enhanced_set = false;
  cors2 = NAN;
  cors3 = NAN;
  noemin = NAN;
  symatom1_resID = -1;
  symatom2_resID = -1;
  symatom1_Name = "";
  symatom2_Name = "";
  seq_start = -1;
}

struct parameters NMRparseparameters( const char* filename ){
  struct parameters data;
  ifstream parafile;
  parafile.open( filename, ios_base::in );
  if( !parafile ){
    cout << "Can't open" << filename << endl;
  }

  char buffer[1024];
  string bufstring;
  bool finished = false;
  string errorstring;
  while( parafile.eof() != true ){
    finished = true;
    errorstring = "";
    parafile.getline( buffer, 1024 );
    bufstring = buffer;
    stringstream parsebuffer( bufstring, ios_base::in );
    string keyword;
    parsebuffer >> keyword;

    if( keyword == "#freq:" ){
      parsebuffer >> data.freq1 >> data.freq2 >> data.freq3;
    }
    if( isnan( data.freq1 ) || isnan( data.freq2 ) || isnan( data.freq3 ) ){
      errorstring.append( "#freq parameters missing\n" );
      finished = false;
    }

    if( keyword == "#spwd:" ){
      parsebuffer >> data.sweepwidth1 >> data.sweepwidth2 >> data.sweepwidth3;
    }
    if( isnan( data.sweepwidth1 ) || 
	isnan( data.sweepwidth2 ) ||
	isnan( data.sweepwidth3 ) ){
      errorstring.append( "#spwd: parameters missing\n" );
      finished = false;
    }

    if( keyword == "#size:" ){
      parsebuffer >> data.size1 >> data.size2 >> data.size3;
    }
    if( data.size1 == -1 || data.size2 == -1 || data.size3 == -1 ){
      errorstring.append( "#size: parameters missing\n" );
      finished = false;
    }

    if( keyword == "#rfpt:" ){
      parsebuffer >> data.refpoint1 >> data.refpoint2 >> data.refpoint3;
    }
    if( isnan( data.refpoint1 ) || 
	isnan( data.refpoint2 ) ||
	isnan( data.refpoint3 ) ){
      errorstring.append( "#rfpt: parameters missing\n" );
      finished = false;
    }

    if( keyword == "#rfpm:" ){
      parsebuffer >> data.refppm1 >> data.refppm2 >> data.refppm3;
    }
    if( isnan( data.refppm1 ) ||
	isnan( data.refppm2 ) ||
	isnan( data.refppm3 ) ){
      errorstring.append( "#rfpm: parameters missing\n" );
      finished = false;
    }

    if( keyword == "#lwfc:" ){
      parsebuffer >> data.linewidthfactor1 
		  >> data.linewidthfactor2
		  >> data.linewidthfactor3;
    }
    if( isnan( data.linewidthfactor1 ) ||
	isnan( data.linewidthfactor2 ) ||
	isnan( data.linewidthfactor3 ) ){
      errorstring.append( "#lwfc: parameters missing\n" );
      finished = false;
    }

    if( keyword == "#seqfile:" ){
      parsebuffer >> data.sequence;
    }
    if( data.sequence.empty() ){
      errorstring.append( "#seqfile: parameter missing\n" );
      finished = false;
    }

    if( keyword == "#pdbfile:" ){
      parsebuffer >> data.pdbfile;
    }
    if( data.pdbfile.empty() ){
      errorstring.append( "#pdbfile: parameter missing\n" );
      finished = false;
    }

    if( keyword == "#shiftfile:" ){
      parsebuffer >> data.shiftfile;
    }
    if( data.shiftfile.empty() ){
      errorstring.append( "#shiftfile: parameter missing\n" );
      finished = false;
    }

    if( keyword == "#label:" ){
      parsebuffer >> data.label;
    }
    if( data.label.empty() ){
      errorstring.append( "#label: parameter missing\n" );
      finished = false;
    }

    if( keyword == "#tumbling:" ){
      parsebuffer >> data.tumbling;
    }
    if( data.tumbling.empty() ){
      errorstring.append( "#tumbling: parameter missing\n" );
      finished = false;
    }

    if( keyword == "#isocorrel_ns:" ){
      parsebuffer >> data.isocorrel;
    }
    if( isnan( data.isocorrel ) ){
      errorstring.append( "#isocorrel_ns: parameter missing\n" );
      finished = false;
    }

    if( keyword == "#relax_delay_sec:" ){
      parsebuffer >> data.Trex;
    }
    if( isnan( data.Trex ) ){
      errorstring.append( "#relax_delay_sec: parameter missing\n" );
      finished = false;
    }

    if( keyword == "#Mixing_times_sec:" ){
      parsebuffer >> data.Nmix;
      data.mixt = (double*) calloc( data.Nmix, sizeof(double) );
      for( unsigned int i = 0; i < data.Nmix; i++ ){
	parsebuffer >> data.mixt[i];
      }
    }
    if( data.Nmix == 0 || data.mixt == NULL ){
      errorstring.append( "#Mixing_times_sec: parameters missing\n" );
      finished = false;
    }

    if( keyword == "#heteromode:" ){
      parsebuffer >> data.heteromode;
    }
    if( data.heteromode.empty() ||
	( data.heteromode != "single" &&
	  data.heteromode != "double" ) ){
      errorstring.append( "#heteromode: parameter missing or != single || double\n" );
      finished = false;
    }
    
    if( keyword == "#heterobound:" ){
      parsebuffer >> data.heterobound;
    }
    if( data.heterobound.empty() ){
      errorstring.append( "#heterobound: parameter missing\n" );
      finished = false;
    }

    if( keyword == "#heterofree:" ){
      parsebuffer >> data.heterofree;
    }
    if( data.heteromode == "double" && data.heterofree.empty() ){
      errorstring.append( "#heterofree: must be set when #heteromode: == double\n" );
      finished = false;
    }

    if( keyword == "#inept_hetbound:" ){
      parsebuffer >> data.inept_hetbound;
    }
    if( isnan( data.inept_hetbound ) ){
      errorstring.append( "#inept_hetbound: parameter missing\n" );
      finished = false;
    }

    if( keyword == "#inept_hetbound_rev:" ){
      parsebuffer >> data.inept_hetbound_rev;
    }
    if( isnan( data.inept_hetbound_rev ) ){
      errorstring.append( "#inept_hetbound_rev: parameter missing\n" );
      finished = false;
    }

    if( keyword == "#inept_hetbound_MQ:" ){
      parsebuffer >> data.inept_hetbound_MQ;
    }
    if( data.hetbound_sensitivity_enhanced_set && 
	data.hetbound_sensitivity_enhanced &&
	isnan( data.inept_hetbound_MQ ) ){
      errorstring.append( "#inept_hetbound_MQ: parameter missing\n" );
      finished = false;
    }

    if( keyword == "#inept_hetfree:" ){
      parsebuffer >> data.inept_hetfree;
    }
    if( data.heteromode == "double" && isnan( data.inept_hetfree ) ){
      errorstring.append( "#inept_hetfree: parameter missing\n" );
      finished = false;
    }

    if( keyword == "#inept_hetfree_rev:" ){
      parsebuffer >> data.inept_hetfree_rev;
    }
    if( data.heteromode == "double" && isnan( data.inept_hetfree_rev ) ){
      errorstring.append( "#inept_hetfree_rev: parameter missing\n" );
      finished = false;
    }

    if( keyword == "#inept_hetfree_MQ:" ){
      parsebuffer >> data.inept_hetfree_MQ;
    }
    if( data.heteromode == "single" && 
	data.hetfree_sensitivity_enhanced_set && 
	data.hetfree_sensitivity_enhanced &&
	isnan( data.inept_hetfree_MQ ) ){
      errorstring.append( "#inept_hetfree_MQ: parameter missing\n" );
      finished = false;
    }

    if( keyword == "#delay_proton:" ){
      parsebuffer >> boolalpha >> data.delay_proton;
      data.delay_proton_set = true;
    }
    if( data.heteromode == "single" && !data.delay_proton_set ){
      errorstring.append( "#delay_proton: parameter missing\n" );
      finished = false;
    }

    if( keyword == "#delay_hetbound:" ){
      parsebuffer >> boolalpha >> data.delay_hetbound;
      data.delay_hetbound_set = true;
    }
    if( !data.delay_hetbound_set ){
      errorstring.append( "#delay_hetbound: parameter missing\n" );
      finished = false;
    }

    if( keyword == "#delay_hetfree:" ){
      parsebuffer >> boolalpha >> data.delay_hetfree;
      data.delay_hetfree_set = true;
    }
    if( data.heteromode == "double" && !data.delay_hetfree_set ){
      errorstring.append( "#delay_hetfree: parameter missing\n" );
      finished = false;
    }

    if( keyword == "#hetbound_sensitivity_enhanced:" ){
      parsebuffer >> boolalpha >> data.hetbound_sensitivity_enhanced;
      data.hetbound_sensitivity_enhanced_set = true;
    }
    if( !data.hetbound_sensitivity_enhanced_set ){
      errorstring.append( "#hetbound_sensitivity_enhanced: parameter missing\n" );
      finished = false;
    }

    if( keyword == "#hetfree_sensitivity_enhanced:" ){
      parsebuffer >> boolalpha >> data.hetfree_sensitivity_enhanced;
      data.hetfree_sensitivity_enhanced_set = true;
    }
    if( data.heteromode == "double" && 
	!data.hetfree_sensitivity_enhanced_set ){
      errorstring.append( "#hetfree_sensitivity_enhanced: parameter missing\n" );
      finished = false;
    }

    if( keyword == "#noemin:" ){
      parsebuffer >> data.noemin;
    }
    if( isnan( data.noemin ) ){
      errorstring.append( "#noemin: parameter missing\n" );
      finished = false;
    }

    if( keyword == "#symatom1_resID:" ){
      parsebuffer >> data.symatom1_resID;
    }
    if( keyword == "#symatom2_resID:" ){
      parsebuffer >> data.symatom2_resID;
    }

    if( keyword == "#symatom1_Name:" ){
      parsebuffer >> data.symatom1_Name;
    }
    if( keyword == "#symatom2_Name:" ){
      parsebuffer >> data.symatom2_Name;
    }

    if( keyword == "#seq_start:" ){
      parsebuffer >> data.seq_start;
    }
    if( data.seq_start == -1 ){
      errorstring.append( "#seq_start: parameter missing\n" );
      finished = false;
    }
  }
  parafile.close();

  if( !finished ){
    cout << "Parameters are missing in the parameter file, exiting." << endl;
    cout << "Details:" << endl;
    cout << errorstring;
    exit(1);
  }

  return data;
}
