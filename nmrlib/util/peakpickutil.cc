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
#include "../include/peakpickutil.h"

peakres::peakres(){
  intensity = NAN;
  accuracy = NAN;
  peakweight = NAN;
  pairweight = NAN;
}

void peakres::tozero(){
  intensity = 0.0;
  accuracy = 0.0;
  peakweight = 0.0;
  pairweight = 0.0;
}

void peakres::operator+=( peakres p ){
  if( peakweight + p.peakweight > 0.0 ){
    intensity = ( intensity * peakweight + p.intensity * p.peakweight ) /
      ( peakweight + p.peakweight );
    accuracy = ( accuracy * peakweight + p.accuracy * p.peakweight ) /
      ( peakweight + p.peakweight );
  }
  peakweight += p.peakweight;
  pairweight += p.pairweight;
}

peakres peakres::operator+ (peakres p){
  peakres temp;
  if( peakweight + p.peakweight > 0.0 ){
    temp.intensity = ( intensity * peakweight + p.intensity * p.peakweight ) /
      ( peakweight + p.peakweight );
    temp.accuracy = ( accuracy * peakweight + p.accuracy * p.peakweight ) /
      ( peakweight + p.peakweight );
  }
  else{
    temp.intensity = 0.0;
    temp.accuracy = 0.0;
  }
  temp.peakweight = peakweight + p.peakweight;
  temp.pairweight = pairweight + p.pairweight;

  return temp;
}

namespace std{
  ostream& operator<<( ostream& os, peakres& pr ){
    return (os << "acc: " << pr.accuracy << ", peakw: " << pr.peakweight);
  }
};
