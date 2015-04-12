# Shine - NMR NOESY simulation program
# Copyright (C) 2015  Michael Riss
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#    
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#    
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

INCLUDE = -Inmrlib/include
LIBS = -lgsl -lgslcblas nmrlib/lib/libpdbimport.a nmrlib/lib/libnmrFileIO.a \
	nmrlib/lib/libnmrutil.a
OBJECTS = NMRparseparameters.o BuildSequence.o NMRdata.o Shiftfileparser.o \
	Shine.o
MOBJECTS = make3D_C++.o

all:	nmrlib Shine make3D_C++ plain2nv

clean:
	$(MAKE) -C nmrlib clean
	$(MAKE) -C plain2nv clean
	rm -f *.o Shine make3D_C++

distclean: clean
	$(MAKE) -C nmrlib distclean
	$(MAKE) -C plain2nv distclean
	rm -f *~

nmrlib:
	$(MAKE) -C nmrlib

plain2nv:
	$(MAKE) -C plain2nv

.PHONY:	nmrlib plain2nv


Shine:	$(OBJECTS)
	c++ -o Shine $(OBJECTS) $(LIBS)

make3D_C++:	$(MOBJECTS)
	c++ -o make3D_C++ $(MOBJECTS) $(LIBS)


%.o:	%.cc
	c++ -Wall -O2 -std=c++11 $(INCLUDE) -c -o $@ $<

