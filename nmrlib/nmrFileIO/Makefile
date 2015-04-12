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

OBJECTS = ReadSequence.o ReadShiftlist.o ReadPeaklist.o aminoshift.o \
	  shiftpair.o ReadNMRViewPar.o ReadBetaSheet.o ReadNogo.o \
	  FileNotFoundException.o

all:	libnmrFileIO.a libnmrFileIO.so

libnmrFileIO.a:	$(OBJECTS)
	ar rs libnmrFileIO.a $(OBJECTS)
	cp $@ ../lib

libnmrFileIO.so: $(OBJECTS)
	c++ -shared -Wl,-soname,libnmrFileIO.so -o $@ $(OBJECTS)
	cp $@ ../lib

.cc.o:
	c++ -Wall -std=c++11 -fPIC -g -c $<

clean:
	rm -f $(OBJECTS) libnmrFileIO.a libnmrFileIO.so

distclean:
	rm -f *~
