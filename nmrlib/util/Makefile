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

LIBOBJECTS = amino.o peakpickutil.o iupac2xplor_pdb.o

all:	libnmrutil.a iupac2xplor_pdb libnmrutil.so

iupac2xplor_pdb:	iupac2xplor_pdb.o libnmrutil.a
	c++ -Wall -g -o $@ $< libnmrutil.a

libnmrutil.a:	$(LIBOBJECTS)
	ar r libnmrutil.a $(LIBOBJECTS)
	cp $@ ../lib

libnmrutil.so: $(LIBOBJECTS)
	c++ -shared -Wl,-soname,libnmrutil.so -o $@ $(LIBOBJECTS)
	cp $@ ../lib

.cc.o:
	c++ -Wall -fPIC -std=c++11 -g -c $<

clean:
	rm -f $(LIBOBJECTS) libnmrutil.a libnmrutil.so \
	iupac2xplor_pdb libutil.a

distclean:
	rm -f *~
