cmake_minimum_required( VERSION 3.2.3 )
project( nmrutil )

set( SOURCES amino.cc peakpickutil.cc iupac2xplor_pdb.cc )

set( CMAKE_CXX_FLAGS "-Wall -std=c++11 -fPIC")

add_library( nmrutil SHARED ${SOURCES} )
add_library( nmrutil_static STATIC ${SOURCES} )
add_executable( iupac2xplor_pdb ${SOURCES} )
