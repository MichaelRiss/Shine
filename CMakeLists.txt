cmake_minimum_required( VERSION 3.2.3 )
project( Shine )

find_package( GSL )

add_subdirectory(nmrlib)
add_subdirectory(plain2nv)

set( SOURCES NMRparseparameters.cc BuildSequence.cc NMRdata.cc
	     Shiftfileparser.cc Shine.cc )

set( M3D_SOURCES make3D_C++.cc )

set( CMAKE_CXX_FLAGS "-Wall -std=c++11 -I${CMAKE_SOURCE_DIR}/nmrlib/include" )

add_executable( Shine ${SOURCES} )
add_executable( make3D_C++ ${M3D_SOURCES} )
target_link_libraries( Shine ${GSL_LIBRARIES} pdbimport_static
		       nmrFileIO_static nmrutil_static )