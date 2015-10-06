cmake_minimum_required( VERSION 3.2.3 )
project( libpdbimport )

set( SOURCES pdbimport.cc )

set( CMAKE_CXX_FLAGS "-Wall -fPIC")

add_library( pdbimport SHARED ${SOURCES} )
add_library( pdbimport_static STATIC ${SOURCES} )
