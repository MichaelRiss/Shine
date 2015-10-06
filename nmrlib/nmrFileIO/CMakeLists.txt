cmake_minimum_required( VERSION 3.2.3 )
project( nmrFileIO )

set( SOURCES ReadSequence.cc ReadShiftlist.cc ReadPeaklist.cc aminoshift.cc 
	     shiftpair.cc ReadNMRViewPar.cc ReadBetaSheet.cc ReadNogo.cc 
	     FileNotFoundException.cc )

set( CMAKE_CXX_FLAGS "-Wall -std=c++11 -fPIC")

add_library( nmrFileIO SHARED ${SOURCES} )
add_library( nmrFileIO_static STATIC ${SOURCES} )
