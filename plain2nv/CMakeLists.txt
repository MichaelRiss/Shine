cmake_minimum_required( VERSION 3.2.3 )
project( plain2nv )

set( CMAKE_C_FLAGS "-D_GNU_SOURCE" )
add_executable( plain2nv plain2nv.c )

install(PROGRAMS plain2nv DESTINATION bin)
