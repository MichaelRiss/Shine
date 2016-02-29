cmake_minimum_required( VERSION 3.2.3 )
project( plain2strip )

set( CMAKE_C_FLAGS "-D_GNU_SOURCE" )
add_executable( plain2strip plain2strip.c )

install(PROGRAMS plain2strip DESTINATION bin)
