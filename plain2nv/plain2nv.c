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

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <errno.h>

extern int errno;

int main( int argc, char** argv ){
  int a, b;
  int numdim;
  int i;
  int* plainres;
  int* tileres;
  int* matrixres;
  int* plaincoord;
  int* matrixcoord;
  int* tilecoord;
  int plainstride;
  int matrixstride;
  int tilestride;
  int plainsize;
  FILE* par_in;
  int data_in;
  int out;
  float* source;
  float* dest;
  void* buf;
  int result;
  int quit;
  int dindex;
  int sindex;
  int tilesize;

  if( argc != 3 ){
    printf( "Usage: plain2nv <input>.plain|.par <output.mat>\n" );
    exit(0);
  }

  char datafile[1024];
  char parfile[1024];
  
  snprintf( parfile, 1024, "%s.par", argv[1] );
  snprintf( datafile, 1024, "%s.plain", argv[1] );
  
  par_in = fopen( parfile, "r" );    
  data_in = open( datafile, O_RDONLY|O_LARGEFILE );
  out = open( argv[2], O_CREAT|O_TRUNC|O_RDWR, 
	      S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH );
  
  fscanf( par_in, "header %i %i\n", &a, &b );
  fscanf( par_in, "dim %i", &numdim );
  plainres = (int*) calloc( numdim, sizeof(int) );
  tileres = (int*) calloc( numdim, sizeof(int) );
  matrixres = (int*) calloc( numdim, sizeof(int) );
  plaincoord = (int*) calloc( numdim, sizeof(int) );
  matrixcoord = (int*) calloc( numdim, sizeof(int) );
  tilecoord = (int*) calloc( numdim, sizeof(int) );
  plainsize = sizeof(float);
  for( i = 0; i < numdim; i++ ){
    fscanf( par_in, " %i", &plainres[i] );
    plainsize *= plainres[i];
  }
  for( i = 0; i < numdim; i++ ){
    fscanf( par_in, " %i", &tileres[i] );
  }
  
  buf = calloc( 1048576, 1 );
  
  for( i = 0; i < plainsize/1048576; i++ ){
    write( out, buf, 1048576 );
  }
  write( out, buf, plainsize % 1048576 );

  source = (float*) mmap( NULL, plainsize, PROT_READ, MAP_PRIVATE, 
			  data_in, 0 );
  if( source == (void*) -1 ){
    printf( "Error: 1. mmap did not succeed.\n" );
    exit( 1 );
  }
  dest = (float*) mmap( NULL, plainsize, PROT_READ|PROT_WRITE, MAP_SHARED, 
			out, 0 );
  if( dest == (void*) -1 ){
    printf( "Error: 2. mmap did not succeed.\n" );
    exit( 1 );
  }
  
  tilesize = 1;
  for( i = 0; i < numdim; i++ ){
    matrixres[i] = plainres[i]/tileres[i];
    tilesize *= tileres[i];
  }

  quit = 0;
  while( !quit ){
    /* translate plain coords in matrix/tile coords */
    for( i = 0; i < numdim; i++ ){
      matrixcoord[i] = plaincoord[i] / tileres[i];
      tilecoord[i] = plaincoord[i] % tileres[i];
    }
    
    /* compute memory addresses */
    sindex = 0;
    plainstride = 1;
    for( i = 0; i < numdim; i++ ){
      sindex += plaincoord[i] * plainstride;
      plainstride *= plainres[i];
    }
    
    dindex = 0;
    tilestride = 1;
    matrixstride = 1;
    for( i = 0; i < numdim; i++ ){
      dindex += matrixcoord[i] * matrixstride * tilesize;
      matrixstride *= matrixres[i];
      dindex += tilecoord[i] * tilestride;
      tilestride *= tileres[i];
    }
    
    /* copy */
    dest[dindex] = source[sindex];
    /*       printf( "sindex: %i, dindex: %i\n", sindex, dindex ); */
    
    /* increase coords */
    plaincoord[0]++;
    for( i = 0; i < numdim; i++ ){
      if( plaincoord[i] >= plainres[i] ){
	plaincoord[i] = 0;
	if( i >= numdim-1 ){
	  quit = 1;
	}
	else{
	  plaincoord[i+1]++;
	}
      }
    }
  }
  
  result = msync( dest, plainsize, MS_SYNC );
  if(result != 0)
    perror( "Da Fehla" );
  munmap( dest, plainsize );
  munmap( source, plainsize );
  
  close( out );
  close( data_in );
  fclose( par_in );
  return 0;
}
