/*  Shine - NMR NOESY simulation program
    Copyright (C) 2016  Murray Coles

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
  int numdim=3;
  int i;
  float* freq;
  float* spwd;
  int* plainres;
  int* matrixres;
  int* plaincoord;
  int* matrixcoord;
  int* tilecoord;
  int plainsize;
  FILE* par_in;
  int data_in;
  FILE* outstrip;
  FILE* outplane;
  int out;
  float* source;
  float* dest;
  void* buf;
  int result;
 
 
  if( argc != 6 ){
    printf( "Usage: plain2strip <input.plain> <simulationparameters> <output.mat> <output.plane> <output.strip>\n" );
    exit(0);
  }

  char datafile[1024];
  char parfile[1024];
  
  snprintf( parfile, 1024, "%s", argv[2] );
  snprintf( datafile, 1024, "%s", argv[1] );
  
  par_in = fopen( parfile, "r" );    
  data_in = open( datafile, O_RDONLY|O_LARGEFILE );
  out = open( argv[3], O_CREAT|O_TRUNC|O_RDWR, 
	      S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH );
  outstrip = fopen(argv[5], "w" );
  outplane = fopen(argv[4], "wb" );


  freq = (float*) calloc( numdim, sizeof(float) );
  spwd = (float*) calloc( numdim, sizeof(float) );

  plainres = (int*) calloc( numdim, sizeof(int) );
  matrixres = (int*) calloc( numdim, sizeof(int) );
  plaincoord = (int*) calloc( numdim, sizeof(int) );
  matrixcoord = (int*) calloc( numdim, sizeof(int) );
  tilecoord = (int*) calloc( numdim, sizeof(int) );
  plainsize = sizeof(float);


  fscanf( par_in, "#freq: %f %f %f\n", &freq[0], &freq[1], &freq[2] );
  fscanf( par_in, "#spwd: %f %f %f\n", &spwd[0], &spwd[1], &spwd[2] );
  fscanf( par_in, "#size: %i %i %i\n", &plainres[0], &plainres[1], &plainres[2] );
  
  
  for( i = 0; i < numdim; i++ ){
    plainsize *= plainres[i];
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


  int j=0;
  int psize=plainres[0]*plainres[1];
  float *plane;
  plane = calloc(psize, sizeof(float));
  int planeidx = 0;
  
  for( i = plainres[0]*plainres[1]*2; i < plainres[0]*plainres[1]*3; i++ ){
    plane[planeidx] = source[i];
          
    if ((i-(plainres[0]/2)) % plainres[0] == 0){
      fprintf(outstrip, "%i %f \n", j, source[i] );
      j++;
    }
    planeidx ++;

  }
  
  unsigned int * fpnt = (unsigned int *) plane;
  for( i = 0; i < plainres[0]*plainres[1]; i++ ){

  /* endian swap */
   
    fpnt[i] = (fpnt[i]<<24) | 
                  ((fpnt[i]<<8) & 0x00ff0000) |
                  ((fpnt[i]>>8) & 0x0000ff00) | 
                  (fpnt[i]>>24);
  }

  int rtn = fwrite(plane, sizeof(float), plainres[0]*plainres[1], outplane);
  printf("size=%d\n",rtn);
  
  fclose(outstrip);
  fclose( outplane );

  free(plane);
  plane=NULL;  
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
