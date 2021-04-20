
/*
	Routines to handle unformatted access to fortran files

	Here the c-interface is given
	For the Fortran interface please use unf_fortran.c
*/

/*
	new:

	get_next_record_32
	get_next_record_64
	read_record_32
	read_record_64
	write_record_32
	write_record_64
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "unformatted.h"

/*****************************************************/
/*       parameters to be changed accordingly        */
/*****************************************************/

/*
  iswap gives number of bytes to swap per variable
  iswap = 0 does not swap
  iswap = 4 swaps regular integer/real
  iswap = 8 swaps double precision

  ihead is length of header/trailer of record
  its contents indicates the length of the record
  on 32 bit architecture it should be 4
  on 64 bit architecture it should be 8

  RSIZE must correspond to the length of ihead
  use int for 4 bytes
  use long long for 8 bytes
*/

#define LONGLONG long long

static int iswap = 0;			/* how many bytes to swap */
static int ihead = 4;			/* length of header of record */
static int debug = 0;

/* #define RSIZE LONGLONG */
#define RSIZE int

/*****************************************************/
/*     do not change anything beyond this point      */
/*****************************************************/

static int ibuffer = 0;			/* size of buffer */
static int iactrec = 0;			/* size of actual record */
static unsigned char *buffer;		/* global buffer for read */

static FILE *filein;
static FILE *fileout;

/*****************************************************/
/*     checks if RSIZE and ihead are compatible      */
/*****************************************************/

void check_rsize( void )
{
  int rsize = sizeof( RSIZE );

  if( rsize != ihead ) {
    Error2(1,"error in size of header: rsize = %d   ihead = %d\n",rsize,ihead);
  }
}

/*****************************************************/

void copy_record_to_array( int *ndim, int *ibuf, void *buf )
{
  /* ndim and ibuf is number of bytes !!! */

  int iret;

  iret = read_record();

  *ibuf = iactrec;

  if( iactrec > *ndim ) {
    Error2(1,"dimension error in: ndim = %d   nsize = %d\n",*ndim,iactrec);
  }

  buf = memcpy(buf, (void *) buffer, (size_t) iactrec);
}

/*****************************************************/

int get_next_record( int *ibuf, void **buf )
{
  int iret;

  iret = read_record();

  *ibuf = iactrec;
  *buf = buffer;
  return iret;
}

int get_next_record_32( int *ibuf, void **buf )
{
  int iret;

  iret = read_record_32();

  *ibuf = iactrec;
  *buf = buffer;
  return iret;
}

int get_next_record_64( int *ibuf, void **buf )
{
  int iret;

  iret = read_record_64();

  *ibuf = iactrec;
  *buf = buffer;
  return iret;
}

int skip_records( int nrec )
{
  while( nrec && read_record() ) {
    nrec--;
  }
  return nrec;		/* is 0 for success */
}

int copy_records( int nrec )
{
  while( nrec && read_record() ) {
    nrec--;
    write_record();
  }
  return nrec;		/* is 0 for success */
}

/*****************************************************/

void set_read_file( char *string )
{
  if( filein && filein != stdin ) {
    if( fclose(filein) ) Error(1,"error closing read file\n");
  }
  filein = fopen(string,"r");
  if( !filein ) ErrorC(1,"error opening file %s\n",string);

  if( debug ) fprintf(stderr,"file %s opened for reading\n",string);
}

void set_write_file( char *string )
{
  if( fileout && fileout != stdout ) {
    if( fclose(fileout) ) Error(1,"error closing write file\n");
  }
  fileout = fopen(string,"w");
  if( !fileout ) ErrorC(1,"error opening file %s\n",string);

  if( debug ) fprintf(stderr,"file %s opened for writing\n",string);
}

void close_read_file( void )
{
  if( filein && filein != stdin ) {
    if( fclose(filein) ) Error(1,"error closing read file\n");
  }
}

void close_write_file( void )
{
  if( fileout && fileout != stdout ) {
    if( fclose(fileout) ) Error(1,"error closing read file\n");
  }
}

void init_read_file( void )
{
  if( filein && filein != stdin ) {
    fseek(filein,0,SEEK_SET);		/* positions at start of file */
  }
}

/*****************************************************/

int *get_int( void *buf )
{
  return (int *) buf;
}

char *get_char( void *buf )
{
  return (char *) buf;
}

float *get_float( void *buf )
{
  return (float *) buf;
}

/*****************************************************/

int read_record( void )
{
  RSIZE nrec1,nrec2;
  int n;
  int n1,n2;

  check_rsize();
  alloc_file_desc();

  // fprintf(stderr,"ihead: %d\n",ihead);

  n = fread( &nrec1, ihead, 1, filein );
  if( n == 0 ) {
    return 0;
  }
  if( iswap ) swap_byte_order(ihead,(void *) &nrec1);

  // fprintf(stderr,"first: %d %d\n",n,nrec1);

  n = nrec1;			/* cast if necessary */
  read_data(n);			/* reads n bytes */

  n = fread( &nrec2, ihead, 1, filein );
  if( iswap ) swap_byte_order(ihead,(void *) &nrec2);

  // fprintf(stderr,"second: %d %d\n",n,nrec1);

  if( n != 1 ) {
    Error1(1,"read error in read_record: %d\n",n);
  } else if( nrec1 != nrec2 ) {
    n1 = nrec1; n2 = nrec2;
    Error2(1,"internal error in unformatted record: %d %d\n",n1,n2);
  }

  iactrec = nrec1;
  return 1;
}

int read_record_32( void )
{
  int nrec1,nrec2;
  int n;
  int n1,n2;
  int ihead32;

  check_rsize();
  alloc_file_desc();

  ihead32 = 4;

  n = fread( &nrec1, ihead32, 1, filein );
  if( n == 0 ) {
    return 0;
  }
  if( iswap ) swap_byte_order(ihead32,(void *) &nrec1);

  /* fprintf(stderr,"first: %d %d\n",n,nrec1); */

  n = nrec1;			/* cast if necessary */
  read_data(n);			/* reads n bytes */

  n = fread( &nrec2, ihead32, 1, filein );
  if( iswap ) swap_byte_order(ihead,(void *) &nrec2);

  if( n != 1 ) {
    Error1(1,"read error in read_record: %d\n",n);
  } else if( nrec1 != nrec2 ) {
    n1 = nrec1; n2 = nrec2;
    Error2(1,"internal error in unformatted record: %d %d\n",n1,n2);
  }

  iactrec = nrec1;
  return 1;
}

int read_record_64( void )
{
  LONGLONG nrec1,nrec2;
  int n;
  int n1,n2;
  int ihead64;

  check_rsize();
  alloc_file_desc();

  ihead64 = 8;

  n = fread( &nrec1, ihead64, 1, filein );
  if( n == 0 ) {
    return 0;
  }
  if( iswap ) swap_byte_order(ihead64,(void *) &nrec1);

  /* fprintf(stderr,"first: %d %d\n",n,nrec1); */

  n = nrec1;			/* cast if necessary */
  read_data(n);			/* reads n bytes */

  n = fread( &nrec2, ihead64, 1, filein );
  if( iswap ) swap_byte_order(ihead,(void *) &nrec2);

  /* fprintf(stderr,"second: %d %d\n",n,nrec2); */

  if( n != 1 ) {
    Error1(1,"read error in read_record: %d\n",n);
  } else if( nrec1 != nrec2 ) {
    n1 = nrec1; n2 = nrec2;
    Error2(1,"internal error in unformatted record: %d %d\n",n1,n2);
  }

  iactrec = nrec1;
  return 1;
}

int write_record( void )
{
  RSIZE nrec;
  int n;

  check_rsize();
  alloc_file_desc();
  nrec = iactrec;

  n = fwrite( &nrec, ihead, 1, fileout );
  if( n != 1 ) Error1(1,"write error in write_record: %d\n",n);

  write_data(iactrec);

  n = fwrite( &nrec, ihead, 1, fileout );
  if( n != 1 ) Error1(1,"write error in write_record: %d\n",n);

  return 1;
}

int write_record_32( void )
{
  int n;
  int nrecc32,ihead32;

  check_rsize();
  alloc_file_desc();

  nrecc32 = iactrec;
  ihead32 = 4;

  n = fwrite( &nrecc32, ihead32, 1, fileout );
  if( n != 1 ) Error1(1,"write error in write_record: %d\n",n);

  write_data(iactrec);

  n = fwrite( &nrecc32, ihead32, 1, fileout );
  if( n != 1 ) Error1(1,"write error in write_record: %d\n",n);

  return 1;
}

int write_record_64( void )
{
  int n;
  LONGLONG nrecc64,ihead64;

  check_rsize();
  alloc_file_desc();

  nrecc64 = iactrec;
  ihead64 = 8;

  n = fwrite( &nrecc64, ihead64, 1, fileout );
  if( n != 1 ) Error1(1,"write error in write_record: %d\n",n);

  write_data(iactrec);

  n = fwrite( &nrecc64, ihead64, 1, fileout );
  if( n != 1 ) Error1(1,"write error in write_record: %d\n",n);

  return 1;
}

/*****************************************************/

void set_swap_byte_order( int nbyte )
{
  iswap = nbyte;
}

void swap_byte_order( int nbyte , void *buffer )
{
  int i,j;
  char unsigned aux;
  char unsigned *cbuffer;

  cbuffer = (char unsigned *) buffer;

  for(i=0,j=nbyte-1;i<j;i++,j--) {
    aux = cbuffer[i];
    cbuffer[i] = cbuffer[j];
    cbuffer[j] = aux;
  }
}

void swap_array( int n , int nbyte , void *buffer )
{
  int i;
  void *p;

  if( n%nbyte != 0 ) {
    Error2(5,"Not a multiple of nbyte: ",n,nbyte);
  }

  p = buffer;
  for(i=0;i<n;i+=nbyte) {
    swap_byte_order( nbyte , p );
    p += nbyte;
  }
}

/*****************************************************/

void skip_bytes( int nbytes )
{
  int ib;
  int n;

  while( nbytes-- ) {
    n = fread( &ib, 1, 1, filein );
    if( n != 1 ) Error1(1,"read error in skip_bytes: %d\n",n);
  }
}

void read_data( int nbytes )
{
  int n;

  alloc_buffer(nbytes);
  n = fread( buffer, 1, nbytes, filein );
  if( n != nbytes ) Error2(1,"read error in read_data: %d %d\n",n,nbytes);
  if( iswap ) swap_array(n,iswap,(void *) buffer);
}

void write_data( int nbytes )
{
  int n;

  n = fwrite( buffer, 1, nbytes, fileout );
  if( n != nbytes ) Error2(1,"write error in write_data: %d %d\n",n,nbytes);
}

/***********************************************************/

void alloc_buffer( int nbytes )
{
  if( nbytes <= ibuffer ) return;

  if( buffer ) {
    free(buffer);
  }
  buffer = (unsigned char *) malloc(nbytes);
  if( !buffer ) {
    fprintf(stderr,"read allocating buffer: %d\n",nbytes);
    exit(1);
  }
  //fprintf(stderr,"Successfully allocated buffer: %d\n",nbytes);
  ibuffer = nbytes;
}

void alloc_file_desc( void )
{
  if( !filein ) filein  = stdin;
  if( !fileout ) fileout = stdout;
}
  
/***********************************************************/

void Error( int ier, char *string )
{
  int dummy;
  fprintf(stderr,string,dummy);
  exit(ier);
}

void Error1( int ier, char *string, int i1 )
{
  fprintf(stderr,string,i1);
  exit(ier);
}

void Error2( int ier, char *string, int i1, int i2 )
{
  fprintf(stderr,string,i1,i2);
  exit(ier);
}

void ErrorC( int ier, char *string, char *s )
{
  fprintf(stderr,string,s);
  exit(ier);
}

/***********************************************************/

