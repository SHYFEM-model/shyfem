
/*
	Fortran interface for unformatted routines

	not all calls have been implemented
*/

#include "unformatted.h"

/*****************************************************/

void get_array__( int *ndim, int *ibuf, void *buf )
{
  /* ndim and ibuf are number of bytes !!! */

  copy_record_to_array( ndim , ibuf , buf );
}

void get_real_array__( int *ndim, int *ibuf, void *buf )
{
  get_array__( ndim, ibuf, buf );
}

void get_int_array__( int *ndim, int *ibuf, void *buf )
{
  get_array__( ndim, ibuf, buf );
}

void set_fortran_read_file__( char *string )
{
  set_read_file( string );
}

void close_fortran_read_file__( void )
{
  close_read_file();
}

/*****************************************************/

