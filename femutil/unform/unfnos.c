
/*
	Manipulation of nos files through c routines
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "unformatted.h"

int get_records( void );

void copy_nos_header( void );
void skip_nos_header( void );
int copy_nos_data( void );
void skip_nos_record( int irec );
void copy_nos_record( int irec );
int count_nos_records( void );

int read_nos( void );
int read_nos_header( void );
int read_nos_data( void );

int Adjust_irec( int irec );

void Usage( void );

int main( int argc, char *argv[] )
{
  int irec = 0;
  int i,ir;
  char *what;
  char *file;
  char *outfile;

  if( argc < 2 ) {
    fprintf(stderr,"arguments %d\n",argc);
    Usage();
  } else {
    what = argv[1];
  }

  if( strcmp(what,"nosinfo") == 0 ) {
    for(i=2;i<argc;i++) {
      file = argv[i];
      set_read_file(file);
      irec = read_nos();
      fprintf(stderr,"A total of %d records read for %s\n",irec,file);
    }
  } else if( strcmp(what,"noscount") == 0 ) {
    for(i=2;i<argc;i++) {
      file = argv[i];
      set_read_file(file);
      irec = count_nos_records();
      fprintf(stderr,"File %s contains %d data records\n",file,irec);
    }
  } else if( strcmp(what,"nosextract") == 0 ) {
    if( argc < 4 ) Usage();
    irec = atoi(argv[2]);
    outfile = argv[3];
    set_write_file(outfile);
    for(i=4;i<argc;i++) {
      file = argv[i];
      set_read_file(file);
      ir = Adjust_irec(irec);	/* must do before handling header */
      if( i == 4 ) {
	copy_nos_header();
      } else {
	skip_nos_header();
      }
      if( ir ) {
        copy_nos_record(ir);
        fprintf(stderr,"record %d has been extracted from %s\n",ir,file);
      }
    }
  } else if( strcmp(what,"noscopy") == 0 ) {
    if( argc < 3 ) Usage();
    outfile = argv[2];
    set_write_file(outfile);
    for(i=3;i<argc;i++) {
      file = argv[i];
      set_read_file(file);
      if( i == 3 ) {
	copy_nos_header();
      } else {
	skip_nos_header();
      }
      irec = copy_nos_data();
      fprintf(stderr,"A total of %d records copied from %s\n",irec,file);
    }
  } else {
    fprintf(stderr,"what: %s\n",what);
    Usage();
  }

  return 0;
}

void Usage( void )
{
  fprintf(stderr,"Usage: unf what options\n");
  fprintf(stderr,"   unf nosinfo [files]\n");
  fprintf(stderr,"   unf noscount [files]\n");
  fprintf(stderr,"   unf noscopy outfile [files]\n");
  fprintf(stderr,"   unf nosextract record outfile [files]\n");
  exit(0);
}

int copy_file( void )
{
  int irec = 0;

  while( copy_records(1) ) {
    irec++;
  }
  return irec;
}

int get_records( void )
{
  int irec = 0;
  int n;
  void *buffer;

  while( get_next_record(&n,&buffer) ) {
    irec++;
    fprintf(stderr,"Record %d  consists of %d Bytes\n",irec,n);
  }

  return irec;
}

/********************************************************/

void copy_nos_header( void )
{
  if( copy_records(6) ) {
    Error(3,"Error copying nos header\n");
  }
}

void skip_nos_header( void )
{
  if( skip_records(6) ) {
    Error(4,"Error skipping nos header\n");
  }
}

int copy_nos_data( void )
{
  int irec = 0;

  while( !copy_records(2) ) {
    irec++;
  }

  return irec;
}

void skip_nos_record( int irec )
{
  while( irec && !skip_records(2) ) {
    irec--;
  }
}

void copy_nos_record( int irec )
{
  skip_nos_record(irec-1);
  copy_records(2);
}

int count_nos_records( void )
{
  int nrec = 0;

  skip_nos_header();
  while( !skip_records(2) ) {
    nrec++;
  }

  return nrec;
}

int Adjust_irec( int irec )
{
  int nrec,ir;

  nrec = count_nos_records();
  init_read_file();

  if( nrec == 0 || irec == 0 ) {
    ir = 0;
  } else if( irec > nrec ) {
    ir = 0;
  } else if( irec < 0 ) {
    ir = -irec;
    if( ir > nrec ) {
      ir = 0;
    } else {
      ir = nrec + 1 - ir;
    }
  } else {
    ir = irec;
  }

  fprintf(stderr,"Adjust_irec: %d %d %d\n",nrec,irec,ir);
  
  return ir;
}

/********************************************************/

int read_nos( void )
{
  int iret;

  iret = read_nos_header();
  if( iret < 0 ) return iret;
  return read_nos_data();
}

int read_nos_header( void )
{
  int irec = 0;
  int n;
  void *buffer;
  int *ibuf;
  char *cbuf;

  if( !get_next_record(&n,&buffer) ) return(-1);
  ibuf = get_int(buffer);
  fprintf(stderr,"Header: %d %d\n",ibuf[0],ibuf[1]);

  if( !get_next_record(&n,&buffer) ) return(-2);
  ibuf = get_int(buffer);
  fprintf(stderr,"Params: %d %d %d %d\n",ibuf[0],ibuf[1],ibuf[2],ibuf[3]);

  if( !get_next_record(&n,&buffer) ) return(-3);
  cbuf = get_char(buffer);
  fprintf(stderr,"Text  : %s\n",cbuf);

  skip_records(3);

  return 0;
}

int read_nos_data( void )
{
  int irec = 0;
  int n;
  int i1,i2;
  void *buffer;
  int *ibuf;
  char *cbuf;

  while( get_next_record(&n,&buffer) ) {
    ibuf = get_int(buffer);
    i1 = ibuf[0];
    i2 = ibuf[1];
    get_next_record(&n,&buffer);
    irec++;
    fprintf(stderr,"Record %d: %d %d %d\n",irec,i1,i2,n);
  }

  return irec;
}

