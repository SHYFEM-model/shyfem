
#include <stdio.h>
#include <stdlib.h>
#include "unformatted.h"

int iswap = 4;
int iwrite = 1;

int mn = 154*248;

int main( int argc, char *argv[] )
{
  int irec = 0;
  int ibuf;
  void *buffer;

  if( iswap ) set_swap_byte_order(iswap);

  while( get_next_record( &ibuf, &buffer ) ) {
    irec++;
    fprintf(stderr,"record = %d    record length = %d\n",irec,ibuf);

    if( irec == 4 ) {
      fprintf(stderr,"...record %d adjusted\n",irec);
      swap_array(mn*2,4,buffer);
      swap_array(mn*2,2,buffer);
    }
    if( irec == 5 ) {
      fprintf(stderr,"...record %d adjusted\n",irec);
      swap_array(mn*4,4,buffer);
      swap_array(mn*4,2,buffer);
    }

    if( iwrite && ! write_record() ) {
      fprintf(stderr,"error in writing record");
      exit(11);
    }
  }

  return 0;
}

