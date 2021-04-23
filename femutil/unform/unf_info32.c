
#include <stdio.h>
#include <stdlib.h>
#include "unformatted.h"

/*--------------- please change parameters ------------------*/

int iswap = 4;		/* swap byte order ? - use 4 for 32 bit */
int iwrite = 1;		/* write transformed records ? */

/* iwrite = 1 is only usefull with iswap != 0 */

/*------------------ end of parameters ----------------------*/

int main( int argc, char *argv[] )
{
  int irec = 0;
  int ibuf;
  void *buffer;

  if( iswap ) set_swap_byte_order(iswap);

  while( get_next_record_32( &ibuf, &buffer ) ) {
    irec++;
    fprintf(stderr,"record = %d    record length = %d\n",irec,ibuf);

    if( iwrite && ! write_record_32() ) {
      fprintf(stderr,"error in writing record");
      exit(11);
    }
  }

  return 0;
}

