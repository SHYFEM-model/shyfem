
#include <stdio.h>
#include <stdlib.h>
#include "unformatted.h"

/*--------------- please change parameters ------------------*/

int iswap = 0;		/* swap byte order ? */
int iwrite = 0;		/* write transformed records ? */

/* iwrite = 1 is only usefull with iswap != 0 */

/*------------------ end of parameters ----------------------*/

/*
unsigned char buffer[10];
FILE *ptr;
ptr = fopen("test.bin","rb");  // r for read, b for binary
fread(buffer,sizeof(buffer),1,ptr); // read 10 bytes to our buffer
*/



int main( int argc, char *argv[] )
{
  int irec = 0;
  int ibuf;
  void *buffer;
  FILE *ptr1, *ptr2;

  if(argc < 3){
    printf("need two files...\n");
    exit(1);
  }
  printf("First file is: %s\n",argv[1]);
  printf("Second file is: %s\n",argv[2]);

  ptr1 = fopen(argv[1],"rb");
  ptr2 = fopen(argv[2],"rb");

  exit(0);	/* program is not finished */

  if( iswap ) set_swap_byte_order(iswap);

  while( get_next_record( &ibuf, &buffer ) ) {
    irec++;
    fprintf(stderr,"record = %d    record length = %d\n",irec,ibuf);

    if( iwrite && ! write_record() ) {
      fprintf(stderr,"error in writing record");
      exit(11);
    }
  }

  return 0;
}

