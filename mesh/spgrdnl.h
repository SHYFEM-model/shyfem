
#include "general.h"
#include "grd.h"

/****************************************************************/

void MakeNodelList( void );
void PrintNodel( int first );
int GetBoxOfNodel( int first , int second );
void MakeLinelList( void );
int CheckNeibBoxes(int number,int *boxes,int *boxesc,int nm);
void WriteNeibBoxes(int box, int *boxes,int *boxesc,int nm, int nbox);
void PrintLine( Line_type *pl );

int *GetNextBox( int nm , int *boxes , int * ind , int *n , int *nbox );
void CheckLinel( void );
void MakeLines( void );
void MakeNodes( void );
