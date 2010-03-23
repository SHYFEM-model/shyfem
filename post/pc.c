
#include "pspost.h"

int main( void )

{
	PsGraphOpen();
	PsStartPage();
	PsLine(1.,1.,10.,10.);
	PsEndPage();
	PsStartPage();
	PsText(2.,2.,"Hello World!");
	PsEndPage();
	PsGraphClose();
	return 0;
}
