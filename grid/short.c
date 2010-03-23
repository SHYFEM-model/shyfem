
#include <stdio.h>
#include "graph.h"

void main( void )

{
        QGraphInit();

        printf("New version 2.10\n");
        press_any_key();

	QText(0.3,0.5,"Hello, world!");
	QLine(0.2,0.3,0.8,0.7);

        press_any_key();
        press_any_key();

        QGraphClose();
}

