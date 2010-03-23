
#include <stdio.h>

int main( void )

{
	int i;
	float x,dx,y,dy;

	x = 10.;
	y = 15.;

	printf("%%!\n");

	printf("28.3 dup scale\n");
	printf("1 setlinecap 1 setlinejoin 0.01 setlinewidth\n");

	printf("-5 %f moveto\n",y);
	printf("35 %f lineto\n",y);
	printf("%f -5 moveto\n",x);
	printf("%f 35 lineto\n",x);
	printf("stroke\n");

	for(i=-5;i<=35;i++) {
	  dy = (i%5 == 0) ? .3 : .1;
	  printf("%d %f moveto\n",i,y);
	  printf("%d %f lineto\n",i,y+dy);
	}
	printf("stroke\n");

	for(i=-5;i<=35;i++) {
	  dx = (i%5 == 0) ? .3 : .1;
	  printf("%f %d moveto\n",x,i);
	  printf("%f %d lineto\n",x+dx,i);
	}
	printf("stroke\n");

	printf("showpage\n");

	return 0;
}
