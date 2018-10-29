
/************************************************************************\
 *
 *    Copyright (C) 1985-2018  Georg Umgiesser
 *
 *    This file is part of SHYFEM.
 *
 *    SHYFEM is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    SHYFEM is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with SHYFEM. Please see the file COPYING in the main directory.
 *    If not, see <http://www.gnu.org/licenses/>.
 *
 *    Contributions to this file can be found below in the revision log.
 *
\************************************************************************/


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
