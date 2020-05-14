
/************************************************************************\
 *
 *    Copyright (C) 2016,2019  Georg Umgiesser
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


/************************************************************************\
 *
 * routines to transform between color spaces
 *
 * revision log :
 *
 * 17.06.2016	ggu	changed VERS_7_5_15
 * 14.02.2019	ggu	changed VERS_7_5_56
 *
\************************************************************************/

#include <stdio.h>
#include "general.h"

#define	hmax	1.			/* h [0-1] */
/* #define	hmax	360. */		/* h [0-360] */

#define hdist	(hmax/6.)
#define hconv	(6./hmax)
#define undef	0.

#define abs(a)	( (a) > 0 ? (a) : (-(a)) )

void equaltb(float a, float b, float eps, int *berror);

/*********************************************************************/

float max3(float a, float b, float c)

{
	if( a > b ) {
	  if( a > c ) {
	    return a;
	  } else {
	    return c;
	  }
	} else {
	  if( c > b ) {
	    return c;
	  } else {
	    return b;
	  }
	}
}

float min3(float a, float b, float c)

{
	if( a < b ) {
	  if( a < c ) {
	    return a;
	  } else {
	    return c;
	  }
	} else {
	  if( c < b ) {
	    return c;
	  } else {
	    return b;
	  }
	}
}

/*********************************************************************/

void rgb2cmy(float r ,float g ,float b ,float *c ,float *m ,float *y)

/*
	rgb -> cmy

	real r,g,b	![0-1]
	real c,m,y	![0-1]
*/

{
	*c = 1. - r;
	*m = 1. - g;
	*y = 1. - b;
}


void cmy2rgb(float c ,float m ,float y ,float *r ,float *g ,float *b)

/*
	cmy -> rgb

	real c,m,y	![0-1]
	real r,g,b	![0-1]
*/

{
	*r = 1. - c;
	*g = 1. - m;
	*b = 1. - y;
}


void rgb2hsv(float r ,float g ,float b ,float *h ,float *s ,float *v)

/*
	rgb -> hsv		note: h is not [0-360] but [0-1]

	real r,g,b	![0-1]
	real h,s,v	![0-1]
*/

{

	float maxv,minv;
	float diff;
	float rdist,gdist,bdist;

	maxv = max3(r,g,b);
	minv = min3(r,g,b);
	diff = maxv - minv;

	*v = maxv;
	*s = 0.;
	if( maxv > 0. ) *s = diff/maxv;

	if( *s == 0. ) {
	  *h = undef;
	} else {
	  rdist = (maxv-r)/diff;
	  gdist = (maxv-g)/diff;
	  bdist = (maxv-b)/diff;
	  if( r == maxv ) {
	    *h = bdist - gdist;
	  } else if( g == maxv ) {
	    *h = 2. + rdist - bdist;
	  } else if( b == maxv ) {
	    *h = 4. + gdist - rdist;
	  } else {
	    Error("error stop rgb2hsv: internal error (1)");
	  }
	  *h = *h * hdist;
	  if( *h < 0. ) *h = *h + hmax;
	}

}

void hsv2rgb(float h, float s, float v, float *r, float *g, float *b)

/*
	hsv -> rgb		note: h is not [0-360] but [0-1]

	real h,s,v	![0-1]
	real r,g,b	![0-1]
*/

{
	int i;
	float p,q;

	i = hconv * h;
	i = i % 6;
	p = v * (h-i*hdist) / hdist;	/* rising */
	q = v - p;			/* falling */

	if( i == 0 ) {
	    *r=v;
	    *g=p;
	    *b=0;
	} else if( i == 1 ) {
	    *r=q;
	    *g=v;
	    *b=0;
	} else if( i == 2 ) {
	    *r=0;
	    *g=v;
	    *b=p;
	} else if( i == 3 ) {
	    *r=0;
	    *g=q;
	    *b=v;
	} else if( i == 4 ) {
	    *r=p;
	    *g=0;
	    *b=v;
	} else if( i == 5 ) {
	    *r=v;
	    *g=0;
	    *b=q;
	} else {
	    Error("error stop hsv2rgb: internal error (1)");
	}

	*r = v + (*r-v) * s;
	*g = v + (*g-v) * s;
	*b = v + (*b-v) * s;

}


#define	m	714025
#define	ia	1366
#define	ic	150889
#define	rm	1.4005112e-6

#define	nrdim	97
#define	nrdim1	(nrdim+1)

float randomtb(int *idump)

{
      static int ir[nrdim1];
      static int iy;
      static int iff=0;

      int j;
      int idum = *idump;
      float randomtb;

      if(idum==0 || iff==0) {
        iff=1;
        idum=(ic-idum)%m;
        for(j=1;j<=nrdim;j++) {
          idum=(ia*idum+ic)%m;
          ir[j]=idum;
	}
        idum=(ia*idum+ic)%m;
        iy=idum;
      }

      j=1+(97*iy)/m;
      if(j>97 || j<1) {
	Error("error stop ran2: internal error");
      }

      iy=ir[j];
      randomtb=iy*rm;
      idum=(ia*idum+ic)%m;
      ir[j]=idum;

      *idump = idum;
      return randomtb;
}


void test_ct()

{

/*
c tests conversion routine
*/

	int i;
	int itot,iseed;
	int berror;
	float r,g,b;
	float rn,gn,bn;
	float h,s,v;
	float hn,sn,vn;
	float eps;

	eps = 0.0001;
	itot = 100;
	iseed = 124357;
	berror = 0;

	for(i=1;i<=itot;i++) {
	      
	  h = randomtb(&iseed);
	  s = randomtb(&iseed);
	  v = randomtb(&iseed);

	  h = h * hmax;
	  hsv2rgb(h,s,v,&r,&g,&b);
	  rgb2hsv(r,g,b,&hn,&sn,&vn);

	  printf("%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n"
			,h,hn,s,sn,v,vn,r,g,b);

	  equaltb(h,hn,eps,&berror);
	  equaltb(s,sn,eps,&berror);
	  equaltb(v,vn,eps,&berror);
	  if( berror ) Error("error in computation h2r...");

	}
	      
	printf("-------------------------------------\n");

	for(i=1;i<=itot;i++) {

	  r = randomtb(&iseed);
	  g = randomtb(&iseed);
	  b = randomtb(&iseed);

	  rgb2hsv(r,g,b,&h,&s,&v);
	  hsv2rgb(h,s,v,&rn,&gn,&bn);

	  printf("%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n"
			,r,rn,g,gn,b,bn,h,s,v);

	  equaltb(r,rn,eps,&berror);
	  equaltb(g,gn,eps,&berror);
	  equaltb(b,bn,eps,&berror);
	  if( berror ) Error("error in computation r2h...");

	}

}


void equaltb(float a, float b, float eps, int *berror)

{
/*
	tests for nearly equality
*/

	float absab = abs(a-b);

	if( absab > eps ) {
	  printf("*** %f %f %f %f\n",a,b,absab,eps);
	  *berror = 1;
	}
}

/*

	subroutine test_ii

	implicit none

	real h,s,v
	real r,g,b
	logical h2r,r2h

	h2r = .false.
	h2r = .true.
	r2h = .true.
	r2h = .false.

	do while(h2r)
	  write(6,*) 'Enter h/s/v: '
	  read(5,*) h,s,v
	  call hsv2rgb(h,s,v,r,g,b)
	  write(6,*) 'r/g/b: ',r,g,b
	end do

	do while(r2h)
	  write(6,*) 'Enter r/g/b: '
	  read(5,*) r,g,b
	  call rgb2hsv(r,g,b,h,s,v)
	  write(6,*) 'h/s/v: ',h,s,v
	end do

	end

*/

int main( void )
{
	test_ct();
	/* test_ii(); */
	return 0;
}
