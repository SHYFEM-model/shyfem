
/************************************************************************\
 *
 *    Copyright (C) 2003  Georg Umgiesser
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
 * colorutil.h - color utilities for color space conversions
 *
 * revision log :
 *
 * 18.08.2003	ggu	adapted for use in psgraph
 *
\************************************************************************/


#ifndef __GUH_COLORUTIL_
#define __GUH_COLORUTIL_


void rgb2cmy(float r ,float g ,float b ,float *c ,float *m ,float *y);
void cmy2rgb(float c ,float m ,float y ,float *r ,float *g ,float *b);
void rgb2hsv(float r ,float g ,float b ,float *h ,float *s ,float *v);
void hsv2rgb(float h, float s, float v, float *r, float *g, float *b);


#endif /* __GUH_COLORUTIL_ */
