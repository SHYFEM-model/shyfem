
/************************************************************************\
 *
 *    Copyright (C) 1995,1997,1999-2000,2002,2011  Georg Umgiesser
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
 * mesh.c - automatic meshing routines
 *
 * revision log :
 *
 * 25.07.1995	ggu	routines written from scratch
 * 08.10.1997	ggu	insert given points in domain delimited by line
 * ...		ggu	introduced pointer to background grid
 * ...		ggu	use new mesh type
 * 15.10.1997	ggu	MarkOuterElements: use InClosedLine() for inside check
 * ...		ggu	use MakeCoordsFromLine() to make coordinate list
 * 16.10.1997	ggu	call CopyBoundaryLine() in main instead RefineBoundary()
 * ...		ggu	new call to RecoverBoundaryNodes() introduced
 * 12.11.1997	ggu	new routine MarkOuterNodes() to get rid of unused nodes
 * ...		ggu	MarkOuterElements() and MarkOuterNodes() are called
 * ...		ggu	after internal refining -> internal nodes can be
 * ...		ggu	inserted even if outer element has to be changed
 * ...		ggu	-> whole routine gets more flexible
 * ...		ggu	RecoverBoundaryNodes() is called after internal refining
 * ...		ggu	in InterpolRes() the standard resolution (OpResolution)
 * ...		ggu	is returned if no background triangle is found
 * 25.11.1997	ggu	new names for .grd files -> M_*.grd
 * 05.05.1999	ggu	check for zero area in SetFEMParam
 * 03.02.2000	ggu	different algorithm to compute resolution in bckgrid
 * 14.08.2002	ggu	InterpolRes uses cubic to compute resolution
 * 05.12.2011	ggu	adapt for changes in AreaConvex
 *
\************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
#include <stdlib.h>
#include <strings.h>
#include <ctype.h>
#include <math.h>
*/


/*
#include "general.h"
#include "gustd.h"
#include "fund.h"
*/

#include "general.h"
#include "gustd.h"

#include "fund.h"
#include "list.h"
#include "nlist.h"
#include "queue.h"
#include "hash.h"
#include "stack.h"

#include "mesh.h"
#include "meshut.h"
#include "meshfi.h"
#include "meshcv.h"
#include "meshck.h"
#include "meshge.h"
#include "meshin.h"
#include "meshbd.h"
#include "meshhs.h"
#include "meshop.h"
#include "meshty.h"


int   NTotNodes  = 0;       /* number of total nodes */
int   NTotElems  = 1;       /* number of total elements (1 is not used) */
int   NTotLines  = 0;       /* number of total lines */

Hashtable_type HNN;
Hashtable_type HEL;
Hashtable_type HLI;
ListTable      NEL;             /* list of new elements */
ListTable      BCK;             /* list of background grid elements */
QueueTable      CM;             /* list of comments */


void MakeMaxiTriang( NodeList hull );
int MarkExternalElements( NodeList hull );
int MarkOuterElements( void );
int MarkOuterNodes( void );
float Resol( float x , float y );
void SetResolution( NodeList list );
void SetFEMParam( Elem_type *pe );
float InterpolRes( Elem_type *pe , float x , float y );




int main(int argc, char *argv[])

{
	int ntotal,external;
	NodeList hull, intern, boundary, given;

	HNN=MakeHashTable();
	HEL=MakeHashTable();
	HLI=MakeHashTable();
	NEL=MakeListTable();
	CM=MakeQueueTable();

 	SetOptions(argc,argv);

	ReadFiles(argc,argv);
	/* CheckBoundary(); */

	ntotal = CheckInput();
	hull = MakeNodeList(ntotal);
	intern = MakeNodeList(ntotal);

	printf("Making convex hull... %d\n",hull->count);
	ConvexHull(hull,intern);
/*
	CheckConvex(hull,intern);
	PrintNodeList("hull",hull);
	PrintNodeList("intern",intern);
*/

	printf("Making maxi elements...\n");
	MakeMaxiTriang(hull);
	CheckNeibor(-1);

	printf("Inserting convex hull... %d\n",hull->count);
	InsertBoundaryNodes(hull);
	CheckCircumCircle();
	CheckNeibor(-2);

	WriteAll("M_hull.grd",hull);
	CheckNeibor(-3);

	printf("Inserting internal boundary points... %d\n",intern->count);
	InsertBoundaryNodes(intern);
	CheckCircumCircle();
	CheckNeibor(-4);
	WriteAll("M_orgbound.grd",NULL);

	SetResolution(hull);
	CopyBoundaryLine();

	printf("Recovering boundary lines 1...\n");
	RecoverBoundary();
	CheckCircumCircle();
	CheckNeibor(-43);
	WriteAll("M_bndrecover.grd",NULL);

	printf("Refining boundary points 1...\n");
	boundary = RefineBoundary();

	if( boundary ) {
	    printf("Inserting new boundary points... %d\n",boundary->count);
	    InsertBoundaryNodes(boundary);
	    CheckCircumCircle();
	    CheckCircumCircleProperty();
	    CheckNeibor(-5);
	}

	TestVersion();

	printf("Marking external elements...");
	external=MarkExternalElements( hull );
	printf(" %d / %d\n",external,NTotElems);
	WriteGrd("M_finebound.grd");

/*
	printf("Marking outer elements...");
	external=MarkOuterElements();
	printf(" %d / %d\n",external,NTotElems);
	WriteGrd("M_test.grd");
*/

	FreeNodeList(hull);
	FreeNodeList(intern);
	FreeNodeList(boundary);

	given = GivenNodes();
	printf("Inserting internal given points... %d\n",given->count);
	InsertNodes(given);
	FreeNodeList(given);
	CheckCircumCircle();
	CheckNeibor(-44);
	WriteAll("M_given.grd",NULL);

	printf("Recovering boundary lines 2...\n");
	RecoverBoundary();
	CheckCircumCircle();
	CheckNeibor(-45);
	WriteAll("M_intrecover.grd",NULL);

	CheckArea();
	printf("Inserting internal points...\n");
	InsertInternalNodes();
	CheckCircumCircle();
	CheckCircumCircleProperty();
	WriteGrd("M_insert.grd");

	TestVersion();

	CheckArea();
	printf("Refining internal points... %f\n",OpAspect);
	RefineInternalNodes();
	CheckArea();
	CheckCircumCircle();
	CheckCircumCircleProperty();
	WriteGrd("M_refine.grd");

	CheckArea();
	printf("Recovering boundary lines 3...\n");
	RecoverBoundary();
	printf("Recovering fault lines...\n");
	RecoverInternalFault();
	CheckCircumCircle();
	CheckNeibor(-48);
	WriteAll("M_intrecover2.grd",NULL);

	printf("Marking outer elements...");
	external=MarkOuterElements();
	printf(" %d / %d\n",external,NTotElems);
	printf("Marking outer nodes...");
	external=MarkOuterNodes();
	printf(" %d / %d\n",external,NTotNodes);
	WriteGrd("M_test.grd");

	TestVersion();

	CheckArea();
	printf("Smoothing internal points... %f\n",OpSmoothOmega);
	SmoothInternalNodes();
	CheckArea();
	WriteGrd("final.grd");

	return 0;
}


void MakeMaxiTriang( NodeList hull )

/*\
 *  Make maxi triangle to include all points
\*/

{
	int ntotal;
	int i,imin;
	int n1,n2,n3,n4;
	int nx1,nx2,nx3,nx4;
	int nex1,nex2,nex3,nex4;
	int nem1,nem2;
	float xmin,ymin,xmax,ymax;
	float x,y,dx,dy;
	float fact=2.1;	/* factor for construction of exterior elements */
	Node_type *pn;

	ntotal = hull->count;

	for(i=0;i<ntotal;i++) {
	    pn = RetrieveByNodeNumber(HNN,hull->index[i]);
	    if( !IsNtype(pn, N_EXTERNAL ) ) break;
	}
	if( i == ntotal )
		Error("MakeMaxiTriang: No node available");

	pn = RetrieveByNodeNumber(HNN,hull->index[i]);
	xmin = pn->coord.x;
	ymin = pn->coord.y;
	xmax = xmin;
	ymax = ymin;
	imin = i+1;

	for(i=imin;i<ntotal;i++) {
		pn = RetrieveByNodeNumber(HNN,hull->index[i]);
	        if( IsNtype(pn, N_EXTERNAL ) ) continue;
		x = pn->coord.x;
		y = pn->coord.y;
		if( x > xmax ) xmax = x;
		if( y > ymax ) ymax = y;
		if( x < xmin ) xmin = x;
		if( y < ymin ) ymin = y;
	}

	dx = xmax - xmin;
	dy = ymax - ymin;

	xmin -= dx/2.;
	ymin -= dy/2.;
	xmax += dx/2.;
	ymax += dy/2.;

	/* xmin,... are coordinates of maxi rectangle */

	/* make nodes from vertices of maxi rectangle */

	n1=InsertNewNode(N_EXTERNAL,xmin,ymin);
	n2=InsertNewNode(N_EXTERNAL,xmax,ymin);
	n3=InsertNewNode(N_EXTERNAL,xmax,ymax);
	n4=InsertNewNode(N_EXTERNAL,xmin,ymax);

	/* make auxiliary nodes for external triangles */

	nx1=InsertNewNode(N_EXTERNAL,(xmin+xmax)/2.,ymin-fact*dy);
	nx2=InsertNewNode(N_EXTERNAL,xmax+fact*dx,(ymin+ymax)/2.);
	nx3=InsertNewNode(N_EXTERNAL,(xmin+xmax)/2.,ymax+fact*dy);
	nx4=InsertNewNode(N_EXTERNAL,xmin-fact*dx,(ymin+ymax)/2.);

	/* make two maxi triangles */

	nem1 = InsertNewElem( E_INTERNAL, n1, n2, n3
				, NB_UNKNOWN, NB_UNKNOWN, NB_UNKNOWN);
	nem2 = InsertNewElem( E_INTERNAL, n1, n3, n4
				, NB_UNKNOWN, NB_UNKNOWN, NB_UNKNOWN);
	
	/* make four external triangles (auxiliary) */

	nex1 = InsertNewElem( E_EXTERNAL, n1, nx1, n2
				, NB_UNKNOWN, NB_UNKNOWN, NB_UNKNOWN);
	nex2 = InsertNewElem( E_EXTERNAL, n2, nx2, n3
				, NB_UNKNOWN, NB_UNKNOWN, NB_UNKNOWN);
	nex3 = InsertNewElem( E_EXTERNAL, n3, nx3, n4
				, NB_UNKNOWN, NB_UNKNOWN, NB_UNKNOWN);
	nex4 = InsertNewElem( E_EXTERNAL, n4, nx4, n1
				, NB_UNKNOWN, NB_UNKNOWN, NB_UNKNOWN);

	/* update neighbor information */

	UpdateNeiborByNumber(nem1,nex2,nem2,nex1);
	UpdateNeiborByNumber(nem2,nex3,nex4,nem1);
	UpdateNeiborByNumber(nex1,0,nem1,0);
	UpdateNeiborByNumber(nex2,0,nem1,0);
	UpdateNeiborByNumber(nex3,0,nem2,0);
	UpdateNeiborByNumber(nex4,0,nem2,0);
	
	/* insert new elements (except extern) in list */

	InsertListTable(NEL,(void *) RetrieveByElemNumber(HEL,nem1));
	InsertListTable(NEL,(void *) RetrieveByElemNumber(HEL,nem2));
}


int MarkExternalElements( NodeList hull )

{
	int i,ncount;
	int external=0;
	int node;
	float *xe, *ye;
	float xm,ym;
	Elem_type *pe;
	Node_type *pn;

	ncount = hull->count;

	xe = (float *) malloc( ncount * sizeof(float) );
	ye = (float *) malloc( ncount * sizeof(float) );
	if( !xe || !ye )
		Error("MarkExternalElements: Cannot allocate coordinate list");

	for(i=0;i<ncount;i++) {
		node = hull->index[i];
		pn = RetrieveByNodeNumber(HNN,node);
		xe[i] = pn->coord.x;
		ye[i] = pn->coord.y;
	}

	ResetHashTable(HEL);
	while( (pe=VisitHashTableE(HEL)) != NULL ) {
		MakeCM(pe,&xm,&ym);
		if( !InConvex(ncount,xe,ye,xm,ym) ) {
			SetEtype(pe, E_EXTERNAL );
			external++;
		}
	}

	free(xe);
	free(ye);

	return external;
}

int MarkOuterElements( void )

{
	int ext;
	int loops;
	int external=0;
	float *xe, *ye;
	float xm,ym;
	Elem_type *pe;
	Line_type *pl;

        ResetHashTable(HLI);
        while( (pl = VisitHashTableL(HLI)) != NULL ) {
	    if( !IsLtype(pl,L_EXTERNAL_REF) && !IsLtype(pl,L_INTERNAL_REF) )
			continue;

	    MakeCoordsFromLine(pl,&xe,&ye);

	    ResetHashTable(HEL);
	    while( (pe=VisitHashTableE(HEL)) != NULL ) {
		if( IsEtype(pe, E_EXTERNAL ) ) continue;
		MakeCM(pe,&xm,&ym);
		loops = InClosedLine(pl->vertex,xe,ye,xm,ym);
		ext = 0;
		if( IsLtype(pl,L_EXTERNAL_REF) ) {
		    if( loops == 0 ) ext = 1;
		} else {
		    if( loops != 0 ) ext = 1;
		}
		if( ext ) {
		    SetEtype(pe, E_EXTERNAL );
                    external++;
		}
	    }

	    free(xe);
	    free(ye);
	}

	return external;
}

int MarkOuterNodes( void )

{
	int ext;
	int loops;
	int external=0;
	float *xe, *ye;
	float xm,ym;
	Node_type *pn;
	Line_type *pl;

        ResetHashTable(HLI);
        while( (pl = VisitHashTableL(HLI)) != NULL ) {
	    if( !IsLtype(pl,L_EXTERNAL_REF) && !IsLtype(pl,L_INTERNAL_REF) )
			continue;

	    MakeCoordsFromLine(pl,&xe,&ye);

	    ResetHashTable(HNN);
	    while( (pn=VisitHashTableN(HNN)) != NULL ) {
		if( IsNtype(pn, N_EXTERNAL ) ) continue;
		if( IsNtype(pn, N_BOUNDARY ) ) continue;
		if( IsNtype(pn, N_ADDBOUND ) ) continue;
		xm = pn->coord.x;
		ym = pn->coord.y;
		loops = InClosedLine(pl->vertex,xe,ye,xm,ym);
		ext = 0;
		if( IsLtype(pl,L_EXTERNAL_REF) ) {
		    if( loops == 0 ) ext = 1;
		} else {
		    if( loops != 0 ) ext = 1;
		}
		if( ext ) {
		    SetNtype(pn, N_EXTERNAL );
                    external++;
		}
	    }

	    free(xe);
	    free(ye);
	}

	return external;
}

float expon( float x , float y )

{
	float r;

	r = 0.3*exp(0.2*x) + 0.*y ;
	return r;
}

float elipsref( float x , float y )

{
	float r;

	r = 0.3*(ABS(x-y)/2. + 1.);	/* for elipsis with refinement */
	r = 0.5*(ABS(x-y)/2. + 1.);	/* for elipsis with refinement */
	return r;
}

float elips( float x , float y )

{
	float r;

	r = 30.-0.*x-0.*y ;	/* also 3 */
	return r;
}

float constant( float x , float y )

{
	float r;

	r = 1.-0.*x-0.*y ;
	return r;
}

float joe( float x , float y )

{
	float r=0.;
	float max=120.;
	float fac=3.;
	float xmin=5.99;

	r = 0.*y;
	if(x>xmin) {
		r += max*exp(-fac*(x-xmin));
	}
	r = max + 1. - r;
	r /= max;

	return r;
}

float Resol( float x , float y )

{
	Elem_type *pe;
	float aux;
	int debug = 0;

	if( OpBackGround >= 0 ) {
		pe = FindXYElement(BCK,x,y);
		aux = InterpolRes(pe,x,y);
		if( debug && pe ) {
		printf("interpol: %f %f %f %f %f\n",x,y,OpResolution,aux,
				aux/OpResolution);
		}
		return aux;
	} else {
		return OpResolution;
	}
/*
	if( OpResolution > 0. ) {
		return OpResolution;
	} else {
		return elipsref(x,y);
	}
*/
}

void SetResolution( NodeList list )

/*\
 *  the factor 0.5 when computing OpResolution is empirical
 *  and accounts for the fact that triangles are not equilateral
\*/

{
	int i;
	int ie=0;
	float area,length;
	float resdef;
	Elem_type *pe;
	Node_type *pn;

	if( OpIntern ) {
		area = AreaConvex(list);
		OpResolution = (4./sqrt(27)) * area/OpIntern;
		printf("SetResolution: (1) %f %f\n",area,OpResolution);
	} else if( OpBoundary ) {
		length = LengthConvex(list);
		OpResolution = (1./sqrt(3.)) * length/OpBoundary;
		OpResolution *= OpResolution;
		printf("SetResolution: (2) %f %f\n",length,OpResolution);
	} else if( OpResolution > 0. ) {
		OpResolution *= OpResolution;
		printf("SetResolution: (3) %f\n",OpResolution);
	}

	resdef = 200.*AreaConvex(list);

	if( OpBackGround >= 0 ) {
		ResetListTable(BCK);
		while( ( pe = (Elem_type *) VisitListTable(BCK) ) != NULL ) {
		    for(i=0;i<3;i++) {
			pn = RetrieveByNodeNumber(HNN,pe->index[i]);
	    		if( IsNtype(pn, N_EXTERNAL ) ) {
			  if( OpResolution > 0. ) {
/*	new: do only after interpolation...
			    pn->depth = OpResolution / (pn->depth*pn->depth);
*/
			  } else {
			    pn->depth = resdef;
			  }
			  SetNtype(pn, N_NONE );
			}
		    }
		    ie++;
		    SetFEMParam(pe);
		}
		ResetListTable(BCK);
		while( ( pe = (Elem_type *) VisitListTable(BCK) ) != NULL ) {
		    for(i=0;i<3;i++) {
			pn = RetrieveByNodeNumber(HNN,pe->index[i]);
			SetNtype(pn, N_EXTERNAL );
		    }
		}
		printf("SetResolution: (background) %d\n",ie);
	}

	if( OpResolution <= 0. ) {
		OpResolution = resdef;
	}
}

void SetFEMParam( Elem_type *pe )

{
	int i,in,ip;
	int debug = 0;
	float f=0.;
	float x[3],y[3],v[3];
	Node_type *pn;
	Backg_type *pb;

	pb = (Backg_type *) malloc( sizeof(Backg_type) );
	if( !pb )
		Error("SetFEMParam: Cannot allocate background mesh");
	((Extra_E_type *)pe->extra)->backg = pb;

	for(i=0;i<3;i++) {
		pn = RetrieveByNodeNumber(HNN,pe->index[i]);
		x[i] = pn->coord.x;
		y[i] = pn->coord.y;
		v[i] = pn->depth;
	}

	for(i=0;i<3;i++) {
		in = (i+1)%3;
		ip = (i+2)%3;
		pb->a[i] = x[in]*y[ip] - x[ip]*y[in];
		pb->b[i] = y[in] - y[ip];
		pb->c[i] = x[ip] - x[in];
		pb->v[i] = v[i];
		f += pb->a[i];
	}

	if( f == 0.0 || debug ) {
	  printf("%d %d %f\n",pe->number,pe->vertex,f);
	  for(i=0;i<3;i++) {
	    pn = RetrieveByNodeNumber(HNN,pe->index[i]);
	    printf("%f %f %f %d\n",x[i],y[i],v[i],pn->number);
	  }
	  if( ! debug ) {
	    Error2("SetFEMParam: Element with zero area: ",itos(pe->number-1));
		/* FIXME -> we have to subtract 1 from element number */
	  }
	}

	pb->fr = 1./f;
}

float InterpolRes( Elem_type *pe , float x , float y )

{
	int i;
	float f;
	float a=0.;
	float vmin,vmax;
	float dv,an;
	Backg_type *pb;

	if( !pe )
		return OpResolution;

	pb = ((Extra_E_type *)pe->extra)->backg;

	vmin = vmax = pb->v[0];

	for(i=0;i<3;i++) {
		f = pb->a[i] + x*pb->b[i] + y*pb->c[i];
		a += pb->v[i] * f * pb->fr;
		if( pb->v[i] > vmax ) vmax = pb->v[i];
		if( pb->v[i] < vmin ) vmin = pb->v[i];
	}

	/* new code */

	dv = vmax - vmin;
	if( dv > 0. ) {
	  x = (a-vmin)/dv;
	  an = vmin + x*x*x * dv;
	} else {
	  an = a;
	}
	/*
	printf("InterpolRes: %f  %f %f   %f %f\n",x,vmin,vmax,a,an);
	*/

	a = an;

	return OpResolution/a;
/*
	return OpResolution/(a*a);
	return a;
*/
}

