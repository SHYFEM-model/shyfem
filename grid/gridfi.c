
/************************************************************************\
 *									*
 * gridfi.c - read/write routines for files				*
 *									*
 * Copyright (c) 1992-1994 by Georg Umgiesser				*
 *									*
 * Permission to use, copy, modify, and distribute this software	*
 * and its documentation for any purpose and without fee is hereby	*
 * granted, provided that the above copyright notice appear in all	*
 * copies and that both that copyright notice and this permission	*
 * notice appear in supporting documentation.				*
 *									*
 * This file is provided AS IS with no warranties of any kind.		*
 * The author shall have no liability with respect to the		*
 * infringement of copyrights, trade secrets or any patents by		*
 * this file or any part thereof.  In no event will the author		*
 * be liable for any lost revenue or profits or other special,		*
 * indirect and consequential damages.					*
 *									*
 * Comments and additions should be sent to the author:			*
 *									*
 *			Georg Umgiesser					*
 *			ISDGM/CNR					*
 *			S. Polo 1364					*
 *			30125 Venezia					*
 *			Italy						*
 *									*
 *			Tel.   : ++39-41-5216875			*
 *			Fax    : ++39-41-2602340			*
 *			E-Mail : georg@lagoon.isdgm.ve.cnr.it		*
 *									*
 * Revision History:							*
 * 16-Feb-2011: write to given file with OpOutFile			*
 * 07-May-1998: type is now integer                                     *
 * 12-Feb-1998: New routine stripgrd() -> strips .grd from file name    *
 *                is used in ReadFiles()                                *
 * 10-Oct-97: New routine SaveFile()                                    *
 * 06-Dec-95: In WriteStandard write for vector                         *
 * 04-Dec-95: ReadVect added, ReadStandard modified                     *
 * 10-Feb-95: closefile calls fclose only if file opened                *
 * 21-Oct-94: Changed introduced -> write only if file is changed       *
 * 08-Oct-94: reading/writing comments (Queuetable routines)            *
 * 07-May-94: new routines for opening and reading file                 *
 * 06-May-94: new file ff created for old read/write                    *
 *              filetype 0 is now the "official" filetype               *
 * 14-Apr-94: use GetActFileType to determine file type                 *
 * 13-Apr-94: use new hash routines                                     *
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "general.h"

#if __GUG_UNIX_

/*
#include <unistd.h>
#include <sys/file.h>
*/

#endif

#if __GUG_DOS_

/*
#include <io.h>
*/

#endif

#include "grid.h"
#include "gustd.h"
#include "keybd.h"

#include "args.h"
#include "queue.h"
#include "list.h"
#include "hash.h"
#include "gridhs.h"
#include "psgraph.h"


/**************************************************************************/

static int FileType = 0;
int GetActFileType( void ) { return FileType; }
void SetActFileType( int type ) { FileType = type; }

/**************************************************************************/

static FILE *FpFile=NULL;
static int   FpLine=0;

static int openfile( char *name , char *mode )

{
	FpLine=0;
	FpFile=fopen(name,mode);
	return FpFile ? 1 : 0;
}

static char *getnewline( void )

{
	char *s;

	s = getlin(FpFile);
	if( s )
		FpLine++;
	return s;
}

static int gettotlines( void ) { return FpLine; }
static void closefile( void ) { if(FpFile) fclose(FpFile); }

static char *stripgrd( char *s )

/* strips .grd from file name if found */

{
	int ls;
	char *t;

	ls = strlen(s);
	ls -= 4;

	if( ls > 0 ) {
		t = &s[ls];
		if( !strcmp(t,".grd") ) {
			t[0] = '\0';
		}
	}

	return s;
}

/**************************************************************************/

void ReadFiles( int argc , char *argv[] )

{
	char sfile[80];
	char *file,*s;
	int nodes,elems;
	int ftype;
	int errors=0;

	FileType = -1;

	if( argc <= OpArgc )
		argc = 0;

	for(;;)
	{
		if( argc == 0 ) {
			printf("Enter filename (<CR> to end) : ");
			file=getlin( stdin );
			s=strcpy(sfile,file);
		} else if( OpArgc < argc ) {
			s=strcpy(sfile,argv[OpArgc++]);
		} else {
			s=NULL;	/* everything read */
		}

		if( !s || !*s ) break;

		s = stripgrd(s);

		if( OpType == -1 ) {
			ftype = 0;
		} else
			ftype = OpType;

		switch(ftype) {
			case 0 :
				s=strcat(s,".grd");
				ReadStandard(s,HNN,HEL,HLI,HVC,CM);
				if(FileType == -1) FileType = ftype;
				nodes=0;
				elems=0;
				break;
			default:
				printf("File type not recognized. ");
				printf("Nothing read.\n");
				nodes=0;
				elems=0;
				break;
		}
		NTotNodes += nodes;
		NTotElems += elems;
	}

	if( errors ) {
		press_any_key();
	}

	printf("Working...\n");
}

void SaveFile( void )

{
	WriteStandard("save.grd",HNN,HEL,HLI,HVC,CM);
}

void WriteFiles( void )

{
	FILE *fp;
	char sfile[80];
	char *file,*s;
	int slength;

	/* files are always written in format 2 */

	if( FileType == 5 || FileType == 6 ) return;
	if( FileType == 4 ) FileType=3;

	ClosePS();

	if( OpOutFile ) {	/* output file name already given -> force */
		s=strcpy(sfile,OpOutFile);
		s=strcat(s,".grd");
		WriteStandard(s,HNN,HEL,HLI,HVC,CM);
		return;
	}

	for(;;) {
		if( Changed ) {
		    printf("Enter output filename (<CR> for no output) : ");
		    file=getlin( stdin );
		    if( ! file ) {
			s = NULL;
		    } else {
		    	s=strcpy(sfile,file);
		    }
		} else {
		    s=NULL;
		}

		if( !s || !*s ) return;

		/* test if file is already existing */

		slength=strlen(s);
		if(FileType == 3)
			s=strcat(s,".geo");
		else
			s=strcat(s,".grd");

		fp=fopen(s,"r");
		if( fp ) {
			fclose(fp);
			printf("File %s already existing. ",s);
			printf("Choose another name.\n");
			continue;
		}

		s[slength]='\0';
		if(FileType == 3)
			s=strcat(s,".dep");
		else
			s=strcat(s,".grd"); /* useless but not hurting */

		fp=fopen(s,"r");
		if( fp ) {
			fclose(fp);
			printf("File %s already existing. ",s);
			printf("Choose another name.\n");
			continue;
		}

		/* if we came here we are sure that file name is unique */

		s[slength]='\0';
		slength=strlen(s);
		s=strcat(s,".grd");
		WriteStandard(s,HNN,HEL,HLI,HVC,CM);

		break;
	}
}

void ReadStandard( char *fname , Hashtable_type HN , Hashtable_type HE
			      , Hashtable_type HL , Hashtable_type HV
				, Queuetable_type C )

{
	int comms=0,nodes=0,elems=0,lines=0,vects=0;
	int nodemax=0,elemmax=0,linemax=0,vectmax=0;
	int err=0;
	int error=FALSE;
	int narg,what,n;
	char *s,*t;

	if( openfile(fname,"r") )
		printf("Reading file %s\n",fname);
	else
		Warning2("ReadStandard : Cannot open file ",fname);

	while( (s=getnewline()) != NULL ) {
		t=firstchar(s);
		if( *t == '0' ) {
			t=savestring(s,-1);
			what=0;
		} else {
			narg = nargs(s);
			if(narg == 0) continue;
			t=readargs();
			what=atoi(t);
		}

		switch(what) {
		case 0 :			/* comment */
			InsertQueueTable(C,(void *)t);
			comms++;
			break;
		case 1 :			/* node */
			n = ReadNode(HN);
			if( n ) nodes++ ; else error=TRUE ;
			if( n > nodemax ) nodemax=n;
			break;
		case 2 :			/* element */
			n = ReadElem(HE);
			if( n ) elems++ ; else error=TRUE ;
			if( n > elemmax ) elemmax=n;
			break;
		case 3 :			/* line */
			n = ReadLine(HL);
			if( n ) lines++ ; else error=TRUE ;
			if( n > linemax ) linemax=n;
			break;
		case 4 :			/* vector */
			n = ReadVect(HV);
			if( n ) vects++ ; else error=TRUE ;
			if( n > vectmax ) vectmax=n;
			break;
		default:
			err++;
			printf("Line %d : ",gettotlines());
			printf("Shape %s not recognized\n",t);
			break;
		}
		if( error ) {
			printf("Read error in line %d :\n",gettotlines());
			error=FALSE;
			err++;
		}
	}

	NTotNodes += nodemax;
	NTotElems += elemmax;
	NTotLines += linemax;
	NTotVects += vectmax;

	printf("%d lines read\n",gettotlines());
	printf("Following shapes read :\n");
	if(comms) printf("Comments : %d ",comms);
	if(nodes) printf("Nodes : %d ",nodes);
	if(elems) printf("Elements : %d ",elems);
	if(lines) printf("Lines : %d ",lines);
	if(vects) printf("Vectors : %d ",vects);
	printf("\n");
	if( err ) {
		printf("Errors : %d\n",err);
		press_any_key();
	}

	closefile();
}

void WriteStandard( char *fname , Hashtable_type HN , Hashtable_type HE
			      , Hashtable_type HL , Hashtable_type HV
				, Queuetable_type C )

{
	FILE *fp;
	Node_type *pn;
	Elem_type *pe;
	Line_type *pl;
	Vect_type *pv;
	int nodes=0,elems=0,lines=0,vects=0,comments=0;
	int i,j;
	char *s;

        fp=fopen(fname,"w");
        if( fp )
                printf("Writing file %s\n",fname);
        else
                Error2("WriteStandard : Cannot open file ",fname);

	ResetQueueTable(C);
	while( (s=(char *)VisitQueueTable(C)) != NULL ) {
		fprintf(fp,"%s\n",s);
		comments++;
	}
        printf("Comments written : %d\n",comments);

	fprintf(fp,"\n");

        for(i=1;i<=NTotNodes;i++) {
          if( (pn=RetrieveByNodeNumber(HN,i)) != NULL ) {
                fprintf(fp,"1 %d %d %f %f"
                        ,pn->number
                        ,(int) pn->type
                        ,pn->coord.x
                        ,pn->coord.y
                        );
		if( pn->depth != NULLDEPTH )
			fprintf(fp," %f\n",pn->depth);
		else
			fprintf(fp,"\n");
                nodes++;
          }
        }
        printf("Nodes written : %d\n",nodes);

	fprintf(fp,"\n");

        for(i=1;i<=NTotElems;i++) {
          if( (pe=RetrieveByElemNumber(HE,i)) != NULL ) {
                fprintf(fp,"2 %d %d %d"
                        ,pe->number
                        ,(int) pe->type
			,pe->vertex
                        );

		for(j=0;j<pe->vertex;j++) {
			if( j%10 == 0 && pe->vertex > 3 )
				fprintf(fp,"\n");
			fprintf(fp," %d",pe->index[j]);
		}

                if( pe->depth != NULLDEPTH )
                        fprintf(fp," %f\n",pe->depth);
                else
                        fprintf(fp,"\n");
                elems++;
          }
        }
        printf("Elements written : %d\n",elems);

	fprintf(fp,"\n");

        for(i=1;i<=NTotLines;i++) {
          if( (pl=RetrieveByLineNumber(HL,i)) != NULL ) {
                fprintf(fp,"3 %d %d %d"
                        ,pl->number
                        ,(int) pl->type
                        ,pl->vertex
                        );

                for(j=0;j<pl->vertex;j++) {
                        if( j%10 == 0 )
                                fprintf(fp,"\n");
                        fprintf(fp," %d",pl->index[j]);
                }

                if( pl->depth != NULLDEPTH )
                        fprintf(fp," %f\n",pl->depth);
                else
                        fprintf(fp,"\n");
                lines++;
          }
        }
        printf("Lines written : %d\n",lines);

	fprintf(fp,"\n");

        for(i=1;i<=NTotVects;i++) {
          if( (pn=RetrieveByNodeNumber(HV,i)) != NULL ) {
                fprintf(fp,"4 %d %d %f %f"
                        ,pn->number
                        ,(int) pn->type
                        ,pn->coord.x
                        ,pn->coord.y
                        );
		pv = pn->extra;
                fprintf(fp," %d",pv->total);
		for(j=0;j<pv->total;j++) {
			if( j%3 == 0 && pv->total > 1 )
				fprintf(fp,"\n");
			fprintf(fp," %f %f",pv->speed[j],pv->dir[j]);
		}
		if( pv->total != 1 )
			fprintf(fp," %d\n",pv->actual+1);
		else
			fprintf(fp,"\n");
                vects++;
          }
        }
        printf("Vectors written : %d\n",vects);

	fclose(fp);
}

int ReadNode( Hashtable_type H )

{
	char *t;
	int number,ntype;
	Point c;
	float depth;
	Node_type *p;

	t=readargs();
	if( !t ) return 0;
	number = atoi(t);

	t=readargs();
	if( !t ) return 0;
	ntype = atoi(t);

	t=readargs();
	if( !t ) return 0;
	c.x = atof(t);

	t=readargs();
	if( !t ) return 0;
	c.y = atof(t);

	t=readargs();
	if( !t )
		depth = NULLDEPTH;
	else
		depth = atof(t);

	p=MakeNode(number+NTotNodes,ntype,&c);
	p->depth = depth;
	InsertByNodeNumber(H,p);

	return number;
}

int ReadElem( Hashtable_type H )

{
	char *t,*s;
	int number,ntype;
	int i,vertex;
	int *index;
	float depth;
	Elem_type *p;

	t=readargs();
	if( !t ) return 0;
	number = atoi(t);

	t=readargs();
	if( !t ) return 0;
	ntype = atoi(t);

	t=readargs();
	if( !t ) return 0;
	vertex = atoi(t);

	index = MakeIndex(vertex);

	i=0;
	while( i<vertex ) {
		t=readargs();
		if( !t ) {
			s=getnewline();
			if( !s ) return 0;
			initargs(s);
		} else {
			index[i++] = atoi(t) + NTotNodes;
		}
	}

	t=readargs();
	if( !t )
		depth = NULLDEPTH;
	else
		depth = atof(t);

	p = MakeElemWithIndex(number+NTotElems,ntype,vertex,index);
	p->depth = depth;
	InsertByElemNumber(H,p);

	return number;
}

int ReadLine( Hashtable_type H )

{
	char *t,*s;
	int number,ntype;
	int i,vertex;
	int *index;
	float depth;
	Line_type *p;

	t=readargs();
	if( !t ) return 0;
	number = atoi(t);

	t=readargs();
	if( !t ) return 0;
	ntype = atoi(t);

	t=readargs();
	if( !t ) return 0;
	vertex = atoi(t);

	index = MakeIndex(vertex);

	i=0;
	while( i<vertex ) {
		t=readargs();
		if( !t ) {
			s=getnewline();
			if( !s ) return 0;
			initargs(s);
		} else {
			index[i++] = atoi(t) + NTotNodes;
		}
	}

	t=readargs();
	if( !t )
		depth = NULLDEPTH;
	else
		depth = atof(t);

	p = MakeLineWithIndex(number+NTotLines,ntype,vertex,index);
	p->depth = depth;
	InsertByLineNumber(H,p);

	return number;
}

int even( int i )

{
	int j=i/2;

	return 2*j == i ? 1 : 0;
}
	
int ReadVect( Hashtable_type H )

{
	char *t,*s;
	int number,ntype;
	int total, actual;
	int i;
	Point c;
	float *speed, *dir;
	Node_type *pn;
	Vect_type *pv;

	t=readargs();
	if( !t ) return 0;
	number = atoi(t);

	t=readargs();
	if( !t ) return 0;
	ntype = atoi(t);

	t=readargs();
	if( !t ) return 0;
	c.x = atof(t);

	t=readargs();
	if( !t ) return 0;
	c.y = atof(t);

	t=readargs();
	while( !t ) {
		s=getnewline();
		if( !s ) return 0;
		initargs(s);
		t=readargs();
	}
	total = atoi(t);

	speed = MakeFloat(total);
	dir = MakeFloat(total);

	i=0;
	while( i<2*total ) {
		t=readargs();
		if( !t ) {
			s=getnewline();
			if( !s ) return 0;
			initargs(s);
		} else {
			if( even(i) )
			    speed[i/2] = atof(t);
			else
			    dir[i/2] = atof(t);
			i++;
		}
	}

	t=readargs();
	if( !t )
		actual = 1;
	else
		actual = atoi(t);

	pn=MakeNode(number+NTotVects,ntype,&c);
	InsertByNodeNumber(H,pn);
	pv=MakeVect(total,actual,speed,dir);
	pn->extra = pv;

	return number;
}
