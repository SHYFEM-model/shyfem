
/************************************************************************\
 *									*
 * grid_fp.h - function prototypes for grid                             *
 *									*
 * Copyright (c) 1992-1994 by Georg Umgiesser				*
 *									*
 * see grid.c for copying information					*
 *									*
 * Revision History:                                                    *
 * 16-Feb-2011: in MakeElem/Line() pass also type			*
 * 13-May-2003: new include menu.h and prototype ExecuteMenuCommand()	*
 * 02-Apr-1998: no SetNewCommand() & Buttons, new ExitEventLoop()       *
 * 09-Feb-1998: ActArgument eliminated, new functions GfZoom, GfShow    *
 * 14-Oct-97: new routines GfRemoveElement, GfRemoveLine                *
 *            new routines GetAct/SetAct... for incapsulation           *
 * 13-Oct-97: new routines AddUseN(), DeleteUseN(), GetUseN()           * 
 *               DeleteLineWithNodes(), DeleteElemWithNodes()           *
 * 10-Oct-97: new prototypes GfSave, GfUnifyNode, SubstituteNode,       *
 *                           SaveFile                                   *
 * 04-Dec-95: new prototypes for Vector routines                        *
 * 11-Mar-95: Prototype for CheckConnections not needed                 *
 *            prototypes for AreaElement, InvertIndex introduced        *
 * 06-Oct-94: new prototypes for GetColor... (no ColTab)                *
 * 07-May-94: new prototypes for ReadNode...                            *
 * 13-Apr-94: splitted up in df, ty, fp, ex                             *
 * 06-Apr-94: copyright notice added to file				*
 * ..-...-92: routines written from scratch				*
 *									*
\************************************************************************/


#ifndef __GUH_GRID_FP_
#define __GUH_GRID_FP_


#include <stdio.h>

/**************************************************************/
/**************************************************************/
/******************** type definitions ************************/
/**************************************************************/
/**************************************************************/

#include "fund.h"
#include "list.h"
#include "queue.h"
#include "hash.h"
#include "grid_ty.h"

#include "menu.h"

/**************************************************************/
/**************************************************************/
/******************** function protoypes **********************/
/**************************************************************/
/**************************************************************/

void GfCancel( void );
void GfRefresh( void );
void GfPrint( void );
void GfSave( void );
void GfExit( void );

void GfZoomWindow( void );
void GfZoomIn( void );
void GfZoomOut( void );
void GfMove( void );
void GfTotalView( void );
void GfMoveRelative( float dx, float dy );

void GfShowNode( void );
void GfShowElement( void );
void GfShowLine( void );
void GfShowVect( void );

void GfMakeNode( void );
void GfDelNode( void );
void GfMoveNode( void );
void GfUnifyNode( void );

void GfMakeElement( void );
void GfDelElement( void );
void GfRemoveElement( void );

void GfMakeLine( void );
void GfDelLine( void );
void GfRemoveLine( void );
void GfJoinLine( void );
void GfSplitLine( void );
void GfDelNodeLine( void );
void GfRemoveNodeLine( void );
void GfInsertNodeLine( void );

void GfDelVect( void );
void GfChangeVect( void );

void GfChangeDepth( void );
void GfChangeType( void );

void InitializeFunctions( void );

void ExitEventLoop( void );
void LoopForInput( void );

int GetActNode( void ) ;
void SetActNode( int node ) ;
int GetActElem( void ) ;
void SetActElem( int elem ) ;
int GetActLine( void ) ;
void SetActLine( int line ) ;
int GetActVect( void ) ;
void SetActVect( int vect ) ;

Rect *GetActPlotWindow( void ) ;

void PlotFieldInput( int horiz , int verti , int button );
void MenuFieldInput( int horiz , int verti );
void KeyboardInput( int c );

void ExecuteMenuCommand( FP fp );

void PlotAll( void );
void RedrawAll( void );
int ResizeWindow( int width , int height );
void WriteToMesWindow( void );
void WriteToComWindow( void );

Node_type *FindClosestNode( Hashtable_type H , float x , float y );
Line_type *FindClosestLine( Hashtable_type HL , Hashtable_type HN
				, float x , float y );
Elem_type *FindElemToNode( Hashtable_type H , Elem_type *p , int node );
Line_type *FindLineToNode( Hashtable_type H , Line_type *p , int node );
Node_type *FindNode( Hashtable_type H , int node );
Elem_type *FindElem( Hashtable_type H , int elem );
Line_type *FindLine( Hashtable_type H , int line );
Elem_type *FindElemToPoint( Hashtable_type HE , Hashtable_type HN
			, float x , float y );

int InElement( Hashtable_type H , Elem_type *pe , float x , float y );
int InConvex( int n , float *xe , float *ye , float x , float y );

int GetUseN( Node_type *pn );
void AddUseN( Node_type *pn );
void DeleteUseN( Node_type *pn );
void AddUseE( Hashtable_type H , Elem_type *pe );
void DeleteUseE( Hashtable_type H , Elem_type *pe );
void AddUseL( Hashtable_type H , Line_type *pl );
void DeleteUseL( Hashtable_type H , Line_type *pl );

void UnActive( void );
void EvidenceNone( void );
void MakeNodeActive( int node );
void MakeElemActive( int elem );
void MakeLineActive( int line );
void MakeVectActive( int line );

void MakeNewCenter( Rect *gp , float *x , float *y , float fact );

void ZoomInOut(Rect *gp , float x , float y , float fact );
void MoveRelative(Rect *gp , float dx , float dy );
void MoveToPoint( float x , float y );

void MakeMidPoint( Line_type *p , float *x , float *y );
void MakeGravityPoint( Elem_type *p , float *x , float *y );

void Working( void );

void DeleteElemWithNodes( Elem_type *pe );
void DeleteLineWithNodes( Line_type *pl );
void SplitLine( Hashtable_type H , Line_type *pl , int node );
void JoinLine( Hashtable_type H , Line_type *p1 ,  Line_type *p2 , int node );
void DelNodeLine( Hashtable_type H , Line_type *pl , int node );
void InsertNodeLine( Hashtable_type H , Line_type *pl , int node );

float MakeDepthFromNodes( Hashtable_type H , Elem_type *p );

void SubstituteNode( Node_type *pna , Node_type *pnu );

/* gridfi */

int GetActFileType( void );
void SetActFileType( int type );

void ReadFiles( int argc , char *argv[] );
void SaveFile( void );
void WriteFiles( void );

void ReadStandard( char *fname , Hashtable_type HN , Hashtable_type HE
                              , Hashtable_type HL , Hashtable_type HV
				, Queuetable_type C);
void WriteStandard( char *fname , Hashtable_type HN , Hashtable_type HE
                              , Hashtable_type HL , Hashtable_type HV
				, Queuetable_type C);

int ReadNode( Hashtable_type H );
int ReadElem( Hashtable_type H );
int ReadLine( Hashtable_type H );
int ReadVect( Hashtable_type H );

/* gridwi */

void MakeTotalWindow( Rect *gbp );
void MakeMenuWindow( Rect *gbp );
void MakePlotWindow( Rect *gbp );
void MakeMessageWindow( Rect *gbp );
void MakeCommandWindow( Rect *gbp );

/* gridpl */

void EvidenceNode( int node , int color );
void EvidenceElem( int elem , int color );
void EvidenceLine( int line , int color );
void EvidenceVect( int line , int color );
void GetScreenCoord( float x , float y , int *horiz , int *verti );
void GetPlotCoord( int horiz , int verti , float *x , float *y );
void GetMenuCoord( int horiz , int verti , float *x , float *y );
void ScalePlotWindow( int xmin , int ymin , int xmax , int ymax , Rect *gb );
void ScaleFactor( Hashtable_type HV );
void PlotRect( Rect *r );
void PlotShadeRect( Rect *r , int col1 , int col2 );
void ErasePlot( void );

void PlotPoints( Hashtable_type HN );
void PlotPoint( Node_type *pn );
void PlotElements( Hashtable_type HE , Hashtable_type HN );
void PlotElem( Hashtable_type H , Elem_type *p );
void FillElem( Hashtable_type H , Elem_type *p , int color );
void PlotLines( Hashtable_type HL , Hashtable_type HN );
void PlotLine( Hashtable_type H , Line_type *p );
void PlotSegment( float x1 , float y1 , float x2 , float y2 , int color );
void PlotVectors( Hashtable_type HV );
void PlotVect( Node_type *p );

/* gridge */

void GetElemMinMax( Hashtable_type HE , Hashtable_type HN , Rect *r );
void InitMinMax( Rect *r );
void SetMinMax( Rect *r );
float GetDepthMinMax( void );
void GetNodeMinMax( Hashtable_type H , Rect *r );
float GetAverLat( Hashtable_type H );
int IsLatLon( Hashtable_type H );
int IsDegenerateRect( Rect *r );
void CopyRect( Rect *d , Rect *s );
void MakeRect ( Point *p , Rect *r );
void MakeRectFromPoints ( Point *p1 , Point *p2 , Rect *r );
void TriMinMax( Hashtable_type H , int *index , Rect *r );
void PolyMinMax( Hashtable_type H , Elem_type *pe , Rect *r );
void PolyMinMaxIndex( Hashtable_type H , int nvert , int *index , Rect *r );
void AdjustBounds( Rect *bounds , Rect *new );
float angle( Hashtable_type H , int k1 , int k2 , int k3 );

/* gridut */

void PrintNodes( Hashtable_type H );
void PrintElems( Hashtable_type H );
void PrintLines( Hashtable_type H );
void PrintLineS( Hashtable_type HL, Hashtable_type HN, int line, int invert);
Node_type *NewNode( void );
Vect_type *MakeVect( int total , int actual , float *s , float *d );
void DeleteVect( Node_type *p );
void ChangeVect( Vect_type *p );

Node_type *MakeNode( int n , int type , Point *c );
void DeleteNode( Node_type *p );
Elem_type *MakeElem( int n , int type , int *c , int vertex );
Elem_type *MakeElemWithIndex( int n , int ntype , int vertex , int *index );
void DeleteElem( Elem_type *p );
Line_type *MakeLine( int n , int type , int *c , int vertex );
Line_type *MakeLineWithIndex( int n , int ntype , int vertex , int *index );
void DeleteLine( Line_type *p );

int *MakeIndex( int vertex );
void InvertIndex( int *index , int nvert );
float *MakeFloat( int total );

float AreaElement( Hashtable_type H , Elem_type *pe );
float Dist2Node( Hashtable_type H , int node1 , int node2 );

/* gridhs */


/* gridop */

void SetOptions(int argc, char *argv[]);
int GetColorTabSize( void );
void SetColors( int coltab );
int GetTypeColor( int type );
int GetDepthColor( float value );
int GetRandomColor( int value );
int GetVelColor( float value );

/* gridps */

void WritePS( Rect *gp , Hashtable_type HN , Hashtable_type HE
                              , Hashtable_type HL , Queuetable_type C );
void ClosePS( void );


#endif
