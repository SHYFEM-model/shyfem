! This file is included in all .F90 files

#include "version.h"

#define PATH_MAX	255

#define stderr		0
#define stdout		6

! Handy for writing
#define STDOUT write(stdout,*)
#define STDERR write(stderr,*)
#define LEVEL0 STDERR
#define LEVEL1 STDERR '   ',
#define LEVEL2 STDERR '       ',
#define LEVEL3 STDERR '           ',
#define LEVEL4 STDERR '               ',
#define FATAL  STDERR 'FATAL ERROR: ',

#define LINE "------------------------------------------------------------------------"

! Shapes of variables
#define POINT           0
#define Z_SHAPE         1
#define T_SHAPE         2
#define XY_SHAPE        3
#define XYT_SHAPE       4
#define XYZT_SHAPE      5
#define OCET_SHAPE      6
#define SURFT_SHAPE     7
#define BOTT_SHAPE      8
#define G_SHAPE         9

! constants for average computations
#define INIT         0
#define MEAN         1
#define ACCUMULATE   10

! To avoid dividing by zero
#define SMALL 1e-8

! What precision will we use in this compilation
#define SINGLE
#undef  SINGLE

#ifdef SINGLE
#define REALTYPE real
#define MPI_REALTYPE	MPI_REAL
#define _ZERO_ 0.0
#define _ONE_  1.0
#else
#define REALTYPE double precision
#define MPI_REALTYPE	MPI_DOUBLE_PRECISION
#define _ZERO_ 0.0d0
#define _ONE_  1.0d0
#endif




