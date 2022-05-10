
!--------------------------------------------------------------------------
!
!    Copyright (C) by the GOTM-team under the GNU Public License - www.gnu.org
!
!    This file is part of SHYFEM. (m)
!
!--------------------------------------------------------------------------

! This file is include in all .F90 files and contains very important
! definitions. Infact GOTM will not compile when this file is not
! in a correct format.
! KBK 20000220

#include "version.h"

#define PATH_MAX	255

!#define stderr		0		!GGU
#define stderr		6
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

! Shapes for variables
#define POINT           0
#define Z_SHAPE         1
#define T_SHAPE         2
#define XY_SHAPE        3
#define XYT_SHAPE       4
#define XYZT_SHAPE      5

#define RAWBINARY       0
#define ASCII           1
#define NETCDF          2
#define GRADS           3
#define OPENDX          4

! For easier reading
#define READING 0
#define WRITING 1

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

! non-local fluxes
#undef NONLOCAL

! KPP turbulence model
#define KPP_SHEAR
#define KPP_INTERNAL_WAVE
#define KPP_CONVEC
#undef KPP_DDMIX
#undef KPP_TWOPOINT_REF
#define KPP_IP_FC
#undef KPP_CLIP_GS
#define KPP_SALINITY

! next are needed by debug GGU (with a trailing _ they are not used)
! in order to use them delete the trailing underscore
! files affected: 
!	cppdefs.h 
!	diff_face.F90 
!	dissipationeq.F90
!	tkeeq.F90 
!	turbulence.F90

#define DEBUG_GGU_
#define GGU_ZERO_DISEQ
#define GGU_ZERO_TKEEQ


