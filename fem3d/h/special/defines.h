
c Chao like velocity profile at inlet -----------------------------------------

#ifndef		BCHAO
#define 	BCHAO				0
#endif

#if BCHAO
#ifndef		BCHAO_IMPOSE_TRANSPORT
#define		BCHAO_IMPOSE_TRANSPORT		1
#endif
#ifndef		BCHAO_ADD_MOMENTUM
#define		BCHAO_ADD_MOMENTUM		0
#endif
#endif BCHAO

c advective terms discretization ----------------------------------------------

#ifndef		ADVECTIVE_NODE_CENTRED
#define		ADVECTIVE_NODE_CENTRED		1
#endif
#ifndef		ADVECTIVE_ELEM_CENTRED
#define		ADVECTIVE_ELEM_CENTRED		0
#endif

c no baroclinicity in certain areas -------------------------------------------

#ifndef		BAROC_AREA0
#define		BAROC_AREA0			0
#endif

c adjustment of vertical austausch coefficient --------------------------------

#ifndef		VERT_AUST_ADJUST
#define 	VERT_AUST_ADJUST		0
#endif

#if VERT_AUST_ADJUST
#ifndef		AFACT
#define		AFACT				1000.
#endif
#endif VERT_AUST_ADJUST

c adjustment of vertical diffusion coefficient --------------------------------

#ifndef		VERT_DIFF_ADJUST
#define 	VERT_DIFF_ADJUST		0
#endif

#if VERT_DIFF_ADJUST
#ifndef		DFACT
#define		DFACT				1000.
#endif
#endif VERT_DIFF_ADJUST

c adjustment of discharge -----------------------------------------------------

#ifndef		ADJUST_DISCHARGE
#define		ADJUST_DISCHARGE		0
#endif

#if ADJUST_DISCHARGE
#ifndef		TADJUST
#define		TADJUST				21600.
#endif
#endif

c end -------------------------------------------------------------------------

