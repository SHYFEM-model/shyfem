
c***************************************************************
c
c look in subcst for initialization
c
c cstinit	M_INIT		start of hp/ht
c cstcheck	M_CHECK		after read of STR file
c cstsetup	M_SETUP		after setup of basic arrays (ev,...)
c
c look in subnsh for read, write and others
c
c prilog	M_PRINT		after everything has been setup (before loop)
c pritst	M_TEST		only for debug
c
c dobefor	M_BEFOR		beginning of each time loop
c doafter	M_AFTER		end of each time loop
c
c nlsh2d	M_READ		read in STR file
c
c***************************************************************
c
c modules still to transform:
c
c wrouta
c resid
c rmsvel
c adjust_chezy (bottom friction)
c prwnds
c prarea (checy)
c prclos
c proxy, prlgr
c prbnds
c pripar
c biocos
c locous, locspc
c
c***************************************************************

	subroutine modules(mode)

c handles module framework

	implicit none

	integer mode		!mode of call

c------------------------------------------------------------
c modules
c------------------------------------------------------------

	call mod_ext(mode)		!extra nodes
	call mod_ets(mode)		!extra time series nodes
c	call mod_box(mode)		!boxes
c	call mod_flx(mode)		!fluxes through sections
c	call mod_vol(mode)		!volumes in areas

c------------------------------------------------------------
c end of routine
c------------------------------------------------------------

	end

c***************************************************************

