
!***************************************************************
!
! look in subcst for initialization
!
! cstinit	M_INIT		start of hp/ht
! cstcheck	M_CHECK		after read of STR file
! cstsetup	M_SETUP		after setup of basic arrays (ev,...)
!
! look in subnsh for read, write and others
!
! prilog	M_PRINT		after everything has been setup (before loop)
! pritst	M_TEST		only for debug
!
! dobefor	M_BEFOR		beginning of each time loop
! doafter	M_AFTER		end of each time loop
!
! nlsh2d	M_READ		read in STR file
!
!***************************************************************
!
! modules still to transform:
!
! wrouta
! resid
! rmsvel
! adjust_chezy (bottom friction)
! prwnds
! prarea (checy)
! prclos
! proxy, prlgr
! prbnds
! pripar
! biocos
! locous, locspc
!
!***************************************************************
!---------------------------------------------------------------------
        module modls
!---------------------------------------------------------------------
        contains
!---------------------------------------------------------------------

	subroutine modules(mode)

! handles module framework

        use extra
        use ets_admin

	implicit none

	integer mode		!mode of call

!------------------------------------------------------------
! modules
!------------------------------------------------------------

	call mod_ext(mode)		!extra nodes
	call mod_ets(mode)		!extra time series nodes
!	call mod_box(mode)		!boxes
!	call mod_flx(mode)		!fluxes through sections
!	call mod_vol(mode)		!volumes in areas

!------------------------------------------------------------
! end of routine
!------------------------------------------------------------

	end

!***************************************************************

!---------------------------------------------------------------------
        end module modls
!---------------------------------------------------------------------
