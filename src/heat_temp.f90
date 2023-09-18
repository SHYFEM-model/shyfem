!
! $Id: subqfxu4.f,v 1.4 2009-11-18 16:50:37 georg Exp $
!
! heat flux module (temperature module)
!
! contents :
!
! subroutine heat2t(dt,dh,qs,qrad,albedo,ts,tsnew)
!               computes new sea temperature
!
! revision log :
!
! 16.08.2004    ggu     heat2t copied from subqfxt.f
! 23.03.2006    ggu     changed time step to double precision
! 11.11.2009    ggu     new routine make_albedo(), pass albedo to heat2t
! 01.06.2011    ggu     use constant albedo
! 04.05.2016    ccf     do not pass albedo into heat2t
!
!*****************************************************************************
!-----------------------------------------------------------------------------
        module heat_temp
!-----------------------------------------------------------------------------
        contains
!-----------------------------------------------------------------------------

        subroutine heat2t(dt,dh,qs,qrad,ts,tsnew)

! computes new sea temperature
!
! radiation is positive if into the water

        implicit none

	include 'subqfxm.h'

        double precision dt                 !time step
        double precision dh                 !layer depth
        double precision qs                 !solar radiation corrected
        double precision qrad               !other radiation
        double precision ts                 !old temperature
        double precision tsnew              !new temperature

        double precision ct
        double precision qseff

!--------------------------------------------------
! general constants
!--------------------------------------------------

        ct = cpw*rhow*dh          !heat capacity

!--------------------------------------------------
! new temperature
!--------------------------------------------------

        tsnew = ts + (qs+qrad)*dt/ct

!--------------------------------------------------
! end of routine
!--------------------------------------------------

        end

!*****************************************************************************

        subroutine make_albedo(temp,albedo)

        use para

        implicit none

        double precision temp       !water temperature
        double precision albedo     !albedo

        double precision albed0,albed4
        save albed0,albed4

        integer icall
        save icall
        data icall /0/

        if( icall .eq. 0 ) then
          albed0 = dgetpar('albedo')
          albed4 = dgetpar('albed4')
          icall = 1
        end if

        albedo = albed0

        if( temp .lt. 4. ) then
          albedo = albed4
          !albedo = 0.5
          !albedo = 0.75
        end if

        end

!*****************************************************************************

!-----------------------------------------------------------------------------
        end module heat_temp
!-----------------------------------------------------------------------------
