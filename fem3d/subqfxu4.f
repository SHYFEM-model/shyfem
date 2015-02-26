c
c $Id: subqfxu4.f,v 1.4 2009-11-18 16:50:37 georg Exp $
c
c heat flux module (temperature module)
c
c contents :
c
c subroutine heat2t(dt,dh,qs,qrad,albedo,ts,tsnew)
c               computes new sea temperature
c
c revision log :
c
c 16.08.2004    ggu     heat2t copied from subqfxt.f
c 23.03.2006    ggu     changed time step to real
c 11.11.2009    ggu     new routine make_albedo(), pass albedo to heat2t
c 01.06.2011    ggu     use constant albedo
c
c*****************************************************************************

        subroutine heat2t(dt,dh,qs,qrad,albedo,ts,tsnew)

c computes new sea temperature
c
c radiation is positive if into the water

        implicit none

	include 'subqfxm.h'

        real dt                 !time step
        real dh                 !layer depth
        real qs                 !solar radiation
        real qrad               !other radiation
	real albedo             !albedo to be used
        real ts                 !old temperature
        real tsnew              !new temperature

        real ct
        real qseff

c--------------------------------------------------
c general constants
c--------------------------------------------------

        ct = cpw*rhow*dh          !heat capacity

c--------------------------------------------------
c new temperature
c--------------------------------------------------

        qseff = qs * ( 1. - albedo )
        tsnew = ts + (qseff+qrad)*dt/ct

c--------------------------------------------------
c end of routine
c--------------------------------------------------

        end

c*****************************************************************************

        subroutine make_albedo(temp,albedo)

        implicit none

        real temp       !water temperature
        real albedo     !albedo

        real albed0,albed4
        save albed0,albed4

        double precision dgetpar

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

c*****************************************************************************

