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
c
c*****************************************************************************

        subroutine heat2t(dt,dh,qs,qrad,albedo,ts,tsnew)

c computes new sea temperature
c
c radiation is positive if into the water

        implicit none

        real dt                 !time step
        real dh                 !layer depth
        real qs                 !solar radiation
        real qrad               !other radiation
	real albedo             !albedo to be used
        real ts                 !old temperature
        real tsnew              !new temperature

        real cw,row,ct
        real qseff

c--------------------------------------------------
c general constants
c--------------------------------------------------

	cw = 4186.		!specific heat (from Ali)
        cw = 3991.              !specific heat
        row = 1026.             !mean density
        ct = cw*row*dh          !heat capacity

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

	real temp	!water temperature
	real albedo	!albedo

        albedo = 0.06

        if( temp .lt. 4. ) then
          albedo = 0.5
          albedo = 0.75
        end if

	end

c*****************************************************************************

