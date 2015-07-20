c
c $Id: supver.f,v 1.4 2009-04-20 08:13:29 georg Exp $
c
c Vertical velocities
c
c revision log :
c
c 26.03.2010    ggu     HACK for bwater -> compute always
c
c******************************************************************

        subroutine vertical

c computes vertical velocities
c
c from sp256w in new3di.F
c
c velocities are computed on S/T points (top and bottom of layer)
c bottom velocity of the whole column is assumed to be 0
c -> maybe change this
c
c computes volume differences and from these computes vertical
c velocities at every time step so that the computed velocities
c satisfy the continuity equation for every single time step
c
c wlnv is computed horizontally at a node and vertically
c it is at the center of the layer -> there are nlv velocities
c computed
c
c b,c are 1/m, (phi is dimensionless)
c aj is m**2
c utlnv... is m**2/s
c dvol is in m**3/s
c vv is m**2 (area)
c
c wlnv is first used to accumulate volume difference -> dvol
c at the end it receives the vertical velocity
c
c wlnv (dvol)   aux array for volume difference
c vv            aux array for area
c
c written on 27.08.91 by ggu  (from scratch)
c 14.08.1998    ggu     w = 0 at open boundary nodes
c 20.08.1998    ggu     some documentation

	use mod_plot2d
	use mod_plot3d
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels
	use basin

        implicit none

c parameters
        include 'param.h'
c common




c local
	logical byaron
        integer k,ie,ii,kk,l
        integer ilevel
	integer inwater
        real aj,wbot,wtop,ff

c initialize

	byaron = .true.
	byaron = .false.

        do k=1,nkn
          do l=0,nlv
            wauxv(l,k)=0.
            wlnv(l,k) = 0.
          end do
        end do

c compute difference of velocities for each layer
c
c f(ii) > 0 ==> flux into node ii

	inwater = 0

        do ie=1,nel
         !if( bwater(ie) ) then           !FIXME	!not working
	  inwater = inwater + 1
          aj=4.*ev(10,ie)               !area of triangle / 3
          ilevel = ilhv(ie)
          do l=1,ilevel
            do ii=1,3
                kk=nen3v(ii,ie)
                ff = utlnv(l,ie)*ev(ii+3,ie) + vtlnv(l,ie)*ev(ii+6,ie)
	if( byaron .and. kk .eq. 2088 .and. l .eq. 8 ) then
		ff = ff + 10./(6.*3.*aj)
	end if
                wlnv(l,kk) = wlnv(l,kk) + 3. * aj * ff
                wauxv(l,kk)=wauxv(l,kk)+aj
            end do
          end do
         !end if
        end do

	!write(6,*) '******** inwater = ',inwater

c from vel difference get absolute velocity (w_bottom = 0)
c       -> wlnv(nlv,k) is already in place !
c       -> wlnv(nlv,k) = 0 + wlnv(nlv,k)
c w of bottom of last layer must be 0 ! -> shift everything up
c wlnv(nlv,k) is always 0
c
c dividing dvol(m**3/s) by area (wauxv) gives vertical velocity

        do k=1,nkn
          wbot = 0.
          do l=nlv,1,-1
            wtop = wlnv(l,k)
            wlnv(l,k) = wbot
            wbot = wbot + wtop
          end do
          wlnv(0,k) = wbot
        end do

        do k=1,nkn
          do l=1,nlv
            if( wauxv(l,k) .gt. 0. ) then
              wlnv(l-1,k) = wlnv(l-1,k) / wauxv(l,k)
            end if
          end do
        end do

c set w to zero at open boundary nodes (new 14.08.1998)
c
c FIXME -> only for ibtyp = 1,2 !!!!

c        do k=1,nkn
c          if( inodv(k) .gt. 0 ) then    !open boundary node
c            do l=0,nlv
c               wlnv(l,k) = 0.
c            end do
c          end if
c        end do

        return
        end

c******************************************************************

