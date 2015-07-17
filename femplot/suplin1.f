c
c $Id: suplin1.f,v 1.1 2009-11-18 17:16:00 georg Exp $
c
c utilitiy routines for section plot (velocities)
c
c revision log :
c
c 13.10.2009    ggu     routines written from scratch
c
c*******************************************************************

	subroutine prepare_vel(pp3)

	use mod_plot2d !COMMON_GGU_SUBST
	use mod_plot3d !COMMON_GGU_SUBST
	use mod_hydro_print !COMMON_GGU_SUBST
	use levels, only : nlvdi,nlv !COMMON_GGU_SUBST
	use basin, only : nkn,nel,ngr,mbw !COMMON_GGU_SUBST

	implicit none

	include 'param.h'

	real pp3(nlvdi,nkn)

COMMON_GGU_DELETED	include 'nbasin.h'
COMMON_GGU_DELETED	include 'nlevel.h'

COMMON_GGU_DELETED	include 'plot_aux.h'
COMMON_GGU_DELETED	include 'plot_aux_3d.h'

COMMON_GGU_DELETED	include 'hydro_print.h'

	integer k,l
	real u,v,w
	real href

	call vertical		!compute wlnv from utlnv,vtlnv

	href = 0.
	call mkht3(nlvdi,het3v,href)

	call make_vel_from_tra(het3v)
	call vel_to_node
	
	do k=1,nkn
	  do l=1,nlv
	    u = uprv(l,k)
	    v = vprv(l,k)
	    w = wprv(l,k)
	    pp3(l,k) = sqrt( u*u + v*v + w*w )
	  end do
	end do

	end

c*******************************************************************

	subroutine make_vel_from_tra(het3v)

	use mod_hydro_vel !COMMON_GGU_SUBST
	use mod_hydro !COMMON_GGU_SUBST
	use levels !COMMON_GGU_SUBST
	use basin, only : nkn,nel,ngr,mbw !COMMON_GGU_SUBST

	implicit none

	include 'param.h'

        real het3v(nlvdi,nel)

COMMON_GGU_DELETED	include 'nbasin.h'

COMMON_GGU_DELETED	include 'hydro.h'
COMMON_GGU_DELETED	include 'hydro_vel.h'
COMMON_GGU_DELETED	include 'nlevel.h'
COMMON_GGU_DELETED	include 'levels.h'

	integer ie,l,lmax
	real h,rh

	do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
	    h = het3v(l,ie)
	    if( h .le. 0. ) goto 99
	    rh = 1. / h
	    ulnv(l,ie) = utlnv(l,ie) * rh
	    vlnv(l,ie) = vtlnv(l,ie) * rh
	  end do
	end do

	return
   99	continue
	write(6,*) 'ie,l,h: ',ie,l,h
	stop 'error stop make_vel_from_tra: zero depth'
	end

c*******************************************************************

	subroutine vel_to_node

	use mod_hydro_print !COMMON_GGU_SUBST
	use mod_hydro_vel !COMMON_GGU_SUBST
	use evgeom !COMMON_GGU_SUBST
	use levels !COMMON_GGU_SUBST
	use basin !COMMON_GGU_SUBST

	implicit none

	include 'param.h'

COMMON_GGU_DELETED	include 'nlevel.h'

COMMON_GGU_DELETED	include 'basin.h'
COMMON_GGU_DELETED	include 'levels.h'

COMMON_GGU_DELETED	include 'hydro_vel.h'

COMMON_GGU_DELETED	include 'hydro_print.h'

COMMON_GGU_DELETED	include 'ev.h'

	integer ie,ii,k,l,lmax
	real aj

        do k=1,nkn
          do l=1,nlv
            uprv(l,k)=0.
            vprv(l,k)=0.
            wprv(l,k)=0.
          end do
        end do

        do ie=1,nel
          aj=ev(10,ie)
	  lmax = ilhv(ie)
          do ii=1,3
            k=nen3v(ii,ie)
	    do l=1,lmax
              wprv(l,k)=wprv(l,k)+aj
              uprv(l,k)=uprv(l,k)+aj*ulnv(l,ie)
              vprv(l,k)=vprv(l,k)+aj*vlnv(l,ie)
            end do
	  end do
	end do

        do k=1,nkn
	  do l=1,nlv
            if(wprv(l,k).gt.0.) then
              uprv(l,k)=uprv(l,k)/wprv(l,k)
              vprv(l,k)=vprv(l,k)/wprv(l,k)
            end if
          end do
        end do

        do k=1,nkn
          do l=1,nlv
            wprv(l,k)=0.5*(wlnv(l,k)+wlnv(l-1,k))
          end do
        end do

	end

c*******************************************************************

