!
! $Id: subcus.f,v 1.58 2010-03-08 17:46:45 georg Exp $
!
! simplified sedimentation module
!
! contents :
!
! subroutine simple_sedi	custom routines
!
! revision log :
!
! 03.02.2017	ggu	old routine copied from subcus.f
!
!******************************************************************

        subroutine simple_sedi

! simplified sedimentation module

	use mod_conz
	use levels
	use basin

        implicit none

	include 'femtime.h'

	real, save, allocatable :: conzs(:)
	real, save, allocatable :: conza(:)
	real, save, allocatable :: conzh(:)
	integer, save, allocatable :: inarea(:)

        integer ie,ii,k,lmax,l,ia
	integer iunit
	integer iout_area
        logical bnoret
        real vol,conz,perc,wsink,dt,sed,h,r,cnew,rhos
        double precision mass,masss
        real volnode,depnode
	real getpar

	integer iu,id,itmcon,idtcon,itstart
	save iu,id,itmcon,idtcon,itstart

        integer, save :: icall = 0

!------------------------------------------------------------
! parameters
!------------------------------------------------------------

	iout_area = -1			!area considered outside, -1 for none
        bnoret = iout_area >= 0

	wsink = 0.
	wsink = 1.e-4
	wsink = 1.e-5
	wsink = 5.e-5
	rhos = 2500.
	call get_timestep(dt)
	call getinfo(iunit)

!------------------------------------------------------------
! initialization
!------------------------------------------------------------

        if( icall .eq. 0 ) then

          write(6,*) 'initialization of routine sedimt: ',wsink

	  allocate(conzs(nkn))
	  allocate(conza(nkn))
	  allocate(conzh(nkn))
	  allocate(inarea(nkn))
	  conzs = 0.
	  conza = 0.
	  conzh = 0.
	  cnv = 0.

	  itstart = nint(getpar('tcust'))

          iu = 55
          itmcon = nint(getpar('itmcon'))
          idtcon = nint(getpar('idtcon'))
          call confop(iu,itmcon,idtcon,1,3,'set')

	  call in_area(iout_area,inarea)

          icall = 1

        end if

!------------------------------------------------------------
! is it time ?
!------------------------------------------------------------

        if( it .lt. itstart ) return

!------------------------------------------------------------
! sinking
!------------------------------------------------------------

	if( wsink .gt. 0. ) then
	  l = 1
          do k=1,nkn
              h = depnode(l,k,+1)
              vol = volnode(l,k,+1)
	      r = 0.
	      if( h .gt. 0. ) r = wsink/h
              conz = cnv(l,k)
              conz = max(0.,conz)
	      cnew = conz * exp(-r*dt)
              cnv(l,k) = cnew
              sed = vol * (conz-cnew) 
	      conzs(k) = conzs(k) + sed
              sed = h * (conz-cnew) 
	      conza(k) = conza(k) + sed
              sed = sed / rhos
	      conzh(k) = conzh(k) + sed
          end do
	end if

!------------------------------------------------------------
! total mass
!------------------------------------------------------------

        mass = 0.
        masss = 0.
        do k=1,nkn
            lmax = ilhkv(k)
            do l=1,lmax
              vol = volnode(l,k,+1)
              conz = cnv(l,k)
              mass = mass + vol*conz
            end do
	    masss = masss + conzs(k)
        end do

!------------------------------------------------------------
! write total mass
!------------------------------------------------------------

        write(6,*) 'sedimt: ',it,mass,masss,mass+masss
        write(iunit,*) 'sedimt: ',it,mass,masss,mass+masss

        id = 22       !for sediment -> [kg]
	call confil(iu,itmcon,idtcon,id,1,conzs)
        id = 23       !for sediment -> [kg/m**2]
	call confil(iu,itmcon,idtcon,id,1,conza)
        id = 24       !for sediment -> [m]
	call confil(iu,itmcon,idtcon,id,1,conzh)

!------------------------------------------------------------
! no return flow
!------------------------------------------------------------

        if( bnoret ) then
          do k=1,nkn
            if( inarea(k) .eq. 0 ) cnv(:,k) = 0.
          end do
        end if

!------------------------------------------------------------
! end of initialization
!------------------------------------------------------------

        end

!*****************************************************************

	subroutine in_area(iout_area,inarea)

	use basin

	implicit none

	integer iout_area
	integer inarea(nkn)

	integer ie,k,ii,ia

        inarea = 0

        do ie=1,nel
          ia = iarv(ie)
          if( ia .ne. iout_area ) then
              do ii=1,3
                k = nen3v(ii,ie)
                inarea(k) = 1
              end do
          end if
        end do

	end

!*****************************************************************

