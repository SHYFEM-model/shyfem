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
	real caux(nlvdi)
	real dc,tcd,tce,eurpar,f

	integer iu,id,itmcon,idtcon,itstart
	save iu,id,itmcon,idtcon,itstart

        integer, save :: icall = 0

	if( icall < 0 ) return

!------------------------------------------------------------
! parameters
!------------------------------------------------------------

	iout_area = -1			!area considered outside, -1 for none
        bnoret = iout_area >= 0

	wsink = 0.
	wsink = 1.e-4
	wsink = 1.e-5
	wsink = 5.e-5

	wsink = 5.e-4		!sinking velocity [m/s]
	wsink = 0.		!if 0 -> do not use module
	rhos = 2500.		!density of sediments [kg/m**3]
	tce = 0.1		!critical threshold for erosion [N/m**2]
	tcd = 0.03		!critical threshold for deposition [N/m**2]
	eurpar = 1.e-3		!erosion parameter [kg/m**2/s]

! erosion rates e,d have units [kg/m**2/s]

	call get_timestep(dt)
	call getinfo(iunit)

	if( tce < tcd ) stop 'error stop simple_sedi: tce < tcd'

!------------------------------------------------------------
! initialization
!------------------------------------------------------------

        if( icall .eq. 0 ) then

          write(6,*) 'initialization of routine sedimt: ',wsink

	  if( wsink <= 0 ) icall = -1
	  if( icall < 0 ) return

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

	  call in_area(iout_area,inarea)	!sets up array inarea

          icall = 1

        end if

!------------------------------------------------------------
! is it time ?
!------------------------------------------------------------

        if( it .lt. itstart ) return

!------------------------------------------------------------
! sinking
!------------------------------------------------------------

          do k=1,nkn
	    lmax = ilhkv(k)
	    caux = 0
	    do l=1,lmax-1
              h = depnode(l,k,+1)
              vol = volnode(l,k,+1)
	      r = 0.
	      if( h .gt. 0. ) r = wsink/h
              conz = max(0.,cnv(l,k))
	      cnew = conz * exp(-r*dt)
	      dc = conz - cnew
	      caux(l) = caux(l) - dc
	      caux(l+1) = caux(l+1) + dc
	    end do
	    !call bottom_flux(k,f)
	    caux(lmax) = caux(lmax) + f
	    cnv(:,k) = cnv(:,k) + caux(:)
	    
	    !  conzs(k) = conzs(k) + vol*dc
	    !  conza(k) = conza(k) + h*dc
	    !  conzh(k) = conzh(k) + (h*dc)/rhos
          end do

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
! write accumulated bottom sediments
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

