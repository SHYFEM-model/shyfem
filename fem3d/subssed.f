
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017-2020  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

! simplified sedimentation module
!
! contents :
!
! subroutine simple_sedi	custom routines
!
! revision log :
!
! 03.02.2017	ggu	old routine copied from subcus.f
! 13.02.2017	ggu	changed VERS_7_5_23
! 31.03.2017	ggu	changed VERS_7_5_24
! 13.04.2017	ggu	changed VERS_7_5_25
! 09.05.2017	ggu	some bugs fixed
! 03.04.2018	ggu	changed VERS_7_5_43
! 31.08.2018	ggu	changed VERS_7_5_49
! 14.02.2019	ggu	changed VERS_7_5_56
! 16.02.2019	ggu	changed VERS_7_5_60
! 20.02.2020	ggu	new routine compute_bottom_flux()
!
! notes :
!
! in order to run the module set issedi=1 in the STR file, $para section
! output frequency is according to itmcon, idtcon
! files written are with extension .ssed.shy
!
!******************************************************************

!==================================================================
	module simple_sediments
!==================================================================

	implicit none

	integer, private, save :: nkn_ssedi = 0

	real, save, allocatable :: conzs(:)	!bottom sediment [kg]
	real, save, allocatable :: conza(:)	!bottom sediment [kg/m**2]
	real, save, allocatable :: conzh(:)	!bottom sediment [m]
	real, save, allocatable :: sedflux(:)	!sediment flux [kg/m**2/s]
	integer, save, allocatable :: inarea(:)	!0 if area out of basin

!------------------------------------------------------------------
! sediment flux is positive from sediment into water column
!------------------------------------------------------------------

	logical, save :: bssedi = .false.	!is running?
	integer, save :: issedi = 0	!1 -> use module (set in STR file)

	double precision, save :: da_out(5)	!index for output file

!------------------------------------------------------------------
! user defined parameters - please customize
!------------------------------------------------------------------

	logical, save :: bonlys = .true. !only resuspend settled sediments
	integer, save :: iout_area = -1	!area considered outside, -1 for none

	real, save :: wsink = 5.e-4	!sinking velocity [m/s]
	real, save :: rhos = 2500.	!density of sediments [kg/m**3]
	real, save :: tce = 0.1		!critical threshold erosion [N/m**2]
	real, save :: tcd = 0.03	!critical threshold deposition [N/m**2]
	real, save :: eurpar = 1.e-3	!erosion parameter [kg/m**2/s]

!==================================================================
	contains
!==================================================================

	subroutine simple_sediments_init(nkn)

        integer ncs
        integer nkn
        integer nlv

        if( nkn == nkn_ssedi ) return

        if( nkn_ssedi > 0 ) then
	  deallocate(conzs)
	  deallocate(conza)
	  deallocate(conzh)
	  deallocate(sedflux)
	  deallocate(inarea)
        end if

        nkn_ssedi = nkn

        if( nkn == 0 ) return

	allocate(conzs(nkn))
	allocate(conza(nkn))
	allocate(conzh(nkn))
	allocate(sedflux(nkn))
	allocate(inarea(nkn))

	conzs = 0.
	conza = 0.
	conzh = 0.
	sedflux = 0.
	inarea = 0

	end subroutine simple_sediments_init

!==================================================================
	end module simple_sediments
!==================================================================

        subroutine simple_sedi

! simplified sedimentation module

	use mod_conz
	use levels
	use basin
	use simple_sediments

        implicit none

        integer ie,ii,k,lmax,l,ia
	integer iunit
        logical bnoret
        real vol,conz,perc,dt,sed,h,r,cnew
        double precision mass,masss
        double precision dtime,dtime0
        real volnode,depnode
	real getpar
	real caux(nlvdi)
	real taubot(nkn)
	real fflux(nkn)
	real dc,f,tau,alpha
	real cmin,cmax
	character*20 aline

	integer iu,id,itmcon,idtcon,itstart
	save iu,id,itmcon,idtcon,itstart

        integer, save :: icall = 0

	if( icall < 0 ) return

!------------------------------------------------------------
! parameters
!------------------------------------------------------------

        bnoret = iout_area >= 0		!set concentrations out of domain to 0

	call get_timestep(dt)
	call getinfo(iunit)
	call get_act_dtime(dtime)
	call get_timeline(dtime,aline)

	if( tce < tcd ) stop 'error stop simple_sedi: tce < tcd'

!------------------------------------------------------------
! initialization
!------------------------------------------------------------

        if( icall .eq. 0 ) then

          issedi = nint(getpar('issedi'))
          if( issedi .le. 0 ) icall = -1
          if( icall .le. -1 ) return
          icall = 1

          write(6,*) 'initialization of routine sedimt: ',issedi

	  if( iconz /= 1 ) then
	    write(6,*) 'cannot run simple sediment module'
	    write(6,*) 'iconz must be == 1'
	    stop 'error stop simple_sedi: iconz /= 1'
	  end if

	  call simple_sediments_init(nkn)

	  call get_first_dtime(dtime0)
	  call simple_sedi_init_output
	  call simple_sedi_write_output(dtime0)

	  call in_area(iout_area,inarea)	!sets up array inarea

          icall = 1
	  bssedi = .true.

          write(6,*) 'finished initialization of routine sedimt'
        end if

!------------------------------------------------------------
! is it time ?
!------------------------------------------------------------

        !if( it .lt. itstart ) return

!------------------------------------------------------------
! sinking
!------------------------------------------------------------

	  !call bottom_stress(taubot)
	  call compute_bottom_flux(dt,fflux)

 	  cmax = maxval(cnv)

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
            h = depnode(lmax,k,+1)
            vol = volnode(lmax,k,+1)
	    tau = taubot(k)
	    r = dt/h
	    !call bottom_flux(k,tau,cnv(lmax,k),r,alpha,f) !f is sediment flux

	    f = fflux(k)
	    if( bonlys .and. f*dt > conza(k) ) then	!limit erosion
	      f = conza(k) / dt
	    end if

	    sedflux(k) = f
	    dc = f * dt / h
	    caux(lmax) = caux(lmax) + dc
	    cnv(:,k) = cnv(:,k) + caux(:)

	    conzs(k) = conzs(k) - vol*dc	! [kg]
	    conza(k) = conza(k) - h*dc		! [kg/m**2]
	    conzh(k) = conzh(k) - (h*dc)/rhos	! [m]

            !write(6,*) conzs(k),conza(k),conzh(k),f,'simple_sed_b_s'
	    !if( k == 100 ) write(6,*) k,tau,cnv(lmax,k),f,dc
          end do

	cmin = minval(cnv)
	cmax = maxval(cnv)
	!write(6,*) cmin,cmax

	where( cnv < 0. ) cnv = 0.

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

        !write(6,*) 'sedimt: ',dtime,mass,masss,mass+masss
        !write(iunit,*) 'sedimt: ',dtime,mass,masss,mass+masss
        write(iunit,1200) ' sedimt: ',aline,mass,masss,mass+masss
 1200	format(a,a,4e14.6)

	call convert_bottom_sediments(.true.)

!------------------------------------------------------------
! write accumulated bottom sediments
!------------------------------------------------------------

	call simple_sedi_write_output(dtime)

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

! computes areas that are considered inside basin

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

	subroutine compute_bottom_flux(dt,fflux)

	use basin
	use levels
	use simple_sediments
	use mod_conz
	use mod_layer_thickness	
	use evgeom

	implicit none

	real dt
	real fflux(nkn)

	integer k,ia,lmax,ie,ii
	real alpha,f,conz
	real area
	real tau,r,h
	real taubot(nkn)
	real aux(nkn)

	call bottom_stress(taubot)
	aux = 0.
	fflux = 0.

	do ie=1,nel
	  area = 4.*ev(10,ie)
	  ia = iarv(ie)
	  lmax = ilhv(ie)
	  h = hdenv(lmax,ie)
	  r = dt/h
	  tce = 0.1
	  !if( ia == 0 ) tce = 10.0	!FIXME
	  do ii=1,3
	    k = nen3v(ii,ie)
	    conz = cnv(lmax,k)
	    tau = taubot(k)
	    call bottom_flux(k,tau,conz,r,alpha,f)
	    fflux(k) = fflux(k) + area*f
	    aux(k) = aux(k) + area
	  end do
	end do

	where( aux > 0. ) fflux = fflux / aux

	end

!*****************************************************************

	subroutine bottom_flux(k,tau,conz,r,alpha,f)

! computes fluxes between bottom and water column

	use simple_sediments
	use mod_conz

	implicit none

	integer k	!node
	real tau	!bottom stress
	real conz	!concentration in last layer
	real alpha	!flux factor [dimensionless]
	real r		!factor for exponential deposition (dt/h)
	real f		!sediment flux, positive into water column [kg/m**2/s]

	real dc

	if( tau < tcd ) then			!deposition (f negative)
	  alpha = - ( 1. - tau/tcd )
	  dc = conz*(exp(alpha*r*wsink)-1.)
	  f = dc / r
	  !f = alpha * wsink * conz
	else if( tau > tce ) then		!erosion (f positive)
	  alpha = ( tau/tce - 1. )
	  f = alpha * eurpar
	else					!nothing
	  f = 0.
	end if

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine simple_sedi_init_output

! this opens two files, one for bottom sediments (2D), and one for
! concentrations in water column (3D)

	use simple_sediments

	implicit none

	integer, save :: nvar2d = 3
	integer, save :: nvar3d = 1
	integer id
	logical has_output_d

	da_out = 0

        call init_output_d('itmcon','idtcon',da_out)
        if( has_output_d(da_out) ) then
          call shyfem_init_scalar_file('ssed',nvar2d,.true.,id)
          da_out(4) = id
          call shyfem_init_scalar_file('csed',nvar3d,.false.,id)
          da_out(5) = id
        end if

	end

!*****************************************************************

	subroutine simple_sedi_write_output(dtime)

	use levels
	use mod_conz, only: cnv
	use simple_sediments

	implicit none

	double precision dtime

	integer id,idcbase
	logical next_output_d

        if( .not. next_output_d(da_out) ) return

	idcbase = 850

        id = nint(da_out(4))

        call shy_write_scalar_record2d(id,dtime,idcbase+1,conzs) ! [kg]
        call shy_write_scalar_record2d(id,dtime,idcbase+2,conza) ! [kg/m**2]
        call shy_write_scalar_record2d(id,dtime,idcbase+3,conzh) ! [m]
	call shy_sync(id)

        id = nint(da_out(5))

        call shy_write_scalar_record(id,dtime,idcbase,nlvdi,cnv)

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine get_sediment_values(k,flux,conz)

	use levels
	use mod_conz
	use simple_sediments

	implicit none

	integer k	 !node
	real flux	 !sediment flux at node k [kg/m**2/s]
	real conz(nlvdi) !sediment concentration in water column [kg/m**3]

	if( .not. bssedi ) then
	  write(6,*) 'bssedi: ',bssedi
	  stop 'error stop get_sediment_values: sediments not running'
	end if

	flux = sedflux(k)
	conz(:) = cnv(:,k)

	end

!*****************************************************************

	subroutine convert_bottom_sediments(bcheck)

! converts conza to conzs and conzh
! if bcheck is true, only checks if consistent

	use basin
	use levels
	use simple_sediments

	implicit none

	logical bcheck

	integer k,lmax
	real h,vol
	real smax,hmax,vmax
	real auxs(nkn),auxh(nkn),auxv(nkn)

	real depnode,volnode

        do k=1,nkn
	  lmax = ilhkv(k)
          h = depnode(lmax,k,+1)
          vol = volnode(lmax,k,+1)

	  !conzs(k) = conzs(k) - vol*dc	! [kg]
	  !conza(k) = conza(k) - h*dc		! [kg/m**2]
	  !conzh(k) = conzh(k) - (h*dc)/rhos	! [m]

	  auxs(k) = conza(k)*vol/h
	  auxh(k) = conza(k)/rhos
	  auxv(k) = vol/h

        end do

	if( bcheck ) then
	  smax = maxval(abs(auxs-conzs))
	  hmax = maxval(abs(auxh-conzh))
	  vmax = maxval(abs(auxv))
	  write(778,*) smax,hmax,vmax
	  call flush(778)
	else
	  conzs = auxs
	  conzh = auxh
	end if

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

        subroutine write_restart_ssedi(iunit)
	use basin
	use simple_sediments
        implicit none
        integer iunit
	integer, save :: nstate = 1
	write(iunit) nstate,nkn
	write(iunit) conza
        end

        subroutine skip_restart_ssedi(iunit)
	use basin
	use simple_sediments
        implicit none
        integer iunit
	read(iunit)
	read(iunit)
        end

        subroutine read_restart_ssedi(iunit)
	use basin
	use simple_sediments
        implicit none
        integer iunit
	integer, save :: nstate = 1
	integer nstate_aux,nkn_aux
	read(iunit) nstate_aux,nkn_aux
	if( nstate /= nstate_aux ) goto 99
	if( nkn /= nkn_aux ) goto 99
	read(iunit) conza
	return
   99	continue
	write(6,*) 'nstate: ',nstate,nstate_aux
	write(6,*) 'nkn   : ',nkn,nkn_aux
	stop 'error stop read_restart_ssedi: incompatible params'
        end

!*****************************************************************

