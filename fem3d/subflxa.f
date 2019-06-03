
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
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

c subroutines for computing discharge / flux
c
c contents :
c
c subroutine inflxa
c subroutine rdflxa
c subroutine ckflxa
c subroutine prflxa
c subroutine tsflxa
c
c subroutine wrflxa				write of flux data
c
c subroutine flxscs(n,kflux,iflux,az,fluxes)	flux through sections
c subroutine flxsec(n,kflux,iflux,az,fluxes)	flux through section
c
c subroutine flxini				initializes flux routines
c subroutine flx_init(kfluxm,kflux,nsect,iflux)	sets up array iflux
c subroutine flxinf(m,kflux,iflux)		sets up one info structure
c function igtnsc(k1,k2)			gets number of internal section
c
c revision log :
c
c 30.04.1998	ggu	newly written routines (subpor deleted)
c 07.05.1998	ggu	check nrdveci on return for error
c 08.05.1998	ggu	restructured with new comodity routines
c 13.09.1999	ggu	type of node computed in own routine flxtype
c 19.11.1999	ggu	iskadj into sublin
c 20.01.2000	ggu	old routines substituted, new routine extrsect
c 20.01.2000	ggu	common block /dimdim/ eliminated
c 20.01.2000	ggu	common block /dimdim/ eliminated
c 01.09.2002	ggu	ggu99 -> bug in flx routines (how to reproduce?)
c 26.05.2003	ggu	in flxnov substituted a,b with b,c
c 26.05.2003	ggu	new routine make_fluxes (for lagrangian)
c 10.08.2003	ggu	do not call setweg, setnod, setkan
c 23.03.2006	ggu	changed time step to real
c 28.09.2007	ggu	use testbndo to determine boundary node in flxtype
c 28.04.2009	ggu	links re-structured
c 23.03.2010	ggu	changed v6.1.1
c 23.02.2011	ggu	new routine call write_node_fluxes() for special output
c 01.03.2011	ggu	changed VERS_6_1_20
c 01.06.2011	ggu	documentation to flxscs() changed
c 14.07.2011	ggu	changed VERS_6_1_27
c 21.09.2011	ggu	some lower-level subroutines copied to subflx.f
c 07.10.2011	ggu	adjusted for 3d flux routines
c 18.10.2011	ggu	changed VERS_6_1_33
c 19.10.2011	ggu	added T/S variables, created fluxes_*() routines
c 19.10.2011	ggu	added conz variables, created fluxes_template()
c 09.12.2011	ggu	changed VERS_6_1_38
c 01.06.2012	ggu	changed VERS_6_1_53
c 10.05.2013	ggu	introduced subflxa.h, common routines to subflxu.f
c 25.10.2013	ggu	changed VERS_6_1_68
c 07.03.2014	ggu	changed VERS_6_1_72
c 26.11.2014	ggu	changed VERS_7_0_7
c 19.12.2014	ggu	changed VERS_7_0_10
c 19.01.2015	ggu	changed VERS_7_1_3
c 20.05.2015	ggu	modules introduced
c 10.07.2015	ggu	changed VERS_7_1_50
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 09.11.2015	ggu	changed VERS_7_3_13
c 12.04.2016	ggu	fluxes_template adjourned
c 15.04.2016	ggu	fluxes_template debugged and finished
c 09.09.2016	ggu	changed VERS_7_5_17
c 31.03.2017	ggu	changed VERS_7_5_24
c 22.09.2017	ccf	added total sediment concentration
c 26.10.2017	ggu	reads itable, chflx and section description
c 04.11.2017	ggu	changed VERS_7_5_34
c 14.11.2017	ggu	changed VERS_7_5_36
c 17.11.2017	ggu	sediment output adapted to new framework
c 17.11.2017	ggu	changed VERS_7_5_38
c 27.03.2018	ggu	new code for generic flux handling (fluxes_generic())
c 03.04.2018	ggu	changed VERS_7_5_43
c 06.07.2018	ggu	changed VERS_7_5_48
c 16.02.2019	ggu	changed VERS_7_5_60
c
c notes :
c
c These routines can also be used internally to compute the flux
c over various sections. The following calling sequence must be respected:
c
c call flx_init(kfluxm,kflux,nsect,iflux)		initializes iflux
c
c call flxscs(kfluxm,kflux,iflux,az,fluxes) computes fluxes 
c
c Initialization can be done anytime.
c
c******************************************************************
c******************************************************************
c******************************************************************

!==================================================================
        module flux
!==================================================================

        implicit none

	integer, save :: nl_flux = 0
	integer, save :: ns_flux = 0

        logical, save :: bflxinit = .false.
        integer, save :: nsect = -1
        integer, save :: kfluxm = 0
        integer, save, allocatable :: kflux(:)
        integer, save, allocatable :: iflux(:,:)
        integer, save, allocatable :: itable(:,:)
        character*80, save, allocatable :: chflx(:)

        integer, save, allocatable :: nlayers(:)
        real, save, allocatable :: fluxes(:,:,:)

        real, save, allocatable :: masst(:,:,:)
        real, save, allocatable :: saltt(:,:,:)
        real, save, allocatable :: tempt(:,:,:)
        real, save, allocatable :: conzt(:,:,:)
        real, save, allocatable :: ssctt(:,:,:)

!==================================================================
        contains
!==================================================================

!==================================================================
        end module flux
!==================================================================

        subroutine flux_read_section(n,ns)

        use flux
        use nls

        integer n,ns

	call nls_init_section

        !n = nls_read_vector()
	call nls_read_isctable(n,ns)
        kfluxm = n
	nsect = ns

        if( n > 0 ) then
          allocate(kflux(n))
          allocate(iflux(3,n))
          allocate(itable(2,ns))
          allocate(chflx(ns))
	  chflx = ' '
	  call nls_copy_isctable(n,ns,kflux,itable,chflx)
          !call nls_copy_int_vect(n,kflux)
        end if

	call nls_finish_section

        end subroutine flux_read_section

c******************************************************************

        subroutine flux_alloc_arrays(nl,ns)

	use flux

	implicit none

	integer nl	!layers
	integer ns	!sections

	if( nl == nl_flux .and. ns == ns_flux ) return

	if( nl > 0 .or. ns > 0 ) then
	  if( nl == 0 .or. ns == 0 ) then
            write(6,*) 'nl,ns: ',nl,ns
            stop 'error stop flux_alloc_arrays: incompatible parameters'
	  end if
	end if

	if( ns_flux > 0 ) then
          deallocate(nlayers)
          deallocate(fluxes)
          deallocate(masst)
          deallocate(saltt)
          deallocate(tempt)
          deallocate(conzt)
          deallocate(ssctt)
	end if

	nl_flux = nl
	ns_flux = ns

	if( ns == 0 ) return

        allocate(nlayers(ns))
        allocate(fluxes(0:nl,3,ns))

        allocate(masst(0:nl,3,ns))
        allocate(saltt(0:nl,3,ns))
        allocate(tempt(0:nl,3,ns))
        allocate(conzt(0:nl,3,ns))
        allocate(ssctt(0:nl,3,ns))

	end

c******************************************************************
c******************************************************************
c******************************************************************

        subroutine mod_flx(mode)
 
        implicit none
 
        integer mode
 
        include 'modules.h'
 
        if( mode .eq. M_AFTER ) then
           call wrflxa
        else if( mode .eq. M_INIT ) then
           call inflxa
        else if( mode .eq. M_READ ) then
           call rdflxa
        else if( mode .eq. M_CHECK ) then
           call ckflxa
        else if( mode .eq. M_SETUP ) then
           call flxini
        else if( mode .eq. M_PRINT ) then
           call prflxa
        else if( mode .eq. M_TEST ) then
           call tsflxa
        else if( mode .eq. M_BEFOR ) then
c          nothing
        else
           write(6,*) 'unknown mode : ', mode
           stop 'error stop mod_flx'
        end if
 
        end

c******************************************************************

        subroutine inflxa

c nsect		total number of sections
c kfluxm	total number of nodes defining sections
c kflux()	node numbers defining sections

        implicit none

        end

c******************************************************************

        subroutine rdflxa

	use flux

        implicit none

	integer n,ns

        call flux_read_section(n,ns)

        if( n .lt. 0 ) then
          write(6,*) 'read error in section $flux'
          stop 'error stop rdflxa'
        end if

        end

c******************************************************************

        subroutine ckflxa

	use flux

        implicit none

	integer k,ii
        logical berror

	berror = .false.
	if( kfluxm > 0 ) call n2int(kfluxm,kflux,berror)

        if( berror ) then
		write(6,*) 'error in section FLUX'
		stop 'error stop: ckflxa'
	end if

c initialize vectors (not strictly necessary)

	if( kfluxm > 0 ) iflux = 0

c the real set up is done in flxini
c but since at this stage we do not have all the arrays set up
c we post-pone it until later

        end

c******************************************************************

	subroutine flxini

c initializes flux routines finally (wrapper for flx_init)

	use flux

	implicit none

	if( kfluxm == 0 ) return

	call flx_init(kfluxm,kflux,nsect,iflux)
	bflxinit = .true.

	end

c******************************************************************

	subroutine prflxa

	use flux

	implicit none

	integer nnode,ifirst,ilast
	integer ntotal,ns
	integer i,ii

	integer ipext
	logical nextline

	write(6,*)
	write(6,*) 'flux section :'
	write(6,*)
	write(6,*) 'nsect,kfluxm ',nsect,kfluxm
	write(6,*)

	if( kfluxm == 0 ) return

	ns = 0
	nnode = 0

	do while( nextline(kflux,kfluxm,nnode,ifirst,ilast) )
	  ns = ns + 1
	  ntotal = ilast - ifirst + 1
	  write(6,*) 'section : ',ns,ntotal,'  ',trim(chflx(ns))
	  do i=ifirst,ilast
	    write(6,*) ipext(kflux(i)),(iflux(ii,i),ii=1,3)
	  end do
	end do

	end

c******************************************************************

	subroutine tsflxa

	use flux

	implicit none

	integer i,ii

	write(6,*) '/kfluxc/'
	write(6,*) nsect,kfluxm
	write(6,*) (kflux(i),i=1,kfluxm)

	write(6,*) '/iflux/'
	write(6,*) ((iflux(ii,i),ii=1,3),i=1,kfluxm)

	end

c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************

	subroutine wrflxa

c administers writing of flux data

	use mod_conz, only : cnv
	use mod_ts
	use mod_sediment, only : tcn
	use levels, only : nlvdi,nlv
	use flux

	implicit none

	include 'simul.h'

	integer j,i,l,lmax,nlmax,ivar,nvers
	integer idtflx,ierr,iv
	integer kext(kfluxm)
	real az,azpar,dt
	double precision atime0,atime
	character*80 title,femver

	integer ifemop,ipext
	real getpar
	double precision dgetpar
	logical has_output_d,next_output_d,is_over_output_d

        real, save :: trm,trs,trt,trc,trsc
        double precision, save :: da_out(4)
        integer, save :: nbflx = 0
	logical, save :: btemp,bsalt,bconz,bsedi
	integer, save :: nvar

c-----------------------------------------------------------------
c start of code
c-----------------------------------------------------------------

        if( nbflx .eq. -1 ) return

c-----------------------------------------------------------------
c initialization
c-----------------------------------------------------------------

        if( nbflx .eq. 0 ) then

		btemp = ( nint(getpar('itemp')) > 0 )
		bsalt = ( nint(getpar('isalt')) > 0 )
		bconz = ( nint(getpar('iconz')) == 1 )
		bsedi = ( nint(getpar('isedi')) > 0 )

		nvar = 1
		if( btemp ) nvar = nvar + 1
		if( bsalt ) nvar = nvar + 1
		if( bconz ) nvar = nvar + 1
		if( bsedi ) nvar = nvar + 1

		call init_output_d('itmflx','idtflx',da_out)
		call increase_output_d(da_out)
                if( .not. has_output_d(da_out) ) nbflx = -1

                if( kfluxm .le. 0 ) nbflx = -1
                if( nsect .le. 0 ) nbflx = -1
                if( nbflx .eq. -1 ) return

        	call flux_alloc_arrays(nlvdi,nsect)
		call get_nlayers(kfluxm,kflux,nlayers,nlmax)

		call fluxes_init(nlvdi,nsect,nlayers,trm,masst)
		if( bsalt ) then
		  call fluxes_init(nlvdi,nsect,nlayers,trs,saltt)
		end if
		if( btemp ) then
		  call fluxes_init(nlvdi,nsect,nlayers,trt,tempt)
		end if
		if( bconz ) then
		  call fluxes_init(nlvdi,nsect,nlayers,trc,conzt)
		end if
		if( bsedi ) then
		  call fluxes_init(nlvdi,nsect,nlayers,trsc,ssctt)
		end if

                nbflx=ifemop('.flx','unform','new')
                if(nbflx.le.0) goto 99
		da_out(4) = nbflx

	        nvers = 5
		idtflx = nint(da_out(1))
		call flx_write_header(nbflx,0,nsect,kfluxm,idtflx
     +                                  ,nlmax,nvar,ierr)
		if( ierr /= 0 ) goto 98

		title = descrp
                call get_shyfem_version_and_commit(femver)
                call get_absolute_ref_time(atime0)

		do i=1,kfluxm
		  kext(i) = ipext(kflux(i))
		end do

        	call flx_write_header2(nbflx,0,nsect,kfluxm
     +                          ,kext,nlayers
     +                          ,atime0,title,femver,chflx,ierr)
		if( ierr /= 0 ) goto 98

c               here we could also compute and write section in m**2

        end if

c-----------------------------------------------------------------
c normal call
c-----------------------------------------------------------------

        if( .not. is_over_output_d(da_out) ) return

	call get_timestep(dt)
	call getaz(azpar)
	az = azpar

c	-------------------------------------------------------
c	accumulate results
c	-------------------------------------------------------

	ivar = 0
	call flxscs(kfluxm,kflux,iflux,az,fluxes,ivar,rhov)
	call fluxes_accum(nlvdi,nsect,nlayers,dt,trm,masst,fluxes)

	if( bsalt ) then
	  ivar = 11
	  call flxscs(kfluxm,kflux,iflux,az,fluxes,ivar,saltv)
	  call fluxes_accum(nlvdi,nsect,nlayers,dt,trs,saltt,fluxes)
	end if
	if( btemp ) then
	  ivar = 12
	  call flxscs(kfluxm,kflux,iflux,az,fluxes,ivar,tempv)
	  call fluxes_accum(nlvdi,nsect,nlayers,dt,trt,tempt,fluxes)
	end if
	if( bconz ) then
	  ivar = 10
	  call flxscs(kfluxm,kflux,iflux,az,fluxes,ivar,cnv)
	  call fluxes_accum(nlvdi,nsect,nlayers,dt,trc,conzt,fluxes)
	end if
	if( bsedi ) then
	  ivar = 800
	  call flxscs(kfluxm,kflux,iflux,az,fluxes,ivar,tcn)
	  call fluxes_accum(nlvdi,nsect,nlayers,dt,trsc,ssctt,fluxes)
	end if

c	-------------------------------------------------------
c	time for output?
c	-------------------------------------------------------

        if( .not. next_output_d(da_out) ) return

        nbflx = nint(da_out(4))
        call get_absolute_act_time(atime)

c	-------------------------------------------------------
c	average and write results
c	-------------------------------------------------------

	ivar = 0
	iv = 1
	call fluxes_aver(nlvdi,nsect,nlayers,trm,masst,fluxes)
	call flx_write_record(nbflx,nvers,atime,nlvdi,nsect,ivar
     +				,nlayers,fluxes,ierr)
	if( ierr /= 0 ) goto 97

	if( bsalt ) then
	  ivar = 11
	  iv = iv + 1
	  call fluxes_aver(nlvdi,nsect,nlayers,trs,saltt,fluxes)
	  call flx_write_record(nbflx,nvers,atime,nlvdi,nsect,ivar
     +				,nlayers,fluxes,ierr)
	  if( ierr /= 0 ) goto 97
	end if
	if( btemp ) then
	  ivar = 12
	  iv = iv + 1
	  call fluxes_aver(nlvdi,nsect,nlayers,trt,tempt,fluxes)
	  call flx_write_record(nbflx,nvers,atime,nlvdi,nsect,ivar
     +				,nlayers,fluxes,ierr)
	  if( ierr /= 0 ) goto 97
	end if
	if( bconz ) then
	  ivar = 10
	  iv = iv + 1
	  call fluxes_aver(nlvdi,nsect,nlayers,trc,conzt,fluxes)
	  call flx_write_record(nbflx,nvers,atime,nlvdi,nsect,ivar
     +				,nlayers,fluxes,ierr)
	  if( ierr /= 0 ) goto 97
	end if
	if( bsedi ) then
	  ivar = 800
	  iv = iv + 1
	  call fluxes_aver(nlvdi,nsect,nlayers,trsc,ssctt,fluxes)
	  call flx_write_record(nbflx,nvers,atime,nlvdi,nsect,ivar
     +				,nlayers,fluxes,ierr)
	  if( ierr /= 0 ) goto 97
	end if

	if( iv /= nvar ) goto 91

c	-------------------------------------------------------
c	reset variables
c	-------------------------------------------------------

	call fluxes_init(nlvdi,nsect,nlayers,trm,masst)

	if( bsalt ) then
	  call fluxes_init(nlvdi,nsect,nlayers,trs,saltt)
	end if
	if( btemp ) then
	  call fluxes_init(nlvdi,nsect,nlayers,trt,tempt)
	end if

	if( bconz ) then
	  call fluxes_init(nlvdi,nsect,nlayers,trc,conzt)
	end if

	if( bsedi ) then
	  call fluxes_init(nlvdi,nsect,nlayers,trsc,ssctt)
	end if

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	return
   91   continue
        write(6,*) 'iv,nvar: ',iv,nvar
        write(6,*) 'iv is different from nvar'
        stop 'error stop wrflxa: internal error (1)'
   99   continue
        write(6,*) 'Error opening FLX file :'
        stop 'error stop wrflxa: opening flx file'
   98   continue
        write(6,*) 'Error writing header of FLX file'
        write(6,*) 'unit,ierr :',nbflx,ierr
        stop 'error stop wrflxa: writing flx header'
   97   continue
        write(6,*) 'Error writing file FLX'
        write(6,*) 'unit,ierr,iv :',nbflx,ierr,iv
        stop 'error stop wrflxa: writing flx record'

	end

!**********************************************************************
!**********************************************************************
!**********************************************************************
!**********************************************************************
!**********************************************************************

	subroutine fluxes_generic(ext,ivbase,nscal,scal)

! administers writing of flux data
!
! serves as a template for new variables
! please adapt to your needs
!
! this routine must be called after the other fluxes have already been set up
!
! here the number of scalars and the scalar values are passed into the routine
! you can also import them through other means (modules, etc..)
!
	use levels, only : nlvdi,nlv
	use basin, only : nkn
	use flux

	implicit none

	include 'simul.h'

	character*(*) ext		!extension of new file (e.g., '.fxw')
	integer ivbase			!base for variable numbers
	integer nscal			!how many tracers to compute/write
	real scal(nlvdi,nkn,nscal)	!tracer

	integer i,nlmax,ivar,nvers
	integer idtflx
	integer nvar,ierr
	integer kext(kfluxm)
	real az,dt
	double precision atime,atime0
	character*80 title,femver

        double precision, save :: da_out(4)
        integer, save :: nbflx = 0

	real, save, allocatable :: trs(:)
	real, save, allocatable :: scalt(:,:,:,:)	!accumulator array

	integer ifemop,ipext
	logical has_output_d,next_output_d,is_over_output_d
	double precision dgetpar

!-----------------------------------------------------------------
! start of code
!-----------------------------------------------------------------

        if( nbflx .eq. -1 ) return

!-----------------------------------------------------------------
! initialization
!-----------------------------------------------------------------

        if( nbflx .eq. 0 ) then

                call init_output_d('itmflx','idtflx',da_out)
                call increase_output_d(da_out)
                if( .not. has_output_d(da_out) ) nbflx = -1

                if( kfluxm .le. 0 ) nbflx = -1
                if( nsect .le. 0 ) nbflx = -1
                if( nbflx .eq. -1 ) return

		if( .not. bflxinit ) goto 94

        	allocate(trs(nscal))
        	allocate(scalt(0:nlvdi,3,nsect,nscal))

        	call flux_alloc_arrays(nlvdi,nsect)
		call get_nlayers(kfluxm,kflux,nlayers,nlmax)

		do i=1,nscal
		  call fluxes_init(nlvdi,nsect,nlayers,trs(i)
     +				,scalt(0,1,1,i))
		end do

                nbflx = ifemop(ext,'unform','new')
                if( nbflx .le. 0 ) goto 99
		write(6,*) 'flux file opened: ',nbflx,' ',ext
		da_out(4) = nbflx

	        nvers = 0
		nvar = nscal
		idtflx = nint(da_out(1))
                call flx_write_header(nbflx,0,nsect,kfluxm,idtflx
     +                                  ,nlmax,nvar,ierr)
                if( ierr /= 0 ) goto 98

                title = descrp
                call get_shyfem_version_and_commit(femver)
                call get_absolute_ref_time(atime0)

                do i=1,kfluxm
                  kext(i) = ipext(kflux(i))
                end do

                call flx_write_header2(nbflx,0,nsect,kfluxm
     +                          ,kext,nlayers
     +                          ,atime0,title,femver,chflx,ierr)
                if( ierr /= 0 ) goto 98

        end if

!-----------------------------------------------------------------
! normal call
!-----------------------------------------------------------------

        if( .not. is_over_output_d(da_out) ) return

	call get_timestep(dt)
	call getaz(az)

!	-------------------------------------------------------
!	accumulate results
!	-------------------------------------------------------

	do i=1,nscal
	  ivar = ivbase + i
	  call flxscs(kfluxm,kflux,iflux,az,fluxes,ivar,scal(1,1,i))
	  call fluxes_accum(nlvdi,nsect,nlayers,dt,trs(i)
     +			,scalt(0,1,1,i),fluxes)
	end do

!	-------------------------------------------------------
!	time for output?
!	-------------------------------------------------------

        if( .not. next_output_d(da_out) ) return

!	-------------------------------------------------------
!	average and write results
!	-------------------------------------------------------

        call get_absolute_act_time(atime)

	do i=1,nscal
	  ivar = ivbase + i
	  call fluxes_aver(nlvdi,nsect,nlayers,trs(i)
     +			,scalt(0,1,1,i),fluxes)
          call flx_write_record(nbflx,nvers,atime,nlvdi,nsect,ivar
     +                          ,nlayers,fluxes,ierr)
          if( ierr /= 0 ) goto 97
	end do

!	-------------------------------------------------------
!	reset variables
!	-------------------------------------------------------

	do i=1,nscal
	  call fluxes_init(nlvdi,nsect,nlayers,trs(i)
     +			,scalt(0,1,1,i))
	end do

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

	return
   94   continue
        write(6,*) 'Flux section has not been initialized'
        stop 'error stop fluxes_template: no initialization'
   97   continue
        write(6,*) 'Error writing data record of FLX file'
        write(6,*) 'unit,ierr :',nbflx,ierr
        stop 'error stop fluxes_template: writing flx record'
   98   continue
        write(6,*) 'Error writing headers of FLX file'
        write(6,*) 'unit,ierr :',nbflx,ierr
        stop 'error stop fluxes_template: writing flx header'
   99	continue
	write(6,*) 'extension: ',ext
        stop 'error stop fluxes_template: cannot open flx file'
	end

!******************************************************************
!******************************************************************
!******************************************************************
!******************************************************************
!******************************************************************

