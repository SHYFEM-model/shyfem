
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998-2000,2002-2003,2006-2007,2009  Georg Umgiesser
!    Copyright (C) 2011,2013-2015,2017-2020  Georg Umgiesser
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

! subroutines for computing discharge / flux for boxes
!
! contents :
!
! subroutine wrboxa(it)				write of flux data
!
! revision log :
!
! 30.04.1998	ggu	newly written routines (subpor deleted)
! 07.05.1998	ggu	check nrdveci on return for error
! 08.05.1998	ggu	restructured with new comodity routines
! 13.09.1999	ggu	type of node computed in own routine flxtype
! 19.11.1999	ggu	iskadj into sublin
! 20.01.2000	ggu	old routines substituted, new routine extrsect
! 20.01.2000	ggu	common block /dimdim/ eliminated
! 20.01.2000	ggu	common block /dimdim/ eliminated
! 01.09.2002	ggu	ggu99 -> bug in flx routines (how to reproduce?)
! 26.05.2003	ggu	in flxnov substituted a,b with b,c
! 26.05.2003	ggu	new routine make_fluxes (for lagrangian)
! 10.08.2003	ggu	do not call setweg, setnod, setkan
! 23.03.2006	ggu	changed time step to real
! 28.09.2007	ggu	use testbndo to determine boundary node in flxtype
! 28.04.2009	ggu	links re-structured
! 23.02.2011	ggu	new routine call write_node_fluxes() for special output
! 01.06.2011	ggu	documentation to flxscs() changed
! 21.09.2011	ggu	some lower-level subroutines copied to subflx.f
! 07.10.2011	ggu	adjusted for 3d flux routines
! 19.10.2011	ggu	added T/S variables, created fluxes_*() routines
! 19.10.2011	ggu	added conz variables, created fluxes_template()
! 10.05.2013	ggu	adapted to boxes computations
! 14.05.2013	ggu	write also OBC sections and contributions
! 13.06.2013	ggu	changed VERS_6_1_65
! 07.03.2014	ggu	changed VERS_6_1_72
! 29.10.2014	ggu	new code and 3d version
! 05.11.2014	ggu	changed VERS_7_0_5
! 26.11.2014	ggu	changed VERS_7_0_7
! 19.12.2014	ggu	changed VERS_7_0_10
! 19.01.2015	ggu	changed VERS_7_1_3
! 10.07.2015	ggu	changed VERS_7_1_50
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 24.07.2015	ggu	changed VERS_7_1_82
! 12.10.2015	ggu	changed VERS_7_3_3
! 09.11.2015	ggu	changed VERS_7_3_13
! 31.03.2017	ggu	changed VERS_7_5_24
! 13.04.2017	ggu	changed VERS_7_5_25
! 15.02.2018	ggu	code checked and error debug (is running)
! 22.02.2018	ggu	changed VERS_7_5_42
! 03.04.2018	ggu	changed VERS_7_5_43
! 06.07.2018	ggu	changed VERS_7_5_48
! 07.02.2019	ggu	code revised, easy addition of other vars
! 14.02.2019	ggu	changed VERS_7_5_56
! 16.02.2019	ggu	changed VERS_7_5_60
! 03.02.2020	ggu	revisted 3d box averaging
! 05.02.2020	ggu	bug in box_3d_aver_vertical() corrected
! 07.02.2020	ggu	final version of box file
! 20.03.2021	ggu	bug fix in box_write_stats()
! 26.06.2021	ggu	wrong units for rain and evaporation
! 27.06.2021	ggu	write header in all files
! 28.06.2021	ggu	flushing of output files
! 16.02.2022	ggu	write nvars in _geom file
! 09.03.2022	ggu	a comment that section nodes in index.txt are external
! 10.05.2022	ggu	reads new version of index file (boxes.txt)
! 30.05.2022	ggu	first changes for mpi use
! 01.06.2022	ggu	adjusted for mpi and 2d
! 08.06.2022	ggu	only master writes to output files
! 16.06.2022	ggu	only write header once (bug)
! 15.10.2022	ggu	shympi_exchange_array substituted with shympi_l2g_array
! 19.01.2023	ggu	adjourned for 3d array averaging and writing
! 29.01.2023	ggu	format of boxes_geom.txt has changed (with box nums)
! 01.04.2023	ggu	debugging code
! 21.04.2023	ggu	changes in writing OB fluxes and boxes_3d_aver()
! 18.05.2023	ggu	tons of changes to get 2D mass check right
! 29.05.2023	ggu	bug fix: all section arrays have same dimension (ns)
! 25.10.2023	ggu	eps2d and eps3d introduced
! 28.10.2023	ggu	bug fix for not existing boundaries
! 02.12.2023	ggu	changes marked with GGU_LAST
! 06.12.2023	ggu	more in assert_values(), 3d writing and check better
! 18.12.2023	ggu	new version 6
! 20.12.2023	ggu	bug fixes for flux computation, no OBC section
!
! notes :
!
! insert parameters idtbox and itmbox into STR file to have box file written
!
! for format of written files please see box_write_*()
! all other information is contained in boxfile (normally boxes.txt)
! boxfile is produced by shybas
!
! still to check: some sections have more layers than adjacent boxes
!
! in order to plot boxes set bbox=.true. in plobas (supsim.f)
!
! still to do
!
! open boundary fluxes in 3d, also compute z boundaries
! no molecular visc/fidd in last interface
! compute and write interface areas
! 
! attention: node numbers in boxes.txt file are external numbers
!            they are changed into internal numbers in box_init
!            just after reading the sections
!
! versions :
!
! 1		first version 
! 2		fairly standard old version
! 3		for MPI	(after 15.03.2022)
! 4		final version of with 3d computations under MPI
! 5		writing of boxes_geom.txt (with boxes and external nums)
!		boxes_geom.txt has different format
!		boxes_stats.txt has different section format
! 6		in 3d file also write act_eta and bstress
!		in sections block write section number
!
!******************************************************************
!
! calling sequence :
!
! init:
!
! 	bdtarea = 0
! 	call boxes_2d_init(nbox,nv2d,val2d)

! loop:
!
!	accumulation
!
!       call boxes_compute_area(barea)
!       bdtarea = bdtarea + dt*barea
!
!       call box_2d_aver_scalar(scalar,aux2d)	!aver over box and domain
!       call boxes_2d_accum(nbox,dt,barea,val2d(:,iv),aux2d)
!
!	writing
!
!       call boxes_2d_aver(nbox,nv2d,bdtarea,val2d)     !collects domains
!
!	reset
!
! 	bdtarea = 0
! 	call boxes_2d_init(nbox,nv2d,val2d)
!
! computation of fluxes:
!
!	call flxscs()			flux for all sections
!	    call flxsec()		flux for one section
!		call flx3d()		flux over one node in section
!			call flx3d_k()	compute divergence free flux transports
!
!	flxscs, flxsec in subflxu.f
!	flx3d, flx3d_k in subflx3d.f
!	flxtype, mkweig in subflx3d.f
!
!******************************************************************

!==================================================================
	module box
!==================================================================

	implicit none

        integer, parameter :: idbox = 473226		!id for box files
        integer, parameter :: nversbox = 6		!newest version

	logical, save :: bextra = .true.		!write extra info
	logical, save :: bflush = .true.		!flush after write
	character*80, save :: boxfile = 'boxes.txt'	!file name of box info
	logical, parameter :: bbox3d = .true.		!write 3d results
	logical, parameter :: bmasserror = .false.	!write mass error
	logical, parameter :: bmass2d = .true.		!check mass error 2d
	logical, parameter :: bmass3d = .true.		!check mass error 3d
	logical, parameter :: bassert = .true.		!check consistency
	logical, parameter :: b777 = .false.		!write debug to 777

        integer, save :: nbxdim = 0	!maximum for nbox
        integer, save :: nscboxdim = 0	!maximum for nsect
        integer, save :: nfxboxdim = 0	!maximum for kfluxm

	real, parameter :: eps2d = 1.e-3	!eps for 2d mass balance
	real, parameter :: eps3d = 1.e-3	!eps for 3d mass balance

        integer, save :: nbox		!total number of boxes
        integer, save :: nbc_ob		!total number of open boundaries
        integer, save :: nsect = -1	!total number of sections
        integer, save :: kfluxm = 0	!total number of nodes in all sections
        integer, save, allocatable :: kflux(:)	!node numbers defining sections
        integer, save, allocatable :: kflux_ext(:)	!external node numbers

        integer, save, allocatable :: iflux(:,:)
        integer, save, allocatable :: iboxes(:)
        integer, save, allocatable :: ikboxes(:)
        integer, save, allocatable :: isects(:,:)
        integer, save, allocatable :: iscbnd(:,:)

!       iflux(3,i)      aux values for nodes in sections
!       iboxes(ie)      index from elements to boxes
!       ikboxes(k)      index from nodes to boxes (<0 if node on box boundary)
!       isects(5,is)    description of section (n,type,box1,box2,ipnt)
!       iscbnd(5,ibc)   description of OB (n,type,box1,box2,is)
!
!       i runs on list of nodes (1:kfluxm)
!       ie runs on elements (1:nel)
!       k runs on nodes (1:nkn)
!       is runs on sections (1:nsect)
!       ibc runs on open boundaries (1:nbc)
!       ib runs on boxes (1:nbox)
!
! file types (ftype):
!
!	0	init
!	1	meteo
!	2	2d
!	3	3d
!	4	vertical
!	7	geom
!	8	stats
!	9	box index (boxes.txt)

!==================================================================
	contains
!==================================================================

	subroutine box_alloc(nkn,nel,nlv,nbc,nb,ns,nn)

	integer nkn,nel,nlv,nbc,nb,ns,nn

	allocate(kflux(nn))
	allocate(kflux_ext(nn))

	allocate(iflux(3,nn))
	allocate(iboxes(nel))
	allocate(ikboxes(nkn))
	allocate(isects(5,ns))
	allocate(iscbnd(5,nbc))

	end subroutine box_alloc

!==================================================================
	end module box
!==================================================================

!==================================================================
	module box_arrays
!==================================================================

	implicit none

	double precision, save, allocatable :: fluxes(:,:,:)	!aux fluxes
	double precision, save, allocatable :: aux2d(:)		!aux 2d vars
	double precision, save, allocatable :: aux3d(:,:)	!aux 3d vars

	real, save, allocatable :: wlaux(:,:)		!aux array
	real, save, allocatable :: wl0aux(:,:)		!aux array for vert w

        integer, save :: nslmax = 0
	integer, save, allocatable :: nslayers(:)	!layers in section
	double precision, save, allocatable :: masst(:,:,:)	!accum mass
	double precision, save, allocatable :: saltt(:,:,:)	!accum salt
	double precision, save, allocatable :: tempt(:,:,:)	!accum temp
	double precision, save, allocatable :: conzt(:,:,:)	!accum conz

	integer, save, allocatable :: nslayers_ob(:)	!layers in section
	double precision, save, allocatable :: fluxes_ob(:,:,:)	!discharges OBC
	double precision, save, allocatable :: masst_ob(:,:,:)	!OBC (dis-accum)

	integer, save, allocatable :: nblayers(:)	!layers in boxes
	double precision, save, allocatable :: barea(:)		!area of boxes
	double precision, save, allocatable :: bvol2d(:)	!vol2d of boxes
	double precision, save, allocatable :: bdepth(:)	!depth of boxes
	double precision, save, allocatable :: bvol3d(:,:)	!vol of boxes
	double precision, save, allocatable :: bdtarea(:)	!area*dt
	double precision, save, allocatable :: bdtvol(:,:)	!vol*dt

	integer, parameter :: nv3d = 4
	double precision, save, allocatable :: val3d(:,:,:)	!3d variables

	integer, parameter :: nvv3d = 4
	double precision, save, allocatable :: valv3d(:,:,:)	!3d vert. vars

	integer, parameter :: nv2d = 2
	double precision, save, allocatable :: val2d(:,:)	!2d variables

	integer, parameter :: nvmet = 8
	double precision, save, allocatable :: valmet(:,:)	!meteo variables

	double precision, save, allocatable :: eta_act(:)	!new water level
	double precision, save, allocatable :: eta_old(:)	!old water level

	double precision, save, allocatable :: fluxesg(:,:,:)	!aux fluxes
	double precision, save, allocatable :: fluxesg_m(:,:,:)	!mass fluxes
	double precision, save, allocatable :: fluxesg_ob(:,:,:)!ob fluxes
	double precision, save, allocatable :: aux3dg(:,:)	!aux 3d vars
	double precision, save, allocatable :: val3dg(:,:,:)	!3d variables
	double precision, save, allocatable :: valv3dg(:,:,:)	!3d variables

!==================================================================
	contains
!==================================================================

	subroutine box_arrays_alloc(nlv,nlvg,nbc,nb,ns,nn)

	integer nlv,nlvg,nbc,nb,ns,nn

	if( nbc > ns ) then
	  write(6,*) 'nbc,ns: ',nbc,ns
	  stop 'error stop box_arrays_alloc: nbc>ns'
	end if

	allocate(fluxes(0:nlv,3,ns))
	allocate(aux2d(nb))
	allocate(aux3d(0:nlv,nb))

	allocate(nslayers(ns))
	allocate(masst(0:nlv,3,ns))
	allocate(saltt(0:nlv,3,ns))
	allocate(tempt(0:nlv,3,ns))
	allocate(conzt(0:nlv,3,ns))

	allocate(nslayers_ob(ns))
	allocate(fluxes_ob(0:nlv,3,ns))
	allocate(masst_ob(0:nlv,3,ns))

	allocate(nblayers(nb))
	allocate(bvol2d(nb))			!static volume - no eta
	allocate(bdepth(nb))

	allocate(barea(nb))
	allocate(bdtarea(nb))			!area*dt

	allocate(bvol3d(0:nlv,nb))
	allocate(bdtvol(0:nlv,nb))		!vol*dt

	allocate(val3d(0:nlv,nb,nv3d))
	allocate(valv3d(0:nlv,nb,nvv3d))
	allocate(val2d(nb,nv2d))
	allocate(valmet(nb,nvmet))

	allocate(fluxesg(0:nlvg,3,ns))		!global 3d values
	allocate(fluxesg_m(0:nlvg,3,ns))	!global 3d mass fluxes
	allocate(fluxesg_ob(0:nlvg,3,ns))		!global 3d values
	allocate(aux3dg(0:nlvg,nb))
	allocate(val3dg(0:nlvg,nb,nv3d))
	allocate(valv3dg(0:nlvg,nb,nvv3d))

	allocate(eta_old(nb))
	allocate(eta_act(nb))

	end subroutine box_arrays_alloc

!==================================================================
	end module box_arrays
!==================================================================

!******************************************************************
!******************************************************************
!******************************************************************
! general box routine (init/accum/aver/write)
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine wrboxa

! administers writing of flux data

	use mod_conz, only : cnv
	use mod_ts
	use mod_diff_visc_fric
	use mod_hydro_vel
	use mod_hydro
	use basin
	use levels, only : nlvdi,nlv,ilhkv
	use box
	use box_arrays
	use mod_debug
	use shympi
	!use simul				!GGU_LAST

	implicit none

	include 'simul.h'			!GGU_LAST

	integer j,i,k,l,lmax,nlmax,ivar,nvers,nk_ob
	integer date,time
	integer idtbox
	integer nvar,ierr,iuaux,ib,iv,iu,is,ibtyp
	real az,azpar,dt,dt1
	double precision atime0,atime,dtime
	double precision daux
	character*80 title,femver,aline

        double precision, save :: trm,trs,trt,trc
	double precision, save :: trob
	double precision, save :: dtbox		!need for time averaging

        integer, save :: nbflx = 0		!unit number for output
        integer, save :: icall = 0
	integer, save :: ibarcl,iconz,ievap
	double precision, save :: da_out(4)

	integer nbnds,nkbnd
	integer ifemop
	integer ipext
	real getpar
	double precision dgetpar
	logical has_output_d,is_over_output_d,next_output_d
	logical is_first_output_d

	logical bfirst
	logical b3d

	real, allocatable, save :: taubot(:)
	integer, allocatable :: kext(:)
	character*80, allocatable :: chflx(:)

!-----------------------------------------------------------------
! start of code
!-----------------------------------------------------------------

        if( icall .eq. -1 ) return

	b3d = nlv > 1

!-----------------------------------------------------------------
! initialization
!-----------------------------------------------------------------

        if( icall .eq. 0 ) then

          	call init_output_d('itmbox','idtbox',da_out)
		call increase_output_d(da_out)  !itbox=itmbox+idtbox
          	if( .not. has_output_d(da_out) ) icall = -1

                if( icall .eq. -1 ) return

		nbc_ob = nbnds()
		nk_ob = nkbnd()
		call box_check(nbox,nsect,kfluxm) !gets info from boxfile

                if( kfluxm .le. 0 ) icall = -1
                if( nsect .le. 0 ) icall = -1
                if( icall .eq. -1 ) return
		icall = 1

		nbxdim = nbox
		nscboxdim = nsect + nbc_ob
		nfxboxdim = kfluxm + nbc_ob + nk_ob + 1
                call box_alloc(nkn,nel,nlvdi,nbc_ob,nbox                &
     &					,nscboxdim,nfxboxdim)
                call box_arrays_alloc(nlv,nlv_global,nbc_ob,nbox        &
     &					,nscboxdim,nfxboxdim)
		allocate(wlaux(nlv,nkn))
		allocate(wl0aux(0:nlv,nkn))

		!------------------------------------------------------
		! call to box_init() reads boxfile and sets up nslayers
		!------------------------------------------------------

		call box_init

		!------------------------------------------------------
		! call to box_make_stats() sets up nblayers
		!------------------------------------------------------

                call box_make_stats(nbox,iboxes,nblayers                &
     &					,barea,bvol2d,bdepth)

		call finalize_nslayers

                if( nsect .gt. nscboxdim ) then
                  stop 'error stop wrboxa: dimension nscboxdim'
                end if

		if( shympi_is_parallel() ) then
		  if( shympi_is_master() ) then
		    write(6,*) 'box model working in mpi mode'
		  end if
		  flush(6)
		  call shympi_syncronize
		  !stop 'error stop wrboxa: no mpi mode'
		end if

		ibarcl = nint(getpar('ibarcl'))
		iconz = nint(getpar('iconz'))
		ievap = nint(getpar('ievap'))

		call fluxes_init_d(nlvdi,nsect,nslayers,trm,masst)
		call fluxes_init_d(nlvdi,nbc_ob,nslayers_ob,trob,masst_ob)
		!call fluxes_init_d(nlvdi,nsect,nslayers,trob,masst_ob)
		nvar = 1
		if( ibarcl .gt. 0 ) then
		  call fluxes_init_d(nlvdi,nsect,nslayers,trs,saltt)
		  call fluxes_init_d(nlvdi,nsect,nslayers,trt,tempt)
		  nvar = nvar + 2
		end if
		if( iconz .eq. 1 ) then
		  call fluxes_init_d(nlvdi,nsect,nslayers,trc,conzt)
		  nvar = nvar + 1
		end if

		dtbox = 0.
		bdtarea = 0.
		bdtvol = 0.
		call boxes_3d_init(nlvdi,nbox,nv3d,val3d)   !3d variables
		call boxes_3d_init(nlvdi,nbox,nvv3d,valv3d) !3d vert variables
		call boxes_2d_init(nbox,nv2d,val2d)	!2d variables
		call boxes_2d_init(nbox,nvmet,valmet)	!meteo variables

		!write(6,*) 'gggguuuu 1',my_id; flush(6)
		call shympi_syncronize
		call box_aver_eta(zeov,aux2d)
		eta_old = aux2d

		!write(6,*) 'gggguuuu 2',my_id; flush(6)
		!call shympi_syncronize
	
        	call box_flx_file_open(nvar,da_out)
		!call shympi_syncronize
		!write(6,*) 'gggguuuu 3',my_id; flush(6)

!               here we could also compute and write section in m**2

                date = nint(dgetpar('date'))
                time = nint(dgetpar('time'))
		idtbox = nint(da_out(1))
                call box_write_stats(date,time,idtbox                   &
     &                                  ,nbox,nsect,iboxes,isects       &
     &                                  ,kfluxm,kflux,kflux_ext         &
     &                                  ,nslayers,nblayers              &
     &                                  ,barea,bvol2d,bdepth            &
     &					,bextra)

		allocate(taubot(nkn))

		!call write_section_info
		call write_matrix_boxes
        end if

	! at this point nsect is total number of sections (intern+OB)

!-----------------------------------------------------------------
! normal call
!-----------------------------------------------------------------

	if( is_first_output_d(da_out) ) then	!initial condition
	  call boxes_compute_area(barea)
	  call box_aver_eta(zenv,aux2d)
	  eta_old = aux2d * barea
	  call boxes_2d_aver(nbox,1,barea,eta_old)	!2d variables
	  call box_write_init(nbox,eta_old,bextra)
	end if

        if( .not. is_over_output_d(da_out) ) return	! it <= itmbox

	dt1 = 1.
	call get_timestep(dt)
	call getaz(azpar)
	az = azpar

!	-------------------------------------------------------
!	accumulate fluxes (is done on local domains)
!	-------------------------------------------------------

	ivar = 0
	call flxscs(kfluxm,kflux,iflux,az,fluxes,ivar,rhov)	!rhov is ignored
	call fluxes_accum_d(nlvdi,nsect,nslayers,dt,trm,masst,fluxes)

	call box_ob_compute(nbc_ob,nlvdi,fluxes_ob)  !open boundary conditions
        call fluxes_accum_d(nlvdi,nbc_ob,nslayers_ob,dt,trob            &
     &				,masst_ob,fluxes_ob)

	if( ibarcl .gt. 0 ) then
	  ivar = 11
	  call flxscs(kfluxm,kflux,iflux,az,fluxes,ivar,saltv)
	  call fluxes_accum_d(nlvdi,nsect,nslayers,dt,trs,saltt,fluxes)
	  ivar = 12
	  call flxscs(kfluxm,kflux,iflux,az,fluxes,ivar,tempv)
	  call fluxes_accum_d(nlvdi,nsect,nslayers,dt,trt,tempt,fluxes)
	end if

	if( iconz .eq. 1 ) then
	  ivar = 10
	  call flxscs(kfluxm,kflux,iflux,az,fluxes,ivar,cnv)
	  call fluxes_accum_d(nlvdi,nsect,nslayers,dt,trc,conzt,fluxes)
	end if

!	-------------------------------------------------------
!	accumulate box variables (is done on local domains)
!	-------------------------------------------------------

	dtbox = dtbox + dt

	call boxes_compute_area(barea)
	bdtarea = bdtarea + dt*barea
	call boxes_compute_volume(bvol3d)
	bdtvol = bdtvol + dt*bvol3d

	!call box_debug_2d('bvol3d',330,nbox,bvol3d(0,:))

	call box_3d_aver_scalar(tempv,aux3d,1)
	call boxes_3d_accum(nlvdi,nbox,dt,bvol3d,val3d(:,:,1),aux3d)
	call box_3d_aver_scalar(saltv,aux3d,1)
	call boxes_3d_accum(nlvdi,nbox,dt,bvol3d,val3d(:,:,2),aux3d)
	call box_3d_aver_vel(aux3d)
	call boxes_3d_accum(nlvdi,nbox,dt,bvol3d,val3d(:,:,3),aux3d)
	!val3d(:,:,4) = bvol3d(:,:)	!GGU_LAST - not needed

	wl0aux = wlnv
	where( wl0aux < 0. ) wl0aux = 0.		!positive vertical flux
	call box_3d_aver_vertical(wl0aux,aux3d,0)
	call boxes_3d_accum(nlvdi,nbox,dt,bvol3d,valv3d(:,:,1),aux3d)
	wl0aux = wlnv
	where( wl0aux > 0. ) wl0aux = 0.		!negative vertical flux
	call box_3d_aver_vertical(wl0aux,aux3d,0)
	call boxes_3d_accum(nlvdi,nbox,dt,bvol3d,valv3d(:,:,2),aux3d)

	call box_3d_aver_vertical(visv,aux3d,1)
	call boxes_3d_accum(nlvdi,nbox,dt,bvol3d,valv3d(:,:,3),aux3d)
	call box_3d_aver_vertical(difv,aux3d,1)
	call boxes_3d_accum(nlvdi,nbox,dt,bvol3d,valv3d(:,:,4),aux3d)

	call box_aver_eta(zenv,aux2d)
	call boxes_2d_accum(nbox,dt,barea,val2d(:,1),aux2d)
	eta_act = aux2d * barea

	call bottom_stress(taubot)
	call box_2d_aver_scalar(taubot,aux2d)
	call boxes_2d_accum(nbox,dt,barea,val2d(:,2),aux2d)

	call boxes_meteo_accum(nbox,dt,nvmet,barea,valmet,ievap)

	if( bassert ) call assert_values

!-----------------------------------------------------------------
! time for output?
!--------------------------------------------------------

        if( .not. next_output_d(da_out) ) return

	nbflx = nint(da_out(4))

!	-------------------------------------------------------
!	time average results and collect domains
!	-------------------------------------------------------

	nvers = 0
	call get_absolute_act_time(atime)
	call get_act_dtime(dtime)

	!------------------------------------------------------
	! average and gather eta and barea
	!------------------------------------------------------

	call box_aver_eta(zenv,aux2d)
	eta_act = aux2d * barea
	call boxes_2d_aver(nbox,1,barea,eta_act)	!2d variables

	call gather_sum_d2(nbox,barea)

	!------------------------------------------------------
	! average, gather, and write fluxes on section
	!------------------------------------------------------

	ivar = 0
	call fluxes_aver_d(nlvdi,nsect,nslayers,trm,masst,fluxes)
	call box_flx_collect(nlvdi,fluxes,fluxesg_m)

        call fluxes_aver_d(nlvdi,nbc_ob,nslayers_ob,trob                &
     &					,masst_ob,fluxes_ob)
	call box_flx_collect(nlvdi,fluxes_ob,fluxesg_ob)

	!write(6,*) 'transfering OB fluxes to section...'

	do i=1,nbc_ob
	  is = iscbnd(5,i)
	  ibtyp = iscbnd(2,i)
	  if( ibtyp == 3 ) fluxesg_m(:,:,is) = fluxesg_ob(:,:,i)
	end do

	call box_flx_write(atime,nbflx,ivar,fluxesg_m)

	if( ibarcl .gt. 0 ) then
	  ivar = 11
	  call fluxes_aver_d(nlvdi,nsect,nslayers,trs,saltt,fluxes)
	  call box_flx_collect(nlvdi,fluxes,fluxesg)
	  call box_flx_write(atime,nbflx,ivar,fluxesg)
	  ivar = 12
	  call fluxes_aver_d(nlvdi,nsect,nslayers,trt,tempt,fluxes)
	  call box_flx_collect(nlvdi,fluxes,fluxesg)
	  call box_flx_write(atime,nbflx,ivar,fluxesg)
	end if

	if( iconz .eq. 1 ) then
	  ivar = 10
	  call fluxes_aver_d(nlvdi,nsect,nslayers,trc,conzt,fluxes)
	  call box_flx_collect(nlvdi,fluxes,fluxesg)
	  call box_flx_write(atime,nbflx,ivar,fluxesg)
	end if

	!------------------------------------------------------
	! average, gather, and write box values
	!------------------------------------------------------

	call boxes_3d_aver(nlvdi,nbox,nv3d,bdtvol,val3d,val3dg)
	call boxes_3d_aver(nlvdi,nbox,nvv3d,bdtvol,valv3d,valv3dg)
	!call boxes_2d_aver(nbox,1,barea,eta_act)	   !eta
	call boxes_2d_aver(nbox,nv2d,bdtarea,val2d)	   !2d variables
	call boxes_2d_aver(nbox,nvmet,bdtarea,valmet)	   !meteo variables

	call boxes_3d_sum(nlvdi,nbox,bdtvol,aux3dg)
	val3dg(:,:,4) = aux3dg / dtbox

	!------------------------------------------------------
	! mass conservation check
	!------------------------------------------------------

	if( bmass2d ) then
	  aux2d = barea*(eta_act-eta_old)
	  call print_2d(nbox,barea,bvol2d,eta_act,eta_old)
	  call print_3d_flx(nsect,nbc_ob,dtbox,fluxesg_m,fluxesg_ob)
          call boxes_mass_balance_2d(dtbox,bvol2d                       &
     &					,fluxesg_m,fluxesg_ob,aux2d)
	end if
	if( b3d .and. bmass3d ) then
	  aux3dg = valv3dg(:,:,1) + valv3dg(:,:,2)	   !vertical velocity
          call boxes_mass_balance_3d(dtbox,bvol2d                       &
     &                                  ,nblayers,nslayers              &
     &                                  ,fluxesg_m,fluxesg_ob           &
     &                                  ,aux2d,aux3dg)
	end if

!	-------------------------------------------------------
!	write results
!	-------------------------------------------------------

	call get_act_timeline(aline)
	if( shympi_is_master() ) then
	  write(6,*) 'writing box results... ',trim(aline)
	end if

        call box_write_2d(dtime,aline,fluxesg_m,fluxesg_ob              &
     &                                  ,nv3d,val3dg                    &
     &                                  ,nv2d,val2d                     &
     &					,eta_act)
	call box_write_meteo(dtime,aline,nvmet,valmet)

	if( bbox3d .and. b3d ) then
          call box_write_3d(dtime,aline,nblayers,nslayers               &
     &                  ,fluxesg_m,fluxesg_ob                           &
     &                  ,nv3d,val3dg                                    &
     &                  ,nv2d,val2d                                     &
     &			,eta_act)
          call box_write_vertical(dtime,aline,nblayers                  &
     &			,nvv3d,valv3d)
	end if

	!stop 'forced stop after write'

!	-------------------------------------------------------
!	reset variables
!	-------------------------------------------------------

	call fluxes_init_d(nlvdi,nsect,nslayers,trm,masst)
	call fluxes_init_d(nlvdi,nbc_ob,nslayers_ob,trob,masst_ob)

	if( ibarcl .gt. 0 ) then
	  call fluxes_init_d(nlvdi,nsect,nslayers,trs,saltt)
	  call fluxes_init_d(nlvdi,nsect,nslayers,trt,tempt)
	end if

	if( iconz .eq. 1 ) then
	  call fluxes_init_d(nlvdi,nsect,nslayers,trc,conzt)
	end if

	dtbox = 0.
	bdtarea = 0.
	bdtvol = 0.
	call boxes_3d_init(nlvdi,nbox,nv3d,val3d)	!3d variables
	call boxes_3d_init(nlvdi,nbox,nvv3d,valv3d)	!3d vert variables
	call boxes_2d_init(nbox,nv2d,val2d)		!2d variables
	call boxes_2d_init(nbox,nvmet,valmet)		!meteo variables

	eta_old = eta_act

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

	return
   97	continue
	write(6,*) 'error writing record of boxes flx file: ',ivar
	stop 'error stop wrboxa: record'
   98	continue
	write(6,*) 'error writing header of boxes flx file'
	stop 'error stop wrboxa: header'
	end

!******************************************************************
!******************************************************************
!******************************************************************
! reading and setup routines
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine box_init

! reads file etc...

	use box
	use box_arrays
	use shympi

	implicit none

	logical berror
	integer nbc
	integer nbnds
	integer nsaux
	integer i,n

	call box_read				!reads boxfile
	call box_elab(iboxes,ikboxes)		!sets up ikboxes

	call convert_nodes(kfluxm,kflux)	!converts ext to int nodes

	call box_ob_init	!initializes OB contributions (internal nodes)

	call shympi_barrier

	! in the next call nslayers are set up

	nslayers = 0
	call flux_initialize(kfluxm,kflux,iflux,nsect,nslayers,nslmax)

	return

!	from here debug...

	write(6,*) nsect
	write(6,*) kfluxm
	write(6,*) kflux(1:kfluxm)
	do i=1,kfluxm
	  write(6,*) iflux(:,i)
	end do
	stop

	end

!******************************************************************

	subroutine box_check(nbx,nsct,kflxm)

! reads boxfile just to return dimensions

	use box

	implicit none

	integer nbx,nsct,kflxm

	integer nb,id,n,ib1,ib2,i
	integer nel,ie
	integer ftype,nvers
	integer iaux

	nbx = 0
	nsct = 0
	kflxm = 0

	open(1,file=boxfile,form='formatted',status='old',err=94)

	read(1,*) id,ftype,nvers

	if( id /= idbox .or. ftype /= 9 ) then	!old read: nel,nbox,nb
	  nel = id
	  nbx = ftype
	  read(1,*) (iaux,ie=1,nel)
	else					!new read
	  read(1,*) nel,nbx,nb
	  do ie=1,nel
	    read(1,*) iaux,iaux
	  end do
	end if

	do
	  read(1,*) id,n,ib1,ib2
	  if( id .le. 0 ) exit
	  read(1,*) (iaux,i=1,n)
	  nsct = nsct + 1
	  if( id .ne. nsct ) goto 97
	  kflxm = kflxm + n + 1
	end do

	close(1)

	return
   94	continue
	write(6,*) 'file: ',trim(boxfile)
	stop 'error stop box_check: cannot read file'
   97	continue
	write(6,*) id,nsect
	stop 'error stop box_check: id'
	end

!******************************************************************

	subroutine box_read

! reads boxfile and returns parameters and arrays
!
! on return the following values are set:
!
! nbox,nsect,kfluxm
! kflux

	use basin
	use box
	use shympi

	implicit none

	integer nb,id,n,ib1,ib2,i
	integer nelbox,ie
	integer ie_int,ie_ext,ia,ic
	integer ftype,nvers
	integer ierr,iaux,itype
	integer iudeb
	integer iauxv(nkn)

	integer ieint

	nsect = 0
	kfluxm = 0
	kflux = 0
	kflux_ext = 0
	ierr = 0
	iudeb = 666
	iudeb = 0

	write(6,*) 'start reading boxfile... ',trim(boxfile)
	if( iudeb > 0 ) write(iudeb,*) 'reading boxfile'
	open(1,file=boxfile,form='formatted',status='old',err=94)

	read(1,*) id,ftype,nvers

	if( id /= idbox .or. ftype /= 9 ) then	!old read: nel,nbox,nb
	  nelbox = id
	  nbox = ftype
	  nb = nvers
	  if( nelbox .ne. nel ) goto 99
	  read(1,*) (iboxes(ie),ie=1,nel)
	else					!new read: id,ftype,nvers
	  read(1,*) nelbox,nbox,nb
	  if( nelbox .ne. nel_global ) goto 99
	  iboxes = -1
	  do ie=1,nelbox
	    read(1,*) ie_ext,ia
	    ie_int = ieint(ie_ext)
	    if( ie_int > 0 ) iboxes(ie_int) = ia
	  end do
	  ic = count( iboxes == -1 )
	  if( ic > 0 ) then
	    write(6,*) 'some elements are not part of a box: ',ic
	    stop 'error stop box_read: incomplete box information'
	  end if
	end if
	if( nbox .gt. nbxdim ) goto 98

	itype = 99		!internal sections

	do
	  read(1,*) id,n,ib1,ib2
	  if( id .le. 0 ) exit
	  read(1,*) (iauxv(i),i=1,n)
	  nsect = nsect + 1
	  if( id .ne. nsect ) goto 97
	  if( nsect .gt. nscboxdim ) goto 95
	  call box_insert_section(id,n,itype,ib1,ib2,iauxv,ierr)
	  if( kfluxm .gt. nfxboxdim ) goto 96	!kfluxm set in sub above
	  if( ierr .gt. 0 ) goto 93		!should be already handled above
	  if( iudeb > 0 ) then
	    write(iudeb,*) id,n,ib1,ib2
	    write(iudeb,*) iauxv(1:n)
	  end if
	end do

	close(1)

	!kfluxm = kfluxm - 1		!we keep last 0
	kflux_ext = kflux

	write(6,*) 'finished reading boxfile: ',trim(boxfile)
	write(6,*) '  boxes    = ',nbox
	write(6,*) '  sections = ',nsect
	write(6,*) '  nodes    = ',kfluxm
	write(6,*) '  sdim     = ',nscboxdim
	write(6,*) '  ndim     = ',nfxboxdim
	!write(6,*) kflux(1:kfluxm)

	if( iudeb == 0 ) return

	write(iudeb,*) 'finished reading boxfile: ',trim(boxfile)
	write(iudeb,*) nbox,nsect,kfluxm,nscboxdim,nfxboxdim
	write(iudeb,*) nel
	do ie=1,nel
	  write(iudeb,*) ie,iboxes(ie)
	end do
	write(iudeb,*) 'finished writing file'
	close(iudeb)

	return
   93	continue
	write(6,*) ierr
	stop 'error stop box_read: unspecified error'
   94	continue
	write(6,*) 'file: ',trim(boxfile)
	stop 'error stop box_read: cannot read file'
   95	continue
	write(6,*) nsect,nscboxdim
	stop 'error stop box_read: nscboxdim'
   96	continue
	write(6,*) kfluxm,nfxboxdim
	stop 'error stop box_read: nfxboxdim'
   97	continue
	write(6,*) id,nsect
	stop 'error stop box_read: id'
   98	continue
	write(6,*) nbox,nbxdim
	stop 'error stop box_read: nbox>nbxdim'
   99	continue
	write(6,*) nel,nelbox
	stop 'error stop box_read: nel/=nelbox'
	end

!******************************************************************

	subroutine box_insert_section(id,n,itype,ib1,ib2,list,ierr)

! inserts section in data structure

	use basin
	use box

	implicit none

	integer id,n,itype
	integer ib1,ib2
	integer list(n)
	integer ierr

	integer i

	if( id .gt. nscboxdim ) ierr = ierr + 1
	if( kfluxm + n + 1 .gt. nfxboxdim ) ierr = ierr + 1

	if( ierr .eq. 0 ) then		!no errors -> insert
	    isects(1,id) = n		!total number of nodes
	    isects(2,id) = itype	!type of section (99 for internal)
	    isects(3,id) = ib1		!from box
	    isects(4,id) = ib2		!to box
	    isects(5,id) = kfluxm+1	!pointer into kflux
	    do i=1,n
	      kflux(kfluxm+i) = list(i)
	    end do
	    kfluxm = kfluxm + n + 1
	    kflux(kfluxm) = 0
	else
	    kfluxm = kfluxm + n + 1
	end if

	end

!******************************************************************

	subroutine box_elab(iboxes,ikboxes)

! sets up ikboxes - index of nodes to boxes
!
! is -1 if node belongs to more than one box

	use basin
	!use box

	implicit none

	integer iboxes(nel)
	integer ikboxes(nkn)

	integer k,ie,ii,ib
	integer nb,nl

	ikboxes = 0

	do ie=1,nel
	  ib = iboxes(ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( ikboxes(k) .eq. 0 ) then
	      ikboxes(k) = ib
	    else if( ikboxes(k) .ne. ib ) then
	      ikboxes(k) = -1
	    end if
	  end do
	end do

	nb = 0
	nl = 0
	do k=1,nkn
	  if( ikboxes(k) .gt. 0 ) then
	    nb = nb + 1
	  else if( ikboxes(k) .lt. 0 ) then
	    nl = nl + 1
	  end if
	end do

	if( nb+nl .ne. nkn ) then
	  write(6,*) nb,nl,nb+nl,nkn
	  stop 'error stop box_elab: internal error'
	end if

	write(6,*) 'box nodes: ',nb,nl,nb+nl,nkn

	end

!******************************************************************

        subroutine box_make_stats(nbox,iboxes,nblayers                  &
     &					,barea,bvol2d,bdepth)

! computes max layers for each box

	use mod_depth
	use evgeom
	use basin
	use levels
	use shympi

	implicit none

	integer nbox
	integer iboxes(nel)			!indicator of boxes
	integer nblayers(nbox)			!number of layers in boxes
	double precision barea(nbox)		!area of boxes
	double precision bvol2d(nbox)		!volume of boxes (static)
	double precision bdepth(nbox)		!depth of boxes

	integer ib,lmax,ie,nel_mpi
	real area,hdep

	nblayers = 0
	barea = 0.
	bvol2d = 0.
	bdepth = 0.
	nel_mpi = nel_unique

	do ie=1,nel_mpi
	  ib = iboxes(ie)
	  hdep = hev(ie)
	  area = 12.*ev(10,ie)
	  barea(ib) = barea(ib) + area
	  bvol2d(ib) = bvol2d(ib) + area * hdep
	  lmax = ilhv(ie)
	  nblayers(ib) = max(nblayers(ib),lmax)
	end do

	call shympi_barrier

	do ib=1,nbox
	  nblayers(ib) = shympi_max(nblayers(ib))
	end do

	call shympi_gather_and_sum(barea)
	call shympi_gather_and_sum(bvol2d)
	bdepth = bvol2d / barea

	!bdepth = 0.
	!where( barea > 0. ) bdepth = bvol2d / barea

	end

!******************************************************************
!******************************************************************
!******************************************************************
! writing routines
!******************************************************************
!******************************************************************
!******************************************************************

        subroutine box_write_stats(date,time,idtbox                     &
     &                                  ,nbox,nsect,iboxes,isects       &
     &                                  ,kfluxm,kflux,kflux_ext         &
     &                                  ,nslayers,nblayers              &
     &                                  ,barea,bvol2d,bdepth            &
     &					,bextra)

! writes statistics to files boxes_stats.txt and boxes_geom.txt

	use basin
	use mod_depth
	use evgeom
	use shympi

	implicit none

	integer date,time
	integer idtbox
	integer nbox,nsect
	integer iboxes(nel)
	integer isects(5,nsect)
	integer kfluxm,kflux(kfluxm),kflux_ext(kfluxm)
	integer nslayers(nsect)		!number of layers in section
	integer nblayers(nbox)		!number of layers in boxes
	double precision barea(nbox)		!area of boxes
	double precision bvol2d(nbox)		!volume of boxes
	double precision bdepth(nbox)		!depth of boxes
	logical bextra

	logical bw
	integer i,ipt,k,ke,n,ndim,ip
	integer ib,iu,iuaux
	integer is,ib1,ib2,itype,ns
	integer nb1,nb2
	integer ie,iext,ia
	integer nvars
	real area,depth
	double precision areatot
	double precision dtime,atime,atime0
	character*20 line
	character*80 string,s1,s2,s3,s4

	integer nsbox(-1:nbox)			!number of nodes in box section
	integer, allocatable :: kbox(:,:)	!nodes of section of box
	integer, allocatable :: kaux(:)
	integer, allocatable :: iboxesg(:)
	integer, allocatable :: ipevg(:)
	real, allocatable :: areav(:)
	real, allocatable :: areag(:)
	real, allocatable :: heg(:)
	real, allocatable :: xxx(:),yyy(:)

	real, parameter :: high = 1.e+30

	integer ifileo,ipext,ipint,ieext
	character*80 file

!-----------------------------------------------------------------
! write boxes_stats.txt file
!-----------------------------------------------------------------

	file = 'boxes_stats.txt'
	if( shympi_is_master() ) then
	  iu = ifileo(0,file,'formatted','new')
	  if( iu <= 0 ) stop 'error stop boxes: opening file'
	  bw = .true.
	else
	  iu = 987
	  bw = .false.
	  bextra = .false.
	end if
	if( bextra ) call box_write_header(iu,8)

	if( bextra ) write(iu,'(a)') '# information on simulation'
	if( bw ) write(iu,*) date,time
	if( bw ) write(iu,*) idtbox
	call get_absolute_ref_time(atime0)
	call get_first_dtime(dtime)
	atime = atime0 + dtime
        call dts_format_abs_time(atime,line)
	if( bw ) write(iu,*) 'start of simulation = ',line
	call get_last_dtime(dtime)
	atime = atime0 + dtime
        call dts_format_abs_time(atime,line)
	if( bw ) write(iu,*) '  end of simulation = ',line

	if( bextra ) write(iu,'(a)') '#      boxes       nvars'
	if( bw ) write(iu,*) nbox,4
	s1 = '#    box  layers'
	!     12345678901234567890123456789012345678901234567890
	s2 = '          area        volume       depth'
	if( bextra ) write(iu,'(a)') trim(s1)//trim(s2)
	do ib=1,nbox
          if( bw ) write(iu,1000) ib,nblayers(ib)                       &
     &				,barea(ib),bvol2d(ib),bdepth(ib)
 1000	  format(2i8,2e14.6,f12.4)
	end do

	nsbox = 0

	if( bextra ) write(iu,'(a)') '#   sections'
        if( bw ) write(iu,*) nsect
	s1 = '#    section        from          to'
	s2 = '         nls        nlb1        nlb2'
	if( bextra ) write(iu,'(a)') trim(s1)//trim(s2)
        do is=1,nsect
          ib1 = isects(3,is)
	  nb1 = 0
	  if( ib1 > 0 ) nb1 = nblayers(ib1)
          ib2 = isects(4,is)
	  nb2 = 0
	  if( ib2 > 0 ) nb2 = nblayers(ib2)
          if( bw ) write(iu,*) is,ib1,ib2,nslayers(is),nb1,nb2
          n = isects(1,is)
	  if( ib1 > 0 ) nsbox(ib1) = nsbox(ib1) + n + 1
	  if( ib2 > 0 ) nsbox(ib2) = nsbox(ib2) + n + 1
        end do

	ndim = maxval(nsbox)
	allocate(kbox(0:ndim,nbox))
	kbox = 0
	allocate(kaux(kfluxm),xxx(kfluxm),yyy(kfluxm))
	kaux(1:kfluxm) = kflux(1:kfluxm)
	call make_kexy(kfluxm,kaux,xxx,yyy,ns)

	iuaux = 440 + my_id
	if( bextra ) write(iu,'(a)') '#  section details'
        if( bw ) write(iu,*) nsect
	s1 = '#    section       nodes'
	s2 = '    from-box      to-box      layers'
	s3 = '#       node node-number'
	s4 = '            x                y'
        do is=1,nsect
	  if( bextra ) write(iu,'(a)') trim(s1)//trim(s2)
          n = isects(1,is)
          itype = isects(2,is)
          ib1 = isects(3,is)
          ib2 = isects(4,is)
          ipt = isects(5,is)
          !write(iuaux,*) is,n,ib1,ib2,nslayers(is)
          if( bw ) write(iu,*) is,n,ib1,ib2,nslayers(is)
	  if( bextra ) write(iu,'(a)') trim(s3)//trim(s4)
	  do i=1,n
	    ip = ipt-1+i
	    k = kflux(ip)
	    ke = kaux(ip)
	    k = ipint(ke)
	    if( bw ) write(iu,*) i,ke,xxx(ip),yyy(ip)
	    call insert_sect_node(ndim,nbox,ib1,kbox,k)
	    call insert_sect_node(ndim,nbox,ib2,kbox,k)
	  end do
	  call insert_sect_node(ndim,nbox,ib1,kbox,0)
	  call insert_sect_node(ndim,nbox,ib2,kbox,0)
	  !write(6,*) 'sssss',my_id,is,n
	  !write(iuaux,*) 'sssss',is,n
	end do

	if( bextra ) write(iu,'(a)') '#  box detail on sections'
	if( bw ) write(iu,*) nbox
	s1 = '#        box       nodes    sections'
	do ib=1,nbox
	  n = kbox(0,ib)
	  kaux(1:n) = kbox(1:n,ib)
	  call make_kexy(n,kaux,xxx,yyy,ns)
	  if( bextra ) write(iu,'(a)') trim(s1)
	  if( bw ) write(iu,*) ib,n,ns
	  if( bextra ) write(iu,'(a)') trim(s3)//trim(s4)
	  do i=1,n
	    k = kaux(i)
	    if( k == 0 ) then
	      if( bw ) write(iu,*) i,0,0.,0.
	    else
	      if( bw ) write(iu,*) i,k,xxx(i),yyy(i)
	    end if
	  end do
	end do

	if( bw ) close(iu)

!-----------------------------------------------------------------
! write boxes_geom.txt file
!-----------------------------------------------------------------

	file = 'boxes_geom.txt'
	if( shympi_is_master() ) then
	  iu = ifileo(0,file,'formatted','new')
	  if( iu <= 0 ) stop 'error stop boxes: opening file'
	  bw = .true.
	else
	  iu = 987
	  bw = .false.
	  bextra = .false.
	end if
	if( bextra ) call box_write_header(iu,7)

	allocate(areav(nel))
	allocate(areag(nel_global))
	allocate(heg(nel_global))
	allocate(iboxesg(nel_global))
	allocate(ipevg(nel_global))

	areatot = 0.
	!nvars = 2						!version 3
	nvars = 3
	if( bextra ) write(iu,'(a)') '#   elements       nvars'
	if( bw ) write(iu,*) nel_global,nvars
	!         12345678901234567890123456789012345678901234567890
	string = '#  element       box            area       depth'
	!string = '#  element            area       depth'	#version 3
	if( bextra ) write(iu,'(a)') trim(string)

	do ie=1,nel
	  area = 12. * ev(10,ie)
	  areav(ie) = area
	end do

	!call shympi_exchange_array(areav,areag)
	!call shympi_exchange_array(hev,heg)
	!call shympi_exchange_array(iboxes,iboxesg)
	call shympi_l2g_array(areav,areag)
	call shympi_l2g_array(hev,heg)
	call shympi_l2g_array(iboxes,iboxesg)
	call shympi_l2g_array(ipev,ipevg)

	do ie=1,nel_global
	  iext = ipevg(ie)
	  ia = iboxesg(ie)
	  area = areag(ie)
	  depth = heg(ie)
	  !if( bw ) write(iu,2000) ie,area,depth		!version 3
	  if( bw ) write(iu,2100) iext,ia,area,depth
 2000	  format(i10,e16.8,f12.4)
 2100	  format(2i10,e16.8,f12.4)
	end do
	write(6,*) 'total area: ',areatot

	if( bw ) close(iu)

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

	end

!******************************************************************

	subroutine insert_sect_node(ndim,nbox,ib,kbox,k)

	implicit none

	integer ndim,nbox,ib,k
	integer kbox(0:ndim,nbox)

	integer i

	if( ib <= 0 ) return

	i = kbox(0,ib)
	i = i + 1
	if( i > ndim ) then
	  write(6,*) ib,k,i,ndim
	  stop 'error stop insert_sect_node: i>ndim'
	end if
	kbox(i,ib) = k
	kbox(0,ib) = i

	end

!******************************************************************

	subroutine make_kexy(n,knodes,x,y,ns)

! computes external nodes and x/y

	use basin
	use shympi

	implicit none

	integer n
	integer knodes(n)
	real x(n),y(n)
	integer ns		!number of sections (return)

	real, parameter :: high = 1.e+30
	integer i,k

	integer ipext

	x = -high
	y = -high

	do i=1,n
	  k = knodes(i)
	  if( k > 0 ) then
	    x(i) = xgv(k)
	    y(i) = ygv(k)
	  end if
	  knodes(i) = ipext(k) 
	end do

	!write(6,*) 'make_kexy',my_id,n

	call shympi_barrier
	call gather_max_i(n,knodes)
	call gather_max_r(n,x)
	call gather_max_r(n,y)
	call shympi_barrier

	ns = 0
	do i=1,n
	  if( knodes(i) == 0 ) ns = ns + 1
	end do

	end

!******************************************************************

	subroutine box_write_init(nbox,eta_act,bextra)

! writes initial conditions boxes_init.txt for eta - still to be done : T/S

	use shympi

	implicit none

	integer nbox
	double precision eta_act(nbox)		!water level of first time step
	logical bextra,bw

	integer ib,iu
	character*80 file

	integer ifileo

	file = 'boxes_init.txt'
	if( shympi_is_master() ) then
	  iu = ifileo(0,file,'formatted','new')
	  if( iu <= 0 ) stop 'error stop boxes: opening file'
	  bw = .true.
	else
	  iu = 987
	  bw = .false.
	  bextra = .false.
	end if
	if( bextra ) call box_write_header(iu,0)

	if( bextra ) write(iu,'(a)') '#      boxes       nvars'
	if( bw ) write(iu,*) nbox,1
	if( bextra ) write(iu,'(a)') '#        box     init_eta'
	do ib=1,nbox
	  !write(iu,*) ib,eta_act(ib)
	  if( bw ) write(iu,'(i12,f13.3)') ib,eta_act(ib)
	end do

	if( bw ) close(iu)

	end

!******************************************************************

        subroutine box_write_2d(dtime,aline,fluxesg,fluxesg_ob          &
     &                          ,nv3d,val3dg                            &
     &                          ,nv2d,val2d                             &
     &				,eta_act)

! writes 2d vertical average box values to file boxes_2d.txt
!
! for boxes variables are:
!	1	temperature
!	2	salinity
!	3	water level
!	4	current velocity

	use basin
	use levels
	use box
	use shympi

	implicit none

	double precision dtime
	character*(*) aline
	double precision fluxesg(0:nlv_global,3,nsect)	!mass fluxes box to box
	double precision fluxesg_ob(0:nlv_global,3,nbc_ob) !mass fluxes at OB
	double precision flux2d(3,nsect)		!mass fluxes box to box
	integer nv3d
	double precision val3dg(0:nlv_global,nbox,nv3d)	!3d variables
	integer nv2d
	double precision val2d(nbox,nv2d)		!2d variables
	double precision eta_act(nbox)			!actual eta

	logical bw
	integer ib,is,ib1,ib2,ii,itype,iv,ns,ic
	integer itypes(nsect)

	integer ifileo
	character*80 file,s1,s2
	integer, save :: iu = 0

	file = 'boxes_2d.txt'
	if( shympi_is_master() ) then
	  if( iu == 0 ) then
	    iu = ifileo(0,file,'formatted','new')
	    if( iu <= 0 ) stop 'error stop boxes: opening file'
	    if( bextra ) call box_write_header(iu,2)
	  end if
	  bw = .true.
	else
	  iu = 987
	  bw = .false.
	  bextra = .false.
	  bflush = .false.
	end if

	if( bextra ) write(iu,'(a)') '# time record'
	if( bw ) write(iu,*) dtime,trim(aline)

	if( bextra ) write(iu,'(a)') '#      boxes       nvars'
	if( bw ) write(iu,*) nbox,nv3d+nv2d+1
	s1 = ' box      temp      salt cur_speed'
	s2 = '   act_eta  aver_eta   bstress        volume'
	!     12345678901234567890123456789012345678901234567890
	if( bextra ) write(iu,'(a)') '#'//trim(s1)//trim(s2)
	do ib=1,nbox
          if( bw ) write(iu,1000) ib                                    &
     &                          ,(val3dg(0,ib,iv),iv=1,nv3d-1)          &
     &                          ,eta_act(ib)                            &
     &                          ,(val2d(ib,iv),iv=1,nv2d)               &
     &				,val3dg(0,ib,nv3d)		!volume
 1000	  format(i5,6f10.5,e14.6)
	end do

	ns = nsect
	itypes(:) = isects(2,:)
	ic = count( itypes == 1 .or. itypes == 99 )
	!write(69,*) ns,ic
	!ns = ic		!only show type 1 (zeta bound) or 99 (internal)

	!if( bextra ) write(iu,'(a)') '# fluxes through sections'
	if( bextra ) write(iu,'(a)') '#   sections'
	if( bw ) write(iu,*) ns
	s1 = '# section  from    to'
	s2 = '         total      positive      negative'
	!     12345678901234567890123456789012345678901234567890
	if( bextra ) write(iu,'(a)') trim(s1)//trim(s2)
	flux2d(:,:) = fluxesg(0,:,:)
	do is=1,ns
	  itype = isects(2,is)
	  !if( itype /= 1 .and. itype /= 99 ) flux2d = 0.
	  !call flx_collect_2d(3*nsect,flux2d)
	  ib1 = isects(3,is)
	  ib2 = isects(4,is)
	  !write(69,*) is,ib1,ib2,itype
	  !write(iu,2000) ib1,ib2,(fluxes(0,ii,is),ii=1,3)
	  if( bw ) write(iu,2000) is,ib1,ib2,(flux2d(ii,is),ii=1,3)
 2000	  format(i9,2i6,3e14.6)
	end do

	if( nversbox < 6 ) then		!only write up to version 5
	ns = nbc_ob
	if( bextra ) write(iu,'(a)') '#        OBC'
	if( bw ) write(iu,*) ns
	if( bextra ) write(iu,'(a)') trim(s1)//trim(s2)

	if( ns > nsect ) stop 'error stop box_write_2d: nbc_ob>nsect'
	flux2d(:,1:nbc_ob) = fluxesg_ob(0,:,1:nbc_ob)

	do is=1,ns
	  itype = iscbnd(2,is)
	  ib1 = iscbnd(3,is)
	  ib2 = iscbnd(4,is)
	  if( itype .eq. 1 ) then	!z-boundary, ignore
	    ib1 = 0
	    ib2 = 0
	  end if
	  if( bw ) write(iu,2000) is,ib1,ib2,(flux2d(ii,is),ii=1,3)
	end do
	end if

	if( bflush ) call file_flush(iu)

	end

!******************************************************************

        subroutine box_write_3d(dtime,aline,nblayers,nslayers           &
     &                  ,fluxesg,fluxesg_ob                             &
     &                  ,nv3d,val3dg                                    &
     &                  ,nv2d,val2d                                     &
     &			,eta_act)

! writes 3d box values to file

	use basin
	use levels
	use box
	use shympi

	implicit none

	double precision dtime
	character*(*) aline
	integer nblayers(nbox)		!number of layers in box
	integer nslayers(nsect)		!number of layers in section
	double precision fluxesg(0:nlv_global,3,nsect) !mass fluxes box to box
	double precision fluxesg_ob(0:nlv_global,3,nbc_ob) !mass fluxes at OB
	integer nv3d
	double precision val3dg(0:nlv_global,nbox,nv3d)	!3d variables
	integer nv2d
	double precision val2d(nbox,nv2d)		!2d variables
	double precision eta_act(nbox)			!actual eta

	logical bw,bcheck
	integer ib,is,ib1,ib2,ii,itype,iv,iss
	integer lmax,l

	integer ifileo
	character*80 file,s1,s2
	integer, save :: iu = 0

	bcheck = .true.
	bw = .false.		!FIXME
	bw = .true.		!FIXME
	if( iu == 0 ) then
	  file = 'boxes_3d.txt'
	  iu = ifileo(0,file,'formatted','new')
	  if( iu <= 0 ) stop 'error stop boxes: opening file'
	  if( bextra ) call box_write_header(iu,3)
	end if

	if( bextra ) write(iu,'(a)') '# time record'
	if( bw ) write(iu,*) dtime,trim(aline)

        if( bextra ) write(iu,'(a)')                                    &
     &		'#      boxes  max_layers       nvars'
	if( bw ) write(iu,*) nbox,nlvdi,nv3d

	do ib=1,nbox
	  lmax = nblayers(ib)
	  s1 = '#   box layers  aver_eta   act_eta   bstress'
	  !                   12345678901234567890123456789012345678901234567890
	  if( bextra ) write(iu,'(a)') trim(s1)
          if( bw ) write(iu,1100) ib,lmax                               &
     &                          ,val2d(ib,1)                            &
     &                          ,eta_act(ib)                            &
     &                          ,val2d(ib,2)                            &
	  s1 = '# layer      temp      salt cur_speed        volume'
	!          12345678901234567890123456789012345678901234567890
	  if( bextra ) write(iu,'(a)') trim(s1)
	  do l=1,lmax
	    if( bw ) write(iu,1000) l,(val3dg(l,ib,iv),iv=1,nv3d)
 1100	    format(2i7,10f10.5)
 1000	    format(i7,3f10.5,e14.6)
	  end do
	end do

	!if( bextra ) write(iu,'(a)') '# fluxes through sections'
	if( bextra ) write(iu,'(a)') '#   sections  max_layers'
	if( bw ) write(iu,*) nsect,nlvdi

	s1 = '#    section        from          to      layers'
	!     12345678901234567890123456789012345678901234567890
	s2 = '# layer         total      positive      negative'

	do is=1,nsect
	  lmax = nslayers(is)
	  ib1 = isects(3,is)
	  ib2 = isects(4,is)
	  if( bextra ) write(iu,'(a)') trim(s1)
	  if( bw ) write(iu,2000) is,ib1,ib2,lmax
	  if( bextra ) write(iu,'(a)') trim(s2)
	  do l=1,lmax
	    if( bw ) write(iu,2100) l,(fluxesg(l,ii,is),ii=1,3)
	  end do
	end do

 2000	format(4i12)
 2100	format(i7,3e14.6)

	if( nversbox < 6 ) then		!only write up to version 5
	if( bextra ) write(iu,'(a)') '#        OBC  max_layers'
	if( bw ) write(iu,*) nbc_ob,nlvdi

	do is=1,nbc_ob
	  itype = iscbnd(2,is)
	  ib1 = iscbnd(3,is)
	  ib2 = iscbnd(4,is)
	  iss = iscbnd(5,is)		!section number
	  lmax = nslayers(iss)
	  if( itype .eq. 1 ) then	!z-boundary, ignore
	    ib1 = 0
	    ib2 = 0
	  end if
	  if( bextra ) write(iu,'(a)') trim(s1)
	  if( bw ) write(iu,2000) is,ib1,ib2,lmax
	  if( bextra ) write(iu,'(a)') trim(s2)
	  do l=1,lmax
	    if( bw ) write(iu,2100) l,(fluxesg_ob(l,ii,is),ii=1,3)
	  end do
	end do
	end if

	if( bflush ) call file_flush(iu)

	if( .not. bcheck ) return

	do is=1,nsect
	  lmax = nslayers(is)
	  if( any( fluxesg(lmax+1:nlv_global,:,is) /= 0. ) ) then
	    write(6,*) 'fluxes below last layer found...'
	    stop 'error stop box_write_3d: fluxes below last layer'
	  end if
	end do

	end

!******************************************************************

        subroutine box_write_vertical(dtime,aline,nblayers              &
     &			,nvv3d,valv3d)

! writes 3d interface box values to file (bottom interfaces)

	use basin
	use levels
	use box

	implicit none

	double precision dtime
	character*(*) aline
	integer nblayers(nbox)		!number of layers in box
	integer nvv3d
	double precision valv3d(0:nlvdi,nbox,nvv3d)

	logical bw
	integer ib,is,ib1,ib2,ii,itype,iv
	integer lmax,l,max

	integer ifileo
	character*80 file,s1
	integer, save :: iu = 0

	!          12345678901234567890123456789012345678901234567890

	bw = .false.		!FIXME
	bw = .true.		!FIXME
	if( iu == 0 ) then
	  file = 'boxes_vertical.txt'
	  iu = ifileo(0,file,'formatted','new')
	  if( iu <= 0 ) stop 'error stop boxes: opening file'
	  if( bextra ) call box_write_header(iu,4)
	end if

	if( bextra ) write(iu,'(a)') '# time record'
	if( bw ) write(iu,*) dtime,trim(aline)

	if( bextra ) write(iu,'(a)') '#      boxes  max_layers   nvars'
	if( bw ) write(iu,*) nbox,nlvdi,nvv3d

	valv3d(:,:,2) = -valv3d(:,:,2)		!make fluxes positive

	do ib=1,nbox
	  lmax = nblayers(ib)
	  s1 = '#   box layers'
	  if( bextra ) write(iu,'(a)') trim(s1)
	  if( bw ) write(iu,1100) ib,lmax
          s1 = '# layer flux_positive flux_negative'//                  &
     &			'     viscosity   diffusivity'
	  if( bextra ) write(iu,'(a)') trim(s1)
	  do l=1,lmax
	    if( bw ) write(iu,1000) l,(valv3d(l,ib,iv),iv=1,nvv3d)
 1100	    format(2i7)
 1000	    format(i7,4e14.6)
	  end do
	end do

	if( bflush ) call file_flush(iu)

	end

!******************************************************************
!******************************************************************
!******************************************************************
! box average routines
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine boxes_compute_area(area)

	use evgeom
	use basin
	use box
	use shympi

	implicit none

	double precision area(nbox)

	double precision vold(nbox)

	integer ib,ie,ii,k
	double precision vol,area3

	vold = 0.

	do ie=1,nel_unique
	  ib = iboxes(ie)
	  area3 = 4.*ev(10,ie)
	  vol = 3.*area3
	  vold(ib) = vold(ib) + vol
	end do

	area = vold

	end

!******************************************************************

	subroutine boxes_compute_volume(vola)

	use mod_layer_thickness
	use evgeom
	use levels
	use basin
	use box
	use shympi

	implicit none

	double precision vola(0:nlvdi,nbox)

	double precision vold(0:nlvdi,nbox)

	logical baver
	integer ib,ie,lmax,l,ii,k
	double precision vol,area3

	vold = 0.

	do ie=1,nel_unique
	  ib = iboxes(ie)
	  area3 = 4.*ev(10,ie)
	  lmax = ilhv(ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    do l=1,lmax
	      vol = area3 * hdknv(l,k)
	      vold(0,ib) = vold(0,ib) + vol
	      vold(l,ib) = vold(l,ib) + vol
	    end do
	  end do
	end do

	vola = vold

	end

!******************************************************************

	subroutine box_aver_eta(zev,val)

! computes average eta values for box

	use mod_geom_dynamic
	use evgeom
	use basin
	use box
	use shympi

	implicit none

        real zev(3,nel)
	double precision val(nbox)

	double precision vald(nbox)
	double precision vold(nbox)
	!double precision, allocatable :: vald(:)
	!double precision, allocatable :: vold(:)

	integer ib,ie,ii
	real vol,area3,z

	!write(6,*) 'gyeywu',my_id,nbox,nel,nel_unique
	!write(6,*) my_id,nbox,nel_unique
	!write(6,*) my_id,size(zev,1),size(zev,2)
	!write(6,*) my_id,size(val)
	!allocate(vald(nbox),vold(nbox))
	!stop

	vald = 0.
	vold = 0.

	do ie=1,nel_unique
	  !if( iwegv(ie) .eq. 0 ) then	!only for wet elements
	    ib = iboxes(ie)
	if( ib <= 0 ) stop 'internal error...'
	    area3 = 4.*ev(10,ie)
	    vol = area3
	    do ii=1,3
	      z = zev(ii,ie)
	      vald(ib) = vald(ib) + vol*z
	      vold(ib) = vold(ib) + vol
	    end do
	  !end if
	end do

	val = 0.
	where( vold > 0. ) val = vald / vold

	end

!******************************************************************

	subroutine box_3d_aver_vel(val)

! computes average velocity values (speed) for box

	use mod_layer_thickness
	use mod_hydro
	use evgeom
	use levels
	use basin
	use box
	use shympi

	implicit none

	double precision val(0:nlvdi,nbox)

	real h,u,v
	real velspeed(nlvdi)

	double precision vald(0:nlvdi,nbox)
	double precision vold(0:nlvdi,nbox)

	integer ib,ie,ii
	integer lmax,l,k
	real vol,area3

	velspeed = 0.

	vald = 0.
	vold = 0.

	do ie=1,nel_unique
	  ib = iboxes(ie)
	  area3 = 4.*ev(10,ie)
	  lmax = ilhv(ie)
	  do l=1,lmax
	    h = hdenv(l,ie)
	    u = utlnv(l,ie)/h
	    v = vtlnv(l,ie)/h
	    velspeed(l) = sqrt(u*u+v*v)
	  end do
	  do ii=1,3
	    k = nen3v(ii,ie)
	    do l=1,lmax
	      vol = area3 * hdknv(l,k)
	      vald(0,ib) = vald(0,ib) + vol*velspeed(l)
	      vold(0,ib) = vold(0,ib) + vol
	      vald(l,ib) = vald(l,ib) + vol*velspeed(l)
	      vold(l,ib) = vold(l,ib) + vol
	    end do
	  end do
	end do

	val = 0.
	where( vold > 0. ) val = vald / vold

	end

!******************************************************************

	subroutine box_3d_aver_vertical(scalar,val,iaver)

! computes average scalar values defined on interface (verticals) for box

	use mod_depth
	use evgeom
	use levels
	use basin
	use box

	implicit none

	real scalar(0:nlvdi,nkn)
	double precision val(0:nlvdi,nbox)
	integer iaver			!average or only accumulate

	double precision vald(0:nlvdi,nbox)
	double precision vold(0:nlvdi,nbox)

	logical baver
	integer ib,ie,lmax,l,ii,k
	double precision vol,area3

	baver = iaver > 0

	vald = 0.
	vold = 0.

	do ie=1,nel
	  ib = iboxes(ie)
	  area3 = 4.*ev(10,ie)
	  vol = area3			!use area for averaging
	  lmax = ilhv(ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    do l=0,lmax
	      vald(l,ib) = vald(l,ib) + vol*scalar(l,k)
	      vold(l,ib) = vold(l,ib) + vol
	    end do
	  end do
	end do

	if( .not. baver ) then
	  val = vald	!only accumulate
	else
	  val = 0.
	  where( vold > 0. ) val = vald / vold
	end if

	end

!******************************************************************

	subroutine box_3d_aver_scalar(scalar,val,iaver)

! computes average scalar values for box

	use mod_layer_thickness
	use evgeom
	use levels
	use basin
	use box
	use shympi

	implicit none

	real scalar(nlvdi,nkn)
	double precision val(0:nlvdi,nbox)
	integer iaver			!average or only accumulate

	double precision vald(0:nlvdi,nbox)
	double precision vold(0:nlvdi,nbox)

	logical baver
	integer ib,ie,lmax,l,ii,k
	double precision vol,area3

	baver = iaver > 0

	vald = 0.
	vold = 0.

	do ie=1,nel_unique
	  ib = iboxes(ie)
	  area3 = 4.*ev(10,ie)
	  lmax = ilhv(ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    do l=1,lmax
	      vol = area3 * hdknv(l,k)
	      vald(0,ib) = vald(0,ib) + vol*scalar(l,k)
	      vold(0,ib) = vold(0,ib) + vol
	      vald(l,ib) = vald(l,ib) + vol*scalar(l,k)
	      vold(l,ib) = vold(l,ib) + vol
	    end do
	  end do
	end do

	if( .not. baver ) then
	  val = vald	!only accumulate
	else
	  val = 0.
	  where( vold > 0. ) val = vald / vold
	end if

	end

!******************************************************************

	subroutine box_2d_aver_scalar(scalar,val)

! computes accumulation of 2d scalar values for box

	use evgeom
	use basin
	use box
	use shympi

	implicit none

	real scalar(nkn)
	double precision val(nbox)

	double precision vald(nbox)
	double precision vold(nbox)

	integer ib,ie,ii,k
	double precision vol,area3

	vald = 0.
	vold = 0.

	do ie=1,nel_unique
	  ib = iboxes(ie)
	  area3 = 4.*ev(10,ie)
	  vol = area3
	  do ii=1,3
	    k = nen3v(ii,ie)
	    !if( k > nkn_unique ) cycle
	    vald(ib) = vald(ib) + vol*scalar(k)
	    vold(ib) = vold(ib) + vol
	  end do
	end do

	val = 0.
	where( vold > 0. ) val = vald / vold
	!val = vald

	end

!******************************************************************

	subroutine box_3d_aver_scalar_box(ibox,scalar,val)

! computes average scalar values for box (only for box ibox - debug)

	use mod_layer_thickness
	use evgeom
	use levels
	use basin
	use box

	implicit none

	integer ibox
	real scalar(nlvdi,nkn)
	real val(0:nlvdi,nbox)

	double precision vald(0:nlvdi,nbox)
	double precision vold(0:nlvdi,nbox)

	integer ib,ie,lmax,l,ii,k
	double precision vol,area3

	vald = 0.
	vold = 0.

	write(167,*) '-----------------------------------'

	do ie=1,nel
	  ib = iboxes(ie)
	  if( ib /= ibox ) cycle
	  area3 = 4.*ev(10,ie)
	  lmax = ilhv(ie)	!FIXME
	  do ii=1,3
	    k = nen3v(ii,ie)
	    do l=1,lmax
	      write(167,*) ib,ie,k,l,scalar(l,k)
	      vol = area3 * hdknv(l,k)	!FIXME
	      vald(0,ib) = vald(0,ib) + vol*scalar(l,k)
	      vold(0,ib) = vold(0,ib) + vol
	      vald(l,ib) = vald(l,ib) + vol*scalar(l,k)
	      vold(l,ib) = vold(l,ib) + vol
	    end do
	  end do
	end do

	val = 0.
	where( vold > 0. ) val = vald / vold
	write(167,*) val(:,ibox)

	end

!******************************************************************
!******************************************************************
!******************************************************************
! time init/accum/aver routines for boxes
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine boxes_2d_init(n,nv,val)

	implicit none

	integer n,nv
	double precision val(n,nv)

	val = 0.

	end

!******************************************************************

	subroutine boxes_2d_accum(n,dt,area,val,value)

	implicit none

	integer n
	real dt
	double precision area(n)	!accumulation area
	double precision val(n)		!accumulator
	double precision value(n)	!value to accumulate

	val = val + value * dt * area

	end

!******************************************************************

	subroutine boxes_2d_aver(n,nv,bdtarea,val)

	implicit none

	integer n,nv
	double precision bdtarea(n)
	double precision val(n,nv)

	integer iv
	double precision area(n)
	double precision vbox(n)

	area = bdtarea
	call gather_sum_d2(n,area)

	do iv=1,nv
	  vbox = val(:,iv)
	  call gather_sum_d2(n,vbox)
	  vbox = vbox / area
	  val(:,iv) = vbox
	end do

	end

!******************************************************************

	subroutine boxes_3d_init(nlvddi,n,nv,val)

	implicit none

	integer nlvddi,n,nv
	double precision val(0:nlvddi,n,nv)

	val = 0.

	end

!******************************************************************

	subroutine boxes_3d_accum(nlvddi,n,dt,bvol3d,val,value)

	implicit none

	integer nlvddi,n
	real dt
	double precision bvol3d(0:nlvddi,n)
	double precision val(0:nlvddi,n)
	double precision value(0:nlvddi,n)

	val = val + value * dt * bvol3d

	end

!******************************************************************

	subroutine boxes_3d_sum(nlvddi,n,val,valg)

	use shympi

	implicit none

	integer, intent(in) :: nlvddi,n
	double precision, intent(in) :: val(0:nlvddi,n)
	double precision, intent(out) :: valg(0:nlv_global,n)

	integer iv
	integer l,i
	double precision vbox(0:nlv_global,n)

	if( .not. bmpi ) then
	  valg = val
	  return
	end if

	vbox = 0.

	vbox(0:nlvddi,:) = val(:,:)
	call gather_sum_d3(nlv_global,n,vbox)
	valg = vbox

	end

!******************************************************************

	subroutine boxes_3d_aver(nlvddi,n,nv,vol,val,valg)

	use shympi

	implicit none

	integer, intent(in) :: nlvddi,n,nv
	double precision, intent(in) :: vol(0:nlvddi,n)
	double precision, intent(in) :: val(0:nlvddi,n,nv)
	double precision, intent(out) :: valg(0:nlv_global,n,nv)

	integer iv
	integer l,i
	double precision volume(0:nlv_global,n)
	double precision vbox(0:nlv_global,n)

	!if( .not. bmpi ) then
	!  valg = 0.
	!  where( vol > 0 ) valg = val / vol
	!  return
	!end if

	volume = 0.
	vbox = 0.

	volume(0:nlvddi,:) = vol
	if( bmpi ) call gather_sum_d3(nlv_global,n,volume)

	do iv=1,nv
	  vbox(0:nlvddi,:) = val(:,:,iv)
	  if( bmpi ) call gather_sum_d3(nlv_global,n,vbox)
	  where( volume > 0 ) vbox = vbox / volume
	  valg(:,:,iv) = vbox
	end do

	end

!******************************************************************
!******************************************************************
!******************************************************************
! open boundary handling routines
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine box_ob_init

! sets up open boundary inputs

	use basin
	use box
	use shympi

	implicit none

	integer nbc,ibc,itype
	integer nk,i,k,ib,ibk
	integer ierr
	integer iauxv(nkn)

	integer nbnds,itybnd,nkbnds,kbnds,ndim
	integer ifrom,ito
	integer ipext

	ierr = 0
	kfluxm = kfluxm + 1			!GGU_LAST
	kflux(kfluxm) = 0
	ndim = nkn

	nbc = nbnds()

	do ibc = 1,nbc
	  !write(6,*) ibc,nbc
	  iscbnd(:,ibc) = 0
	  itype = itybnd(ibc)
	  nk = nkbnds(ibc)
	  if( nk <= 0 ) cycle			!no boundary
	  ib = 0
	  call irbnds(ibc,ndim,nk,iauxv)
	  where( iauxv == 0 ) iauxv = -1	!flag non existing nodes
	  do i=1,nk
	    k = iauxv(i)
	    if( k <= 0 ) cycle
	    ibk = ikboxes(k)
	    if( ibk .le. 0 ) goto 99	!node not in unique box
	    if( ib .eq. 0 ) ib = ibk
	    if( ibk .ne. ib ) goto 98	!boundary nodes in different boxes
	  end do

	  ifrom = -ibc
	  ito = ib
	  ito = shympi_max(ito)
	  itype = shympi_max(itype)

	  nsect = nsect + 1

	  iscbnd(1,ibc) = nk		!total number of boundary nodes
	  iscbnd(2,ibc) = itype		!type of boundary
	  iscbnd(3,ibc) = ifrom		!from box
	  iscbnd(4,ibc) = ito		!to box
	  iscbnd(5,ibc) = nsect		!section number

	  if( nsect .gt. nscboxdim ) goto 95
	  call box_insert_section(nsect,nk,itype,ifrom,ito,iauxv,ierr)
	  if( kfluxm .gt. nfxboxdim ) goto 96
	  if( ierr /= 0 ) goto 93
	end do

	!write(630+my_id,'(5i5)') (ibc,iscbnd(:,ibc),ibc=1,nbc)
	!flush(630+my_id)

	!kfluxm = kfluxm - 1

	return
   93	continue
	write(6,*) ierr
	stop 'error stop box_ob_init: unspecified error'
   95	continue
	write(6,*) nsect,nscboxdim
	stop 'error stop box_ob_init: nscboxdim'
   96	continue
	write(6,*) kfluxm,nfxboxdim
	stop 'error stop box_ob_init: nfxboxdim'
   98	continue
	write(6,*) 'boxes: ',ib,ibk,'  node: ',ipext(k)
	stop 'error stop box_ob_init: different boxes for nodes'
   99	continue
	write(6,*) 'boxes: ',ibk,'  node: ',ipext(k)
	stop 'error stop box_ob_init: not unique box for node'
	end

!******************************************************************

	subroutine box_ob_compute(nbc_ob,nlvddi,fluxes_ob)

! computes open boundary inputs

	use shympi

	implicit none

	integer nbc_ob
	integer nlvddi
	double precision fluxes_ob(0:nlvddi,3,nbc_ob) !discharges at open bound

	integer ibc,nbc,l,lmax
	real val
	real flux3d(nlvddi)
	real flux2d

	nbc = nbc_ob
	lmax = nlvddi
	fluxes_ob = 0.

	do ibc = 1,nbc
	  call get_discharge_3d(ibc,flux3d,flux2d)
	  val = flux2d
	  if( val .ge. 0 ) then
	    fluxes_ob(0,1,ibc) = val
	    fluxes_ob(0,2,ibc) = val
	    fluxes_ob(0,3,ibc) = 0.
	  else
	    fluxes_ob(0,1,ibc) = val
	    fluxes_ob(0,2,ibc) = 0.
	    fluxes_ob(0,3,ibc) = -val
	  end if
	  do l=1,lmax
	    val = flux3d(l)
	    if( val .ge. 0 ) then
	      fluxes_ob(l,1,ibc) = val
	      fluxes_ob(l,2,ibc) = val
	      fluxes_ob(l,3,ibc) = 0.
	    else
	      fluxes_ob(l,1,ibc) = val
	      fluxes_ob(l,2,ibc) = 0.
	      fluxes_ob(l,3,ibc) = -val
	    end if
	  end do
	end do

	end

!******************************************************************
!******************************************************************
!******************************************************************
! meteo routines
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine boxes_meteo_accum(nbox,dt,nvmet,barea,valmet,ievap)

! meteo parameter accumulation
!
! rain and evaporation are in [m/s]
! before writing they are converted to [mm/day]

	use mod_meteo

	implicit none

	integer nbox
	real dt
	integer nvmet
	double precision barea(nbox)
	double precision valmet(nbox,nvmet)
	integer ievap

        real zconv
        parameter( zconv = 1.e-3 / 86400. )	!convert from [mm/day] to [m/s]

	integer i,ip
	real econv,rconv
	double precision val(nbox)

	rconv = 1. / zconv
	econv = 1. / zconv
	if( ievap .le. 0 ) econv = 0.

	ip = 1
	call box_2d_aver_scalar(metrad,val)
	call boxes_2d_accum(nbox,dt,barea,valmet(:,ip),val)

	ip = 2
	call box_2d_aver_scalar(mettair,val)
	call boxes_2d_accum(nbox,dt,barea,valmet(:,ip),val)

	ip = 3
	call box_2d_aver_scalar(methum,val)
	call boxes_2d_accum(nbox,dt,barea,valmet(:,ip),val)

	ip = 4
	call box_2d_aver_scalar(metws,val)
	call boxes_2d_accum(nbox,dt,barea,valmet(:,ip),val)

	ip = 5
	call box_2d_aver_scalar(metcc,val)
	call boxes_2d_accum(nbox,dt,barea,valmet(:,ip),val)

	ip = 6
	call box_2d_aver_scalar(metrain,val)
	val = val * rconv
	call boxes_2d_accum(nbox,dt,barea,valmet(:,ip),val)

	ip = 7
	call box_2d_aver_scalar(evapv,val)
	val = val * econv
	call boxes_2d_accum(nbox,dt,barea,valmet(:,ip),val)

	ip = 8
	call box_2d_aver_scalar(metice,val)
	call boxes_2d_accum(nbox,dt,barea,valmet(:,ip),val)

	if( ip /= nvmet ) then
	  write(6,*) 'ip,nvmet: ',ip,nvmet
	  stop 'error stop boxes_meteo_accum: parameter mismatch'
	end if

	end

!******************************************************************

	subroutine box_write_meteo(dtime,aline,nvmet,valmet)

! parameters are:
!	1	srad
!	2	tair
!	3	rhum
!	4	wspeed
!	5	cc
!	6	rain
!	7	evap
!	8	ice

	use basin
	use box
	use shympi

	implicit none

	double precision dtime
	character*(*) aline
	integer nvmet
	double precision valmet(nbox,nvmet)			!meteo variables

	logical bw
	integer ib,iv

	integer ifileo
	character*80 file,dline,s1,s2
	double precision aux2d(nbox)			!meteo variables
	integer, save :: iu = 0

	file = 'boxes_meteo.txt'
	if( shympi_is_master() ) then
	  if( iu == 0 ) then
	    iu = ifileo(0,file,'formatted','new')
	    if( iu <= 0 ) stop 'error stop boxes: opening file'
	    if( bextra ) call box_write_header(iu,1)
	  end if
	  bw = .true.
	else
	  iu = 987
	  bw = .false.
	  bextra = .false.
	  bflush = .false.
	end if

	if( bextra ) write(iu,'(a)') '# time record'
	if( bw ) write(iu,*) dtime,trim(aline)

	if( bextra ) write(iu,'(a)') '#      boxes       nvars'
	if( bw ) write(iu,*) nbox,nvmet
	s1 = ' box     srad    tair    rhum  wspeed      cc'
	s2 = '    rain    evap     ice' 
	if( bextra ) write(iu,'(a)') '#'//trim(s1)//trim(s2)

	do ib=1,nbox
	  if( bw ) write(iu,1002) ib,(valmet(ib,iv),iv=1,nvmet)
	end do

	if( bflush ) call file_flush(iu)

 1000	format(i14,20g12.4)
 1001	format(i4,1x,5f8.3,2g12.4,5f8.3)
 1002	format(i5,1x,20f8.3)
	end

!******************************************************************

	subroutine box_write_header(iu,ftype)

	use box

	implicit none

	integer iu,ftype

!                             1234567890 
	write(iu,'(a,3i10)') '#      idbox file_type     nvers'
	write(iu,'(a,3i10)') '# ',idbox,ftype,nversbox

	end

!******************************************************************
!******************************************************************
!******************************************************************
! mass balance routines
!******************************************************************
!******************************************************************
!******************************************************************

        subroutine boxes_mass_balance_2d(dtbox,bvol2d                   &
     &					,fluxesg,fluxesg_ob,voldif)

! computes mass balance

	use basin
	use levels
	use box
	use shympi

	implicit none

	double precision dtbox
	double precision bvol2d(nbox)
	double precision fluxesg(0:nlv_global,3,nsect)
	double precision fluxesg_ob(0:nlv_global,3,nbc_ob)
	double precision voldif(nbox)		!volume difference due to eta

	character*80 string
	logical bdebug,bw
	integer ib,is,ib1,ib2,iu,l
	real volf,vole,volb,vola,vold,volr
	real errmax,errv(nbox),errd(nbox)
	real eps
	double precision vol(nbox)
	double precision flux

	eps = eps2d
	if( eps == 0. ) return

	bdebug = .false.
	bdebug = .true.
	bw = ( my_id == 0 )

	vol = 0.

	!write(601,*) nsect
	do is=1,nsect
	  ib1 = isects(3,is)
	  ib2 = isects(4,is)
	  flux = dtbox * fluxesg(0,1,is)
	  !write(601,*) ib1,ib2,flux,fluxes(0,1,is)
	  if( ib1 > 0 ) vol(ib1) = vol(ib1) - flux
	  if( ib2 > 0 ) vol(ib2) = vol(ib2) + flux
	end do

	!write(601,*) nbc_ob
	do is=1,nbc_ob
	  ib1 = iscbnd(3,is)
	  ib2 = iscbnd(4,is)
	  flux = dtbox * fluxesg_ob(0,1,is)
	  !write(601,*) ib1,ib2,flux,fluxes_ob(0,1,is)
	  !if( ib1 > 0 ) vol(ib1) = vol(ib1) - flux
	  !if( ib2 > 0 ) vol(ib2) = vol(ib2) + flux
	end do

	errmax = 0.
	do ib=1,nbox
	  volb = bvol2d(ib)			!total volume of box
	  if( volb <= 0. ) cycle
	  volf = vol(ib)			!volume change by fluxes
	  vole = voldif(ib)			!volume change by eta
	  vola = max(abs(volf),abs(vole))
	  vold = abs(volf-vole)
	  !write(601,*) ib,volf,vole,vold
	  volr = 0.
	  !volr = vold/vola
	  volr = vold/volb	!use total volume
	  errv(ib) = volr
	  errd(ib) = vold
	  if( volr .gt. eps ) then
	    if(bw) write(6,*) '*** mass balance error in box ',ib,my_id
	    if(bw) write(6,*) '***',ib,volf,vole,vold,volr
	  end if
	  errmax = max(errmax,volr)
	  !write(115,*) ib,volf,vole,vold,volr
	end do

	string = 'boxes mass balance errmax 2d: '

	if( errmax .gt. eps ) then
	  if( bw ) write(6,*) '*** ',trim(string),errmax,eps
	else if( bdebug .and. bw ) then
	  write(6,*) trim(string),errmax,eps
	end if

	if( bmasserror .and. bw ) then
	  write(601,*) trim(string),errmax,eps
	  write(601,*) nbox
	  do ib=1,nbox
	    write(601,*) ib,errv(ib),errd(ib)
	  end do
	end if

	end

!******************************************************************

        subroutine boxes_mass_balance_3d(dtbox,bvol2d                   &
     &                                  ,nblayers,nslayers              &
     &                                  ,fluxesg,fluxesg_ob             &
     &					,voldif,wflux)

! computes mass balance

	use basin
	use levels
	use box
	use shympi

	implicit none

	double precision dtbox
	double precision bvol2d(nbox)
	integer nblayers(nbox)
	integer nslayers(nsect)		!number of layers in section
	double precision fluxesg(0:nlv_global,3,nsect)
	double precision fluxesg_ob(0:nlv_global,3,nbc_ob)
	double precision voldif(nbox)			!volume difference due to eta
	double precision wflux(0:nlv_global,nbox)	!vertical fluxes [m**3/s]

	character*80 string
	logical bdebug,bw
	integer ierr
	integer ib,is,ib1,ib2
	integer nb1,nb2,nbs
	integer l,lmax,ltot
	real volf,vole,volb,vola,vold,volr,volv
	real wtop
	real errmax,errv(nbox),errd(nbox)
	real eps
	double precision vol(nlvdi,nbox)
	double precision flux

	eps = eps3d
	if( eps == 0. ) return

	bdebug = .false.
	bdebug = .true.
	bw = ( my_id == 0 )

	ltot = 0
	do ib=1,nbox
	  ltot = max(ltot,nblayers(ib))
	end do

	vol = 0.
	ierr = 0
	string = '   is  ib1  ib2  nb1  nb2    l          flux'

	do is=1,nsect
	  nbs = nslayers(is)
	  ib1 = isects(3,is)
	  ib2 = isects(4,is)
	  nb1 = 0
	  if( ib1 > 0 ) nb1 = nblayers(ib1)
	  nb2 = nblayers(ib2)
	  lmax = nbs
	  do l=1,lmax
	    flux = dtbox * fluxesg(l,1,is)
	    if( ib1 > 0 ) vol(l,ib1) = vol(l,ib1) - flux
	    if( ib2 > 0 ) vol(l,ib2) = vol(l,ib2) + flux
	    if( bdebug .and. abs(flux) > 0. ) then
	      if( l > nb1 .or. l > nb2 ) then
	       if( ib1 > 0 .and. nb1 > 0 ) then		
		if( ierr == 0 ) write(6,'(19x,a)') trim(string)
                write(6,'(a,6i5,g14.4)') 'non existing flux: '          &
     &                          ,is,ib1,ib2                             &
     &				,nb1,nb2,l,flux
		ierr = ierr + 1
	       end if
	      end if
	    end if
	  end do
	end do

	do is=1,nbc_ob
	  ib1 = iscbnd(3,is)
	  ib2 = iscbnd(4,is)
	  flux = dtbox * fluxesg_ob(1,1,is)
	  !if( ib1 > 0 ) vol(1,ib1) = vol(1,ib1) - flux
	  !if( ib2 > 0 ) vol(1,ib2) = vol(1,ib2) + flux
	end do

	do ib=1,nbox
	  wtop = 0.
	  wflux(ltot,ib) = wtop		!bottom
	  do l=ltot,1,-1
	    wtop = wtop + vol(l,ib)
	    wflux(l-1,ib) = -wtop / dtbox
	  end do
	end do

	errmax = 0.
	do ib=1,nbox
	  lmax = nblayers(ib)
	  do l=lmax+1,ltot
	    if( bdebug .and. vol(l,ib) /= 0. ) then
              write(6,'(a,3i5,g14.4)') 'non existing box: '             &
     &				,ib,lmax,l,vol(l,ib)
	    end if
	  end do
	  
	  volb = bvol2d(ib)
	  if( volb <= 0. ) cycle

	  do l=1,ltot
	    volf = vol(l,ib)
	    volv = dtbox * (wflux(l-1,ib) - wflux(l,ib))
	    vold = abs(volf + volv)
	    volr = vold/volb
	    !write(116,'(2i5,4g14.4)') ib,l,volf,volv,vold,volr
	    errmax = max(errmax,abs(volr))
	  end do

	  wtop = -dtbox * wflux(0,ib)
	  vole = voldif(ib)
	  vold = abs(wtop - vole)
	  volr = vold/volb
	  errv(ib) = volr
	  errd(ib) = vold
	  !write(116,'(a,i5,4g14.4)') 'first layer: ',ib,wtop,vole,vold,volr
	  errmax = max(errmax,abs(volr))

	end do

	string = 'boxes mass balance errmax 3d: '
	if( errmax .gt. eps ) then
	  if( bw ) write(6,*) '*** ',trim(string),errmax,eps
	else if( bdebug .and. bw ) then
	  write(6,*) trim(string),errmax,eps
	end if

	if( bmasserror .and. bw ) then
	  write(601,*) trim(string),errmax,eps
	  write(601,*) nbox
	  do ib=1,nbox
	    write(601,*) ib,errv(ib),errd(ib)
	  end do
	end if

	end

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine box_flx_file_open(nvar,da_out)

        use box
	use box_arrays

        implicit none

        integer nvar
	double precision da_out(4)

	integer i
	!character*80 chflx(nsect)
	character*80, allocatable :: chflx(:)
	character*80 string

	allocate(chflx(nsect))

	do i=1,nsect
	  write(string,'(a,i5)') 'internal box section ',i
	  chflx(i) = string
	end do

!-----------------------------------------------------------------
! normal call
!-----------------------------------------------------------------

        call flux_file_open('.box.flx',da_out,nvar,nsect,kfluxm         &
     &                          ,kflux_ext,nslayers,chflx)

	end

!******************************************************************

	subroutine box_flx_collect(nlvddi,flux_local,flux_global)

	use box
	use box_arrays
	use shympi

	implicit none

	integer nlvddi
        double precision flux_local(0:nlvddi,3,nsect)
        double precision flux_global(0:nlv_global,3,nsect)

        integer n

	n = 3*nsect

	call flx_collect_3d(nlvddi,n,flux_local,flux_global)

	end

!******************************************************************

	subroutine box_flx_write(atime,nbflx,ivar,flux_global)

	use box
	use box_arrays
	use shympi

	implicit none

        double precision atime
	integer nbflx
        integer ivar
        double precision flux_global(0:nlv_global,3,nsect)

        integer nl,ns

        nl = nlv_global
        ns = nsect

        call flux_write(nbflx,atime,ivar,nl,ns                          &
     &                          ,nslayers,flux_global)

	end

!******************************************************************

	subroutine box_debug_2d(text,iubase,n,val)

	use shympi

	implicit none

	character*(*) text
	integer iubase,n
	double precision val(n)

	integer iuaux,ib
	double precision valtot

	!return

	valtot = shympi_sum(val)
	
	iuaux = iubase + my_id
	write(iuaux,*) '--- ',trim(text),' ---',valtot
	do ib=1,n
	  write(iuaux,*) ib,val(ib)
	end do
	flush(iuaux)

	end

!******************************************************************

	subroutine box_debug_3d(text)

	use box
	use box_arrays
	use levels

	implicit none

	character*(*) text

	integer iu,ib,iv
	integer lmax,l

	return

	iu = 456

	write(iu,*) '==================================='
	write(iu,*) trim(text)
	write(iu,*) '==================================='

	do ib=1,nbox
	  lmax = nblayers(ib)
	  write(iu,*) ib,lmax
	  do l=0,lmax
	    write(iu,'(i4,4e15.4)') l,(val3d(l,ib,iv),iv=1,nv3d)
	  end do
	end do

	write(iu,*) 'nblayers: ',nblayers
	write(iu,*) 'nslayers: ',nslayers

	end

!******************************************************************

	subroutine assert_values

! assert correctness of values

	use box
	use box_arrays

	implicit none

	logical bw
	integer ib,lmax,l,iv,is
	integer iu
	double precision val,val0,dval,rval
	double precision fmax,fsum,ftot
	double precision vals(3),vals0(3),avals(3),dvals(3),rvals(3)
	double precision eps
	character*80 text

	bw = .true.
	bw = .false.
	iu = 456

!-------------------------------------------------
! check volume
!-------------------------------------------------

	eps = 1.e-10
	text = 'checking volume'
	iv = 0

	if( bw ) write(iu,*) 'checking variable ',iv
	do ib=1,nbox
	  lmax = nblayers(ib)
	  val = 0.
	  do l=1,lmax
	    val = val + bvol3d(l,ib)
	  end do
	  val0 = bvol3d(0,ib)
	  dval = abs(val-val0)
	  rval = dval / val
	  if( bw ) write(iu,2000) iv,ib,val,val0,dval,rval
	  if( rval > eps ) goto 99
	end do

!-------------------------------------------------
! check 3d scalars
!-------------------------------------------------

	eps = 1.e-8
	text = 'checking val3d'

	do iv=1,nv3d
	  if( bw ) write(iu,*) 'checking variable ',iv
	  do ib=1,nbox
	    lmax = nblayers(ib)
	    val = 0.
	    do l=1,lmax
	      val = val + val3d(l,ib,iv)
	    end do
	    val0 = val3d(0,ib,iv)
	    dval = abs(val-val0)
	    rval = 0.
	    if( val > 0. ) rval = dval / val
	    if( bw ) write(iu,2000) iv,ib,val,val0,dval,rval
	    if( rval > eps ) goto 99
	  end do
	end do

!-------------------------------------------------
! check fluxes
!-------------------------------------------------

	eps = 1.e-8
	text = 'checking fluxes'

	do is=1,nsect
	  if( bw ) write(iu,*) 'checking flux ',is
	  lmax = nslayers(is)
	  do l=1,lmax
	    fmax = max(fluxes(l,2,is),fluxes(l,3,is))
	    fsum = fluxes(l,2,is)-fluxes(l,3,is)
	    ftot = fluxes(l,1,is)
	    dval = abs(fsum-ftot)
	    rval = 0.
	    if( fmax > 0. ) rval = dval / fmax
	    !if( bw ) write(iu,2000) l,is,fsum,ftot,dval,rval
	    if( rval > eps ) goto 98
	  end do
	end do

!-------------------------------------------------
! check sections
!-------------------------------------------------

	eps = 1.e-8
	text = 'checking sections'

	do is=1,nsect
	  if( bw ) write(iu,*) 'checking section ',is
	  lmax = nslayers(is)
	  vals = 0.
	  do l=1,lmax
	    vals(:) = vals(:) + fluxes(l,:,is)
	  end do
	  vals0(:) = fluxes(0,:,is)
	  dvals = abs(vals-vals0)
	  avals = abs(vals)
	  rvals = 0.
	  where( avals > 0 ) rvals = dvals / avals 
	  val = maxval( vals )
	  val0 = maxval( vals0 )
	  dval = maxval( dvals )
	  rval = maxval( rvals )
	  if( bw ) write(iu,2000) 0,is,val,val0,dval,rval
	  if( rval > eps ) goto 97
	end do

!-------------------------------------------------
! end of routine
!-------------------------------------------------

	return
 1000	format(i5,4e14.5)
 2000	format(2i5,4e14.5)
   97	continue
	write(6,*) trim(text)
	write(6,*) 'assert error for eps = ',eps
	write(6,2000) 0,is,val,val0,dval,rval
	stop 'error stop assert_values'
   98	continue
	write(6,*) trim(text)
	write(6,*) 'assert error for eps = ',eps
	write(6,2000) l,is,fsum,ftot,dval,rval
	stop 'error stop assert_values'
   99	continue
	write(6,*) trim(text)
	write(6,*) 'assert error for eps = ',eps
	write(6,2000) iv,ib,val,val0,dval,rval
	stop 'error stop assert_values'
	end

!******************************************************************

	subroutine write_matrix_boxes

	use basin
	use levels
	use box
	use box_arrays
	use shympi

	implicit none

	integer ia,ie,ib,iu
	integer nblayers_domain(nbox,n_threads)
	integer ibox(nbox),lbox(nbox),header(nbox)
	integer ibox_domain(nbox,n_threads)
	integer lbox_domain(nbox,n_threads)
	integer nlv_domain(n_threads)

	do ib=1,nbox
	  header(ib) = mod(ib,10)
	end do

	nlv_domain = 0
	call shympi_gather(nlv,nlv_domain)

	nblayers_domain = 0
	call shympi_gather(nblayers,nblayers_domain)

	ibox = 0
	lbox = 0
	do ie=1,nel
	  ib = iboxes(ie)
	  if( ib < 1 .or. ib > nbox ) stop 'error stop dddddd'
	  ibox(ib) = 1
	  lbox(ib) = max(lbox(ib),ilhv(ie))
	end do

	ibox_domain = 0
	call shympi_gather(ibox,ibox_domain)
	lbox_domain = 0
	call shympi_gather(lbox,lbox_domain)

	if( my_id /= 0 ) return

	open(iu,file='boxes_matrix.txt',status='unknown',form='formatted')

	write(iu,1000) 0,header(:)
	write(iu,*) '---------------------------------'

	do ia=1,n_threads
	  write(iu,1000) ia,nblayers_domain(:,ia)
	end do
	
	write(iu,*) '---------------------------------'
	do ia=1,n_threads
	  write(iu,1000) ia,ibox_domain(:,ia)
	end do

	write(iu,*) '---------------------------------'
	do ia=1,n_threads
	  write(iu,1000) ia,lbox_domain(:,ia)
	end do

	write(iu,*) '---------------------------------'
	do ia=1,n_threads
	  write(iu,1000) ia,nlv_domain(ia)
	end do

	write(iu,*) '---------------------------------'
	close(iu)

 1000	format(i5,20i3)
	end

!******************************************************************

	subroutine flx_print(iu,text,nsect,nlvddi,fluxes)

	implicit none

	integer iu
	character*(*) text
	integer nsect,nlvddi
	double precision fluxes(0:nlvddi,3,nsect)

	integer is,l

	write(iu,*) '-------------'
	write(iu,*) trim(text),nsect,nlvddi
	write(iu,*) '-------------'

	do is=1,nsect
	  write(iu,1000) is,0,fluxes(0,:,is)
	  do l=1,nlvddi
	    !write(iu,1000) is,l,fluxes(l,:,is)
	  end do
	end do

	flush(iu)

 1000	format(2i5,3g16.8)
	end

!******************************************************************

	subroutine print_2d(nbox,barea,bvol2d,eta_act,eta_old)

	use box, only : b777
	use shympi

	implicit none

	integer nbox
	double precision barea(nbox)
	double precision bvol2d(nbox)
	double precision eta_act(nbox)
	double precision eta_old(nbox)

	integer iu,ib

	if( .not. b777 ) return
	if( my_id /= 0 ) return

	iu = 777
	write(iu,*) 3
	write(iu,*) '==========================='
	write(iu,*) 'print_2d:'
	write(iu,*) '==========================='
	write(iu,*) nbox,2,5

	do ib=1,nbox
	  write(iu,1000) ib,barea(ib),bvol2d(ib),eta_act(ib),eta_old(ib)
	end do

	flush(iu)

	return
 1000	format(i5,4g18.10)
	end

!******************************************************************

	subroutine print_3d_flx(nsect,nbc_ob,dtbox,fluxesg,fluxesg_ob)

	use box, only : b777
	use shympi

	implicit none

	integer nsect,nbc_ob
	double precision dtbox
	double precision fluxesg(0:nlv_global,3,nsect)
	double precision fluxesg_ob(0:nlv_global,3,nbc_ob)

	logical b3d
	integer is,l,iu,nm

	b3d = .false.
	nm = 1
	if( b3d ) nm = nlv_global + 1

	if( .not. b777 ) return
	if( my_id /= 0 ) return

	iu = 777
	write(iu,*) 3
	write(iu,*) '==========================='
	write(iu,*) 'print_3d_flx: ',dtbox
	write(iu,*) '==========================='
	write(iu,*) nsect*nm,3,5

	do is=1,nsect
	  write(iu,1000) is,0,fluxesg(0,:,is)
	  do l=1,nlv_global
	    if( b3d ) write(iu,*) is,l,fluxesg(l,:,is)
	  end do
	end do

	write(iu,*) 1
	write(iu,*) '---------------------------'
	write(iu,*) nbc_ob*nm,3,5

	do is=1,nbc_ob
	  write(iu,1000) is,0,fluxesg_ob(0,:,is)
	  do l=1,nlv_global
	    if( b3d ) write(iu,*) is,l,fluxesg_ob(l,:,is)
	  end do
	end do

	flush(iu)

	return
 1000	format(2i5,3g20.10)
	end

!******************************************************************

	subroutine finalize_nslayers

	use box
	use box_arrays
	use levels
	use basin

	implicit none

	integer is,nbs,ib1,ib2,nb1,nb2
	integer k,ie,i,ii,lmax
	integer kk(2)

!       isects(5,is)    description of section (n,type,box1,box2,ipnt)

	do is=1,nsect
	  nbs = nslayers(is)
	  ib1 = isects(3,is)
	  ib2 = isects(4,is)
	  nb1 = nblayers(ib1)
	  nb2 = nblayers(ib2)
	  !write(6,*) is,ib1,ib2,nbs,nb1,nb2
	  if( ib1 <= 0 ) nb1 = nb2	!do not use nb1
	  nbs = min(nbs,nb1,nb2)
	  nslayers(is) = nbs
	end do

	return

	write(6,*) 'ggg1 debug_nslayers: ',nblayers

	kk(1) = 8189
	kk(2) = 7545
	do ii=1,2
	k = kk(ii)
	do ie=1,nel
	  !i = findloc(nen3v(:,ie),k,1)		!compiler error
	  do i=1,3
	    if( nen3v(i,k) == k ) exit
	  end do
	  if( i > 3 ) cycle
	  if( i <= 0 ) cycle
	  lmax = ilhv(ie)
	  write(6,'(7i8)') k,nen3v(:,ie),ie,lmax,iboxes(ie)
	end do
	end do

	stop 'stop debug_nslayers'
	end

!******************************************************************

	subroutine write_section_info

!       isects(5,is)    description of section (n,type,box1,box2,ipnt)

	use box
	use box_arrays
	use levels
	use basin

	implicit none

	integer is,nbs

	write(6,*) 'total number of sections: ',nsect
	write(6,'(a)') '    is itype   ib1   ib2  lmax'

	do is=1,nsect
	  nbs = nslayers(is)
	  write(6,'(5i6)') is,isects(2:4,is),nbs
	end do

	!stop 'write_section_info'

	end

!******************************************************************

