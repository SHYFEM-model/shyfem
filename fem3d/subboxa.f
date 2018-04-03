c
c $Id: subboxa.f,v 1.25 2009-05-21 09:24:00 georg Exp $
c
c subroutines for computing discharge / flux for boxes
c
c contents :
c
c subroutine wrboxa(it)				write of flux data
c
c revision log :
c
c 30.04.1998    ggu	newly written routines (subpor deleted)
c 07.05.1998    ggu	check nrdveci on return for error
c 08.05.1998    ggu	restructured with new comodity routines
c 13.09.1999    ggu	type of node computed in own routine flxtype
c 19.11.1999    ggu	iskadj into sublin
c 20.01.2000	ggu	old routines substituted, new routine extrsect
c 20.01.2000    ggu     common block /dimdim/ eliminated
c 20.01.2000    ggu     common block /dimdim/ eliminated
c 01.09.2002	ggu	ggu99 -> bug in flx routines (how to reproduce?)
c 26.05.2003	ggu	in flxnov substituted a,b with b,c
c 26.05.2003	ggu	new routine make_fluxes (for lagrangian)
c 10.08.2003	ggu	do not call setweg, setnod, setkan
c 23.03.2006    ggu     changed time step to real
c 28.09.2007    ggu     use testbndo to determine boundary node in flxtype
c 28.04.2009    ggu     links re-structured
c 23.02.2011    ggu     new routine call write_node_fluxes() for special output
c 01.06.2011    ggu     documentation to flxscs() changed
c 21.09.2011    ggu     some lower-level subroutines copied to subflx.f
c 07.10.2011    ggu     adjusted for 3d flux routines
c 19.10.2011    ggu     added T/S variables, created fluxes_*() routines
c 19.10.2011    ggu     added conz variables, created fluxes_template()
c 10.05.2013    ggu     adapted to boxes computations
c 14.05.2013    ggu     write also OBC sections and contributions
c 29.10.2014    ggu     new code and 3d version
c 15.02.2018    ggu     code checked and error debug (is running)
c
c notes :
c
c insert parameters idtbox and itmbox into STR file to have box file written
c
c format of written files  please see box_write_*()
c
c still to check: some sections have more layers than adjacent boxes
c
c******************************************************************

!==================================================================
	module box
!==================================================================

        integer, save :: nbxdim = 0	!maximum for nbox
        integer, save :: nscboxdim = 0	!maximum for nsect
        integer, save :: nfxboxdim = 0	!maximum for kfluxm

        integer, save :: nbox		!total number of boxes
        integer, save :: nbc_ob		!total number of open boundaries
        integer, save :: nsect = -1	!total number of sections
        integer, save :: kfluxm = 0	!total number of nodes in all sections
        integer, save, allocatable :: kflux(:)	!node numbers defining sections

        integer, save, allocatable :: iflux(:,:)
        integer, save, allocatable :: iboxes(:)
        integer, save, allocatable :: ikboxes(:)
        integer, save, allocatable :: isects(:,:)
        integer, save, allocatable :: iscbnd(:,:)

!       iflux(3,i)      aux values for nodes in sections
!       iboxes(ie)      index from elements to boxes
!       ikboxes(k)      index from nodes to boxes (<0 if node on box boundary)
!       isects(4,is)    description of section (n,ipnt,box1,box2)
!       iscbnd(4,ibc)   description of OB (n,type,box1,box2)
!
!       i runs on list of nodes (1:kfluxm)
!       ie runs on elements (1:nel)
!       k runs on nodes (1:nkn)
!       is runs on sections (1:nsect)
!       ibc runs on open boundaries (1:nbc)
!       ib runs on boxes (1:nbox)

!==================================================================
	contains
!==================================================================

	subroutine box_alloc(nkn,nel,nlv,nbc,nb,ns,nn)

	integer nkn,nel,nlv,nbc,nb,ns,nn

	allocate(kflux(nn))

	allocate(iflux(3,nn))
	allocate(iboxes(nel))
	allocate(ikboxes(nkn))
	allocate(isects(4,ns))
	allocate(iscbnd(4,nbc))

	end subroutine box_alloc

!==================================================================
	end module box
!==================================================================

!==================================================================
	module box_arrays
!==================================================================

	implicit none

	real, save, allocatable :: fluxes(:,:,:)	!aux array for fluxes
	real, save, allocatable :: val2d(:)		!aux array for 2d vars
	real, save, allocatable :: val3d(:,:)		!aux array for 3d vars

	real, save, allocatable :: fluxes_m(:,:,:)	!mass fluxes

	integer, save, allocatable :: nslayers(:)	!layers in section
	double precision, save, allocatable :: masst(:,:,:)	!accum mass
	double precision, save, allocatable :: saltt(:,:,:)	!accum salt
	double precision, save, allocatable :: tempt(:,:,:)	!accum temp
	double precision, save, allocatable :: conzt(:,:,:)	!accum conz

	integer, save, allocatable :: nslayers_ob(:)	!layers in section
	real, save, allocatable :: fluxes_ob(:,:,:)	!discharges OBC (aux)
	double precision, save, allocatable :: masst_ob(:,:,:)	!OBC (dis-accum)

	integer, save, allocatable :: nblayers(:)	!layers in boxes
	real, save, allocatable :: barea(:)		!area of boxes
	real, save, allocatable :: bvolume(:)		!volume of boxes
	real, save, allocatable :: bdepth(:)		!depth of boxes

	real, save, allocatable :: valt(:,:)		!temperature of boxes
	real, save, allocatable :: vals(:,:)		!salinity of boxes
	real, save, allocatable :: valv(:,:)		!current speed of boxes
	real, save, allocatable :: vale(:)		!water level of boxes
	real, save, allocatable :: valm(:,:)		!meteo of boxes

	real, save, allocatable :: eta_act(:)		!new water level
	real, save, allocatable :: eta_old(:)		!old water level

	real, save, allocatable :: valw(:,:)		!vertical velocities
	real, save, allocatable :: valvis(:,:)		!viscosities
	real, save, allocatable :: valdif(:,:)		!diffusivities

!==================================================================
	contains
!==================================================================

	subroutine box_arrays_alloc(nlv,nbc,nb,ns,nn)

	integer nlv,nbc,nb,ns,nn

	allocate(fluxes(0:nlv,3,ns))
	allocate(val2d(nb))
	allocate(val3d(0:nlv,nb))

	allocate(fluxes_m(0:nlv,3,ns))

	allocate(nslayers(ns))
	allocate(masst(0:nlv,3,ns))
	allocate(saltt(0:nlv,3,ns))
	allocate(tempt(0:nlv,3,ns))
	allocate(conzt(0:nlv,3,ns))

	allocate(nslayers_ob(nbc))
	allocate(fluxes_ob(0:1,3,nbc))
	allocate(masst_ob(0:1,3,nbc))

	allocate(nblayers(nb))
	allocate(barea(nb))
	allocate(bvolume(nb))
	allocate(bdepth(nb))

	allocate(valt(0:nlv,nb))
	allocate(vals(0:nlv,nb))
	allocate(valv(0:nlv,nb))
	allocate(vale(nb))
	allocate(valm(7,nb))

	allocate(eta_old(nb))
	allocate(eta_act(nb))
	allocate(valw(0:nlv,nb))
	allocate(valvis(0:nlv,nb))
	allocate(valdif(0:nlv,nb))

	end subroutine box_arrays_alloc

!==================================================================
	end module box_arrays
!==================================================================

c******************************************************************
c******************************************************************
c******************************************************************
c general box routine (init/accum/aver/write)
c******************************************************************
c******************************************************************
c******************************************************************

	subroutine wrboxa

c administers writing of flux data

	use mod_conz, only : cnv
	use mod_ts
	use mod_diff_visc_fric
	use mod_hydro_vel
	use mod_hydro
	use basin
	use levels, only : nlvdi,nlv
	use box
	use box_arrays

	implicit none

	include 'simul.h'

	integer j,i,l,lmax,nlmax,ivar,nvers,nk_ob
	integer date,time
	integer idtbox
	integer nvar,ierr
	real az,azpar,dt
	double precision atime0,atime,dtime
	character*80 title,femver,aline

        double precision, save :: trm,trs,trt,trc
	double precision, save :: trob
	real, save :: dtbox			!need for time averaging

        integer, save :: nbbox = 0		!unit number for output
	integer, save :: ibarcl,iconz,ievap
	double precision, save :: da_out(4)

	integer nbnds,nkbnd
	integer ifemop
	integer ipext
	real getpar
	double precision dgetpar
	logical has_output_d,is_over_output_d,next_output_d

	logical bfirst
	logical bbox3d,b3d

	integer, allocatable :: kext(:)
	character*80, allocatable :: chflx(:)


c-----------------------------------------------------------------
c start of code
c-----------------------------------------------------------------

        if( nbbox .eq. -1 ) return

	bbox3d = .true.		!write 3d results if in 3d
	bbox3d = .false.	!write 3d results if in 3d

	b3d = nlv > 1

c-----------------------------------------------------------------
c initialization
c-----------------------------------------------------------------

        if( nbbox .eq. 0 ) then

          	call init_output_d('itmbox','idtbox',da_out)
		call increase_output_d(da_out)  !itbox=itmbox+idtbox
          	if( .not. has_output_d(da_out) ) nbbox = -1

                if( nbbox .eq. -1 ) return

		nbc_ob = nbnds()
		nk_ob = nkbnd()
		call box_check(nbox,nsect,kfluxm)
		nbxdim = nbox
		nscboxdim = nsect + nbc_ob
		nfxboxdim = kfluxm + nbc_ob + nk_ob
		call box_alloc(nkn,nel,nlvdi,nbc_ob,nbox,nscboxdim,nfxboxdim)
		call box_arrays_alloc(nlv,nbc_ob,nbox,nscboxdim,nfxboxdim)

		call box_init		!reads boxes.txt and sets up things
		call box_make_stats(nbox,iboxes,nblayers
     +					,barea,bvolume,bdepth)

                if( kfluxm .le. 0 ) nbbox = -1
                if( nsect .le. 0 ) nbbox = -1
                if( nbbox .eq. -1 ) return

                if( nsect .gt. nscboxdim ) then
                  stop 'error stop wrboxa: dimension nscboxdim'
                end if

		ibarcl = nint(getpar('ibarcl'))
		iconz = nint(getpar('iconz'))
		ievap = nint(getpar('ievap'))
		!ibarcl = 0
		!iconz = 0

		call get_nlayers(kfluxm,kflux,nslayers,nlmax)

		call fluxes_init_d(nlvdi,nsect,nslayers,trm,masst)
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
		call boxes_3d_init(nlvdi,nbox,valt)	!temp in box
		call boxes_3d_init(nlvdi,nbox,vals)	!salt in box
		call boxes_3d_init(nlvdi,nbox,valv)	!vel speed in box
		call boxes_2d_init(nbox,vale)		!zeta in box
		call boxes_multi_init(nbox,7,valm)	!meteo in box

		call boxes_3d_init(nlvdi,nbox,valw)	!vertical vel
		call boxes_3d_init(nlvdi,nbox,valvis)	!viscosity
		call boxes_3d_init(nlvdi,nbox,valdif)	!diffusivity

		call box_aver_eta(zeov,val2d)
		eta_old = val2d

		nslayers_ob = 1
		call fluxes_init_d(1,nbc_ob,nslayers_ob,trob,masst_ob)

                nbbox=ifemop('.box','unform','new')
                if(nbbox.le.0) then
        	   stop 'error stop wrboxa : Cannot open BOX file'
		end if
		da_out(4) = nbbox

	        nvers = 0
		idtbox = nint(da_out(1))
                call flx_write_header      (nbbox,nvers
     +                          ,nsect,kfluxm,idtbox,nlmax
     +                          ,nvar
     +                          ,ierr
     +                          )
		if( ierr /= 0 ) goto 98

                title = descrp
                call get_shyfem_version(femver)
                call get_absolute_ref_time(atime0)

		allocate(kext(kfluxm),chflx(nsect))
                do i=1,kfluxm
                  kext(i) = ipext(kflux(i))
                end do
		do i=1,nsect
		  write(chflx(i),'(a,i5)') 'internal box section ',i
		end do

                call flx_write_header2(nbbox,0,nsect,kfluxm
     +                          ,kext,nslayers
     +                          ,atime0,title,femver,chflx,ierr)
                if( ierr /= 0 ) goto 98

c               here we could also compute and write section in m**2

                date = nint(dgetpar('date'))
                time = nint(dgetpar('time'))
		call box_write_stats(date,time,idtbox
     +					,nbox,nsect,isects
     +					,kfluxm,kflux
     +					,nslayers,nblayers
     +					,barea,bvolume,bdepth)

        end if

c-----------------------------------------------------------------
c normal call
c-----------------------------------------------------------------

	call is_time_first(bfirst)
	if( bfirst ) then	!initial condition
	  call box_aver_eta(zenv,val2d)
	  eta_old = val2d
	  call box_write_init(nbox,eta_old)
	  return		!initial call to wrboxa
	end if

        if( .not. is_over_output_d(da_out) ) return	! it <= itmbox

	call get_timestep(dt)
	call getaz(azpar)
	az = azpar

c	-------------------------------------------------------
c	accumulate results
c	-------------------------------------------------------

	ivar = 0
	call flxscs(kfluxm,kflux,iflux,az,fluxes,ivar,rhov)
	call fluxes_accum_d(nlvdi,nsect,nslayers,dt,trm,masst,fluxes)

	call box_ob_compute(nbc_ob,fluxes_ob)	!open boundary conditions
	call fluxes_accum_d(1,nbc_ob,nslayers_ob,dt,trob
     +				,masst_ob,fluxes_ob)

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

	dtbox = dtbox + dt

	call box_3d_aver_scalar(tempv,val3d)
	call boxes_3d_accum(nlvdi,nbox,dt,valt,val3d)
	call box_3d_aver_scalar(saltv,val3d)
	call boxes_3d_accum(nlvdi,nbox,dt,vals,val3d)
	call box_3d_aver_vel(val3d)
	call boxes_3d_accum(nlvdi,nbox,dt,valv,val3d)
	call box_aver_eta(zenv,val2d)
	call boxes_2d_accum(nbox,dt,vale,val2d)
	eta_act = val2d

	call box_3d_aver_vertical(wlnv,val3d,0)
	call boxes_3d_accum(nlvdi,nbox,dt,valw,val3d)
	call box_3d_aver_vertical(visv,val3d,1)
	call boxes_3d_accum(nlvdi,nbox,dt,valvis,val3d)
	call box_3d_aver_vertical(difv,val3d,1)
	call boxes_3d_accum(nlvdi,nbox,dt,valdif,val3d)

	call boxes_meteo_accum(nbox,dt,valm,val2d,ievap)

c	-------------------------------------------------------
c	time for output?
c	-------------------------------------------------------

        if( .not. next_output_d(da_out) ) return

c	-------------------------------------------------------
c	time average results
c	-------------------------------------------------------

	nvers = 0
	call get_absolute_act_time(atime)
	call get_act_dtime(dtime)

	ivar = 0
	call fluxes_aver_d(nlvdi,nsect,nslayers,trm,masst,fluxes)
        call flx_write_record(nbbox,nvers,atime,nlvdi,nsect,ivar
     +                          ,nslayers,fluxes,ierr)
	if( ierr /= 0 ) goto 97

	fluxes_m = fluxes
	call fluxes_aver_d(1,nbc_ob,nslayers_ob,trob,masst_ob,fluxes_ob)

	if( ibarcl .gt. 0 ) then
	  ivar = 11
	  call fluxes_aver_d(nlvdi,nsect,nslayers,trs,saltt,fluxes)
          call flx_write_record(nbbox,nvers,atime,nlvdi,nsect,ivar
     +                          ,nslayers,fluxes,ierr)
	  if( ierr /= 0 ) goto 97
	  ivar = 12
	  call fluxes_aver_d(nlvdi,nsect,nslayers,trt,tempt,fluxes)
          call flx_write_record(nbbox,nvers,atime,nlvdi,nsect,ivar
     +                          ,nslayers,fluxes,ierr)
	  if( ierr /= 0 ) goto 97
	end if

	if( iconz .eq. 1 ) then
	  ivar = 10
	  call fluxes_aver_d(nlvdi,nsect,nslayers,trc,conzt,fluxes)
          call flx_write_record(nbbox,nvers,atime,nlvdi,nsect,ivar
     +                          ,nslayers,fluxes,ierr)
	  if( ierr /= 0 ) goto 97
	end if

	call boxes_3d_aver(nlvdi,nbox,dtbox,valt)	!temp in box
	call boxes_3d_aver(nlvdi,nbox,dtbox,vals)	!salt in box
	call boxes_3d_aver(nlvdi,nbox,dtbox,valv)	!vel speed in box
	call boxes_2d_aver(nbox,dtbox,vale)		!zeta in box

	call boxes_3d_aver(nlvdi,nbox,dtbox,valw)
	call boxes_3d_aver(nlvdi,nbox,dtbox,valvis)
	call boxes_3d_aver(nlvdi,nbox,dtbox,valdif)

	call boxes_multi_aver(nbox,7,dtbox,valm)	!meteo in box

	val2d = barea*(eta_act-eta_old)
	call boxes_mass_balance_2d(dtbox,bvolume,fluxes_m,fluxes_ob,val2d)
	if( b3d  ) then
	  call boxes_mass_balance_3d(dtbox,bvolume
     +					,nblayers,nslayers
     +					,fluxes_m,fluxes_ob,val2d,valw)
	end if

c	-------------------------------------------------------
c	write results
c	-------------------------------------------------------

	call get_act_timeline(aline)
	write(6,*) 'writing box results... ',aline

	call box_write_2d(dtime,aline,fluxes_m,valt,vals,eta_act,valv
     +					,fluxes_ob)
	call box_write_meteo(dtime,aline,valm)

	if( bbox3d .and. b3d ) then
	  call box_write_3d(dtime,aline,nblayers,nslayers
     +			,fluxes,valt,vals,vale,valv,fluxes_ob)
	  call box_write_vertical(dtime,aline,nblayers
     +			,valw,valvis,valdif)
	end if

c	-------------------------------------------------------
c	reset variables
c	-------------------------------------------------------

	call fluxes_init_d(nlvdi,nsect,nslayers,trm,masst)
	call fluxes_init_d(1,nbc_ob,nslayers_ob,trob,masst_ob)

	if( ibarcl .gt. 0 ) then
	  call fluxes_init_d(nlvdi,nsect,nslayers,trs,saltt)
	  call fluxes_init_d(nlvdi,nsect,nslayers,trt,tempt)
	end if

	if( iconz .eq. 1 ) then
	  call fluxes_init_d(nlvdi,nsect,nslayers,trc,conzt)
	end if

	dtbox = 0.
	call boxes_3d_init(nlvdi,nbox,valt)	!temp in box
	call boxes_3d_init(nlvdi,nbox,vals)	!salt in box
	call boxes_2d_init(nbox,vale)		!zeta in box
	call boxes_3d_init(nlvdi,nbox,valv)	!vel speed in box

	call boxes_multi_init(nbox,7,valm)	!meteo

	call boxes_3d_init(nlvdi,nbox,valw)	!vertical vel
	call boxes_3d_init(nlvdi,nbox,valvis)	!viscosity
	call boxes_3d_init(nlvdi,nbox,valdif)	!diffusivity

	eta_old = eta_act

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	return
   97	continue
	write(6,*) 'error writing record of boxes flx file: ',ivar
	stop 'error stop wrboxa: record'
   98	continue
	write(6,*) 'error writing header of boxes flx file'
	stop 'error stop wrboxa: header'
	end

c******************************************************************
c******************************************************************
c******************************************************************
c reading and setup routines
c******************************************************************
c******************************************************************
c******************************************************************

	subroutine box_init

c reads file etc...

	use box

	implicit none

	logical berror
	integer nbc
	integer nbnds
	integer nsaux
	integer i

	call box_read				!reads boxes.txt
	call box_elab(iboxes,ikboxes)		!sets up ikboxes

	call n2int(kfluxm,kflux,berror)		!converts ext to int nodes
        if( berror ) stop 'error stop box_init: converting node numbers'

	call box_ob_init	!initializes OB contributions (internal nodes)

	iflux = 0
	call flx_init(kfluxm,kflux,nsaux,iflux)
	if( nsaux /= nsect ) then
	  write(6,*) 'sections: ',nsaux,nsect
	  stop 'error stop box_init: different number of sections read'
	end if

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

c******************************************************************

	subroutine box_check(nbox,nsect,kfluxm)

! reads boxes.txt just to return dimensions

	implicit none

	integer nbox,nsect,kfluxm

	integer nb,id,n,ib1,ib2,i
	integer nel,ie
	integer iaux

	nbox = 0
	nsect = 0
	kfluxm = 0

	open(1,file='boxes.txt',form='formatted',status='old',err=94)

	read(1,*) nel,nbox,nb
	read(1,*) (iaux,ie=1,nel)

	do
	  read(1,*) id,n,ib1,ib2
	  if( id .le. 0 ) exit
	  read(1,*) (iaux,i=1,n)
	  nsect = nsect + 1
	  if( id .ne. nsect ) goto 97
	  kfluxm = kfluxm + n + 1
	end do

	close(1)

	return
   94	continue
	write(6,*) 'file: boxes.txt'
	stop 'error stop box_check: cannot read file'
   97	continue
	write(6,*) id,nsect
	stop 'error stop box_check: id'
	end

c******************************************************************

	subroutine box_read

! reads boxes.txt and returns parameters and arrays
!
! on return the following values are set:
!
! nbox,nsect,kfluxm
! kflux

	use basin
	use box

	implicit none

	integer nb,id,n,ib1,ib2,i
	integer nelaux,ie
	integer ierr,iaux
	integer iauxv(nkn)

	nsect = 0
	kfluxm = 0
	kflux = 0
	ierr = 0

	write(6,*) 'start reading boxes.txt... '
	open(1,file='boxes.txt',form='formatted',status='old',err=94)

	read(1,*) nelaux,nbox,nb
	if( nelaux .ne. nel ) goto 99
	if( nbox .gt. nbxdim ) goto 98
	read(1,*) (iboxes(ie),ie=1,nel)

	do
	  read(1,*) id,n,ib1,ib2
	  if( id .le. 0 ) exit
	  read(1,*) (iauxv(i),i=1,n)
	  nsect = nsect + 1
	  if( id .ne. nsect ) goto 97
	  if( nsect .gt. nscboxdim ) goto 95
	  call box_insert_section(id,n,ib1,ib2,iauxv,ierr)
	  if( kfluxm .gt. nfxboxdim ) goto 96
	  if( ierr .gt. 0 ) goto 93		!should be already handled above
	end do

	close(1)

	kfluxm = kfluxm - 1

	write(6,*) 'finished reading boxes.txt: '
	write(6,*) nbox,nsect,kfluxm,nscboxdim,nfxboxdim

	return
   93	continue
	write(6,*) ierr
	stop 'error stop box_read: unspecified error'
   94	continue
	write(6,*) 'file: boxes.txt'
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
	stop 'error stop box_read: nbxdim'
   99	continue
	write(6,*) nel,nelaux
	stop 'error stop box_read: nel'
	end

c******************************************************************

	subroutine box_insert_section(id,n,ib1,ib2,list,ierr)

c inserts section in data structure

	use basin
	use box

	implicit none

	integer id,n,ib1,ib2
	integer list(n)
	integer ierr

	integer i

	if( id .gt. nscboxdim ) ierr = ierr + 1
	if( kfluxm + n + 1 .gt. nfxboxdim ) ierr = ierr + 1

	if( ierr .eq. 0 ) then		!no errors -> insert
	    isects(1,id) = n		!total number of nodes
	    isects(2,id) = kfluxm+1	!pointer into kflux
	    isects(3,id) = ib1		!from box
	    isects(4,id) = ib2		!to box
	    do i=1,n
	      kflux(kfluxm+i) = list(i)
	    end do
	    kfluxm = kfluxm + n + 1
	    kflux(kfluxm) = 0
	else
	    kfluxm = kfluxm + n + 1
	end if

	end

c******************************************************************

	subroutine box_elab(iboxes,ikboxes)

c sets up ikboxes - index of nodes to boxes
c
c is -1 if node belongs to more than one box

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

	if( nb+nl .ne. nkn ) stop 'error stop box_elab: internal error'

	write(6,*) 'box nodes: ',nb,nl,nb+nl,nkn

	end

c******************************************************************

	subroutine box_make_stats(nbox,iboxes,nblayers
     +					,barea,bvolume,bdepth)

c computes max layers for each box

	use mod_depth
	use evgeom
	use basin
	use levels

	implicit none

	integer nbox
	integer iboxes(nel)		!number of layers in boxes
	integer nblayers(nbox)		!number of layers in boxes
	real barea(nbox)		!area of boxes
	real bvolume(nbox)		!volume of boxes
	real bdepth(nbox)		!depth of boxes

	integer ib,lmax,ie
	real area,hdep

	nblayers = 0
	barea = 0.
	bvolume = 0.
	bdepth = 0.

	do ie=1,nel
	  ib = iboxes(ie)
	  hdep = hev(ie)
	  area = 12.*ev(10,ie)
	  barea(ib) = barea(ib) + area
	  bvolume(ib) = bvolume(ib) + area * hdep
	  lmax = ilhv(ie)
	  nblayers(ib) = max(nblayers(ib),lmax)
	end do

	bdepth = bvolume / barea

	end

c******************************************************************
c******************************************************************
c******************************************************************
c writing routines
c******************************************************************
c******************************************************************
c******************************************************************

	subroutine box_write_stats(date,time,idtbox
     +					,nbox,nsect,isects
     +					,kfluxm,kflux
     +					,nslayers,nblayers
     +					,barea,bvolume,bdepth)

c writes statistics to file

	use basin
	use mod_depth
	use evgeom

	implicit none

	integer date,time
	integer idtbox
	integer nbox,nsect
	integer i,ipt,k,ke,n,ndim
	integer isects(4,nsect)
	integer kfluxm,kflux(kfluxm)
	integer nslayers(nsect)		!number of layers in section
	integer nblayers(nbox)		!number of layers in boxes
	real barea(nbox)		!area of boxes
	real bvolume(nbox)		!volume of boxes
	real bdepth(nbox)		!depth of boxes

	integer ib,iu
	integer is,ib1,ib2
	integer nb1,nb2
	integer ie
	real area,depth
	double precision areatot
	double precision dtime,atime,atime0
	character*20 line

	integer nsbox(-1:nbox)			!number of nodes in box section
	integer, allocatable :: kbox(:,:)	!nodes of section of box

	integer ifileo,ipext
	character*80 file

	file = 'boxes_stats.txt'
	iu = ifileo(0,file,'formatted','new')
	if( iu <= 0 ) stop 'error stop boxes: opening file'

	write(iu,*) date,time
	write(iu,*) idtbox
	call get_absolute_ref_time(atime0)
	call get_first_dtime(dtime)
	atime = atime0 + dtime
        call dts_format_abs_time(atime,line)
	write(iu,*) 'start of simulation = ',line
	call get_last_dtime(dtime)
	atime = atime0 + dtime
        call dts_format_abs_time(atime,line)
	write(iu,*) 'end of simulation = ',line

	write(iu,*) nbox
	do ib=1,nbox
	  write(iu,*) ib,nblayers(ib),barea(ib),bvolume(ib),bdepth(ib)
	end do

	nsbox = 0

        write(iu,*) nsect
        do is=1,nsect
          ib1 = isects(3,is)
	  nb1 = 0
	  if( ib1 > 0 ) nb1 = nblayers(ib1)
          ib2 = isects(4,is)
	  nb2 = 0
	  if( ib2 > 0 ) nb2 = nblayers(ib2)
          write(iu,*) ib1,ib2,nslayers(is),nb1,nb2
          n = isects(1,is)
	  if( ib1 > 0 ) nsbox(ib1) = nsbox(nb1) + n + 1
	  if( ib2 > 0 ) nsbox(ib2) = nsbox(nb2) + n + 1
        end do

	ndim = maxval(nsbox)
	allocate(kbox(0:ndim,nbox))
	kbox = 0

        write(iu,*) nsect
        do is=1,nsect
          n = isects(1,is)
          ipt = isects(2,is)
          ib1 = isects(3,is)
          ib2 = isects(4,is)
          write(iu,*) is,n,ib1,ib2,nslayers(is)
	  do i=1,n
	    k = kflux(ipt-1+i)
	    ke = ipext(k)
	    write(iu,*) i,ke,xgv(k),ygv(k)
	    call insert_sect_node(ndim,nbox,ib1,kbox,k)
	    call insert_sect_node(ndim,nbox,ib2,kbox,k)
	  end do
	  call insert_sect_node(ndim,nbox,ib1,kbox,0)
	  call insert_sect_node(ndim,nbox,ib2,kbox,0)
	end do

	write(iu,*) nbox
	do ib=1,nbox
	  n = kbox(0,ib)
	  write(iu,*) ib,n
	  do i=1,n
	    k = kbox(i,ib)
	    ke = ipext(k)
	    if( k == 0 ) then
	      write(iu,*) i,0,0.,0.
	    else
	      write(iu,*) i,ke,xgv(k),ygv(k)
	    end if
	  end do
	end do

	close(iu)

	file = 'boxes_geom.txt'
	iu = ifileo(0,file,'formatted','new')
	if( iu <= 0 ) stop 'error stop boxes: opening file'

	areatot = 0.
	write(iu,*) nel
	do ie=1,nel
	  area = 12. * ev(10,ie)
	  areatot = areatot + area
	  depth = hev(ie)
	  write(iu,*) ie,area,depth
	end do
	write(6,*) 'total area: ',areatot

	close(iu)

	end

c******************************************************************

	subroutine insert_sect_node(ndim,nbox,ib,kbox,k)

	implicit none

	integer ndim,nbox,ib,k
	integer kbox(0:ndim,nbox)

	integer i

	if( ib <= 0 ) return

	i = kbox(0,ib)
	i = i + 1
	kbox(i,ib) = k
	kbox(0,ib) = i

	end

c******************************************************************

	subroutine box_write_init(nbox,eta_act)

c writes initial conditions for eta - still to be done : T/S

	implicit none

	integer nbox
	real eta_act(nbox)			!water level of first time step

	integer ib,iu
	character*80 file

	integer ifileo

	file = 'boxes_init.txt'
	iu = ifileo(135,file,'formatted','new')
	if( iu <= 0 ) stop 'error stop boxes: opening file'

	write(iu,*) nbox
	do ib=1,nbox
	  write(iu,*) ib,eta_act(ib)
	end do

	close(iu)

	end

c******************************************************************

	subroutine box_write_2d(dtime,aline,fluxes
     +				,valt,vals,vale,valv,fluxes_ob)

c writes 2d vertical average box values to file
c
c for boxes variables are:
c	1	temperature
c	2	salinity
c	3	water level
c	4	current velocity

	use basin
	use levels
	use box

	implicit none

	double precision dtime
	character*(*) aline
	real fluxes(0:nlvdi,3,nsect)
	real valt(0:nlvdi,nbox)			!temperature
	real vals(0:nlvdi,nbox)			!salinity
	real vale(nbox)				!zeta
	real valv(0:nlvdi,nbox)			!velocity speed
	real fluxes_ob(0:1,3,nbc_ob)

	integer ib,is,ib1,ib2,ii,itype

	integer ifileo
	character*80 file
	integer, save :: iu = 0

	if( iu == 0 ) then
	  file = 'boxes_2d.txt'
	  iu = ifileo(0,file,'formatted','new')
	  if( iu <= 0 ) stop 'error stop boxes: opening file'
	end if

	write(iu,*) dtime,trim(aline)

	write(iu,*) nbox
	do ib=1,nbox
	  write(iu,*) ib,valt(0,ib),vals(0,ib),vale(ib),valv(0,ib)
	end do

	write(iu,*) nsect
	do is=1,nsect
	  ib1 = isects(3,is)
	  ib2 = isects(4,is)
	  write(iu,*) ib1,ib2,(fluxes(0,ii,is),ii=1,3)
	end do

	write(iu,*) nbc_ob
	do is=1,nbc_ob
	  itype = iscbnd(2,is)
	  ib1 = iscbnd(3,is)
	  ib2 = iscbnd(4,is)
	  if( itype .eq. 1 ) then	!z-boundary, ignore
	    ib1 = 0
	    ib2 = 0
	  end if
	  write(iu,*) ib1,ib2,(fluxes_ob(0,ii,is),ii=1,3)
	end do

	end

c******************************************************************

	subroutine box_write_3d(dtime,aline,nblayers,nslayers
     +			,fluxes,valt,vals,vale,valv,fluxes_ob)

c writes 3d box values to file

	use basin
	use levels
	use box

	implicit none

	double precision dtime
	character*(*) aline
	integer nblayers(nbox)		!number of layers in box
	integer nslayers(nsect)		!number of layers in section
	real fluxes(0:nlvdi,3,nsect)
	real valt(0:nlvdi,nbox)			!temperature
	real vals(0:nlvdi,nbox)			!salinity
	real vale(nbox)				!zeta
	real valv(0:nlvdi,nbox)			!velocity speed
	real fluxes_ob(0:1,3,nbc_ob)

	integer ib,is,ib1,ib2,ii,itype
	integer lmax,l

	integer ifileo
	character*80 file
	integer, save :: iu = 0

	if( iu == 0 ) then
	  file = 'boxes_3d.txt'
	  iu = ifileo(0,file,'formatted','new')
	  if( iu <= 0 ) stop 'error stop boxes: opening file'
	end if

	write(iu,*) dtime,trim(aline)

	write(iu,*) nbox
	do ib=1,nbox
	  lmax = nblayers(ib)
	  write(iu,*) ib,lmax,vale(ib)
	  do l=1,lmax
	    write(iu,*) l,valt(l,ib),vals(l,ib),valv(l,ib)
	  end do
	end do

	write(iu,*) nsect
	do is=1,nsect
	  lmax = nslayers(is)
	  ib1 = isects(3,is)
	  ib2 = isects(4,is)
	  write(iu,*) ib1,ib2,lmax
	  do l=1,lmax
	    write(iu,*) l,(fluxes(l,ii,is),ii=1,3)
	  end do
	end do

	write(iu,*) nbc_ob
	do is=1,nbc_ob
	  itype = iscbnd(2,is)
	  ib1 = iscbnd(3,is)
	  ib2 = iscbnd(4,is)
	  if( itype .eq. 1 ) then	!z-boundary, ignore
	    ib1 = 0
	    ib2 = 0
	  end if
	  write(iu,*) ib1,ib2,1
	  write(iu,*) 0,(fluxes_ob(0,ii,is),ii=1,3)
	end do

	end

c******************************************************************

	subroutine box_write_vertical(dtime,aline,nblayers
     +			,valw,valvis,valdif)

c writes 3d interface box values to file

	use basin
	use levels
	use box

	implicit none

	double precision dtime
	character*(*) aline
	integer nblayers(nbox)		!number of layers in box
	real valw(0:nlvdi,nbox)		!vertical velocities
	real valvis(0:nlvdi,nbox)		!viscosity
	real valdif(0:nlvdi,nbox)		!diffusivity

	integer ib,is,ib1,ib2,ii,itype
	integer lmax,l,max

	integer ifileo
	character*80 file
	integer, save :: iu = 0

	if( iu == 0 ) then
	  file = 'boxes_vertical.txt'
	  iu = ifileo(0,file,'formatted','new')
	  if( iu <= 0 ) stop 'error stop boxes: opening file'
	end if

	write(iu,*) dtime,trim(aline)

	write(iu,*) nbox
	do ib=1,nbox
	  lmax = nblayers(ib)
	  write(iu,*) ib,lmax
	  do l=1,lmax
	    write(iu,*) l,valw(l,ib),valvis(l,ib),valdif(l,ib)
	  end do
	end do

	end

c******************************************************************
c******************************************************************
c******************************************************************
c box average routines
c******************************************************************
c******************************************************************
c******************************************************************

	subroutine box_aver_eta(zev,val)

c computes average zeta values for box

	use mod_geom_dynamic
	use evgeom
	use basin
	use box

	implicit none

        real zev(3,nel)
	real val(nbox)

	double precision vald(nbox)
	double precision vold(nbox)

	integer ib,ie,ii
	real vol,area3,z

	vald = 0.
	vold = 0.

	do ie=1,nel
	  !if( iwegv(ie) .eq. 0 ) then	!only for wet elements
	    ib = iboxes(ie)
	    area3 = 4.*ev(10,ie)
	    do ii=1,3
	      z = zev(ii,ie)
	      vol = area3
	      vald(ib) = vald(ib) + vol*z
	      vold(ib) = vold(ib) + vol
	    end do
	  !end if
	end do

	val = 0.
	where( vold > 0. ) val = vald / vold

	end

c******************************************************************

	subroutine box_3d_aver_vel(val)

c computes average velocity values (speed) for box

	use mod_layer_thickness
	use mod_hydro
	use evgeom
	use levels
	use basin
	use box

	implicit none

	real val(0:nlvdi,nbox)

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

	do ie=1,nel
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

c******************************************************************

	subroutine box_3d_aver_vertical(scalar,val,iaver)

c computes average scalar values defined on interface (verticals) for box

	use mod_depth
	use evgeom
	use levels
	use basin
	use box

	implicit none

	real scalar(0:nlvdi,nkn)
	real val(0:nlvdi,nbox)
	integer iaver			!average or only accumulate

	double precision vald(0:nlvdi,nbox)
	double precision vold(0:nlvdi,nbox)

	logical baver
	integer ib,ie,lmax,l,ii,k
	real vol,area3

	baver = iaver > 0

	vald = 0.
	vold = 0.

	do ie=1,nel
	  ib = iboxes(ie)
	  area3 = 4.*ev(10,ie)
	  lmax = ilhv(ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    do l=0,lmax
	      vol = area3
	      vald(l,ib) = vald(l,ib) + vol*scalar(l,k)
	      vold(l,ib) = vold(l,ib) + vol
	    end do
	  end do
	end do

	if( .not. baver ) return	!only accumulate

	val = 0.
	where( vold > 0. ) val = vald / vold

	end

c******************************************************************

	subroutine box_3d_aver_scalar(scalar,val)

c computes average scalar values for box

	use mod_layer_thickness
	use evgeom
	use levels
	use basin
	use box

	implicit none

	real scalar(nlvdi,nkn)
	real val(0:nlvdi,nbox)

	double precision vald(0:nlvdi,nbox)
	double precision vold(0:nlvdi,nbox)

	integer ib,ie,lmax,l,ii,k
	real vol,area3

	vald = 0.
	vold = 0.

	do ie=1,nel
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

	val = 0.
	where( vold > 0. ) val = vald / vold

	end

c******************************************************************

	subroutine box_2d_aver_scalar(scalar,val)

c computes average 2d scalar values for box

	use evgeom
	use basin
	use box

	implicit none

	real scalar(nkn)
	real val(nbox)

	double precision vald(nbox)
	double precision vold(nbox)

	integer ib,ie,ii,k
	real vol,area3

	vald = 0.
	vold = 0.

	do ie=1,nel
	  ib = iboxes(ie)
	  area3 = 4.*ev(10,ie)
	  vol = area3
	  do ii=1,3
	    k = nen3v(ii,ie)
	    vald(ib) = vald(ib) + vol*scalar(k)
	    vold(ib) = vold(ib) + vol
	  end do
	end do

	val = 0.
	where( vold > 0. ) val = vald / vold

	end

c******************************************************************
c******************************************************************
c******************************************************************
c time init/accum/aver routines for boxes
c******************************************************************
c******************************************************************
c******************************************************************

	subroutine boxes_2d_init(n,val)

	implicit none

	integer n
	real val(n)

	val = 0.

	end

c******************************************************************

	subroutine boxes_2d_accum(n,dt,val,value)

	implicit none

	integer n
	real dt
	real val(n)
	real value(n)

	val = val + value * dt

	end

c******************************************************************

	subroutine boxes_2d_aver(n,dt,val)

	implicit none

	integer n
	real dt
	real val(n)

	if( dt .le. 0. ) return
	val = val / dt

	end

c******************************************************************

	subroutine boxes_3d_init(nlvddi,n,val)

	implicit none

	integer nlvddi,n
	real val(0:nlvddi,n)

	val = 0.

	end

c******************************************************************

	subroutine boxes_3d_accum(nlvddi,n,dt,val,value)

	implicit none

	integer nlvddi,n
	real dt
	real val(0:nlvddi,n)
	real value(0:nlvddi,n)

	val = val + value * dt

	end

c******************************************************************

	subroutine boxes_3d_aver(nlvddi,n,dt,val)

	implicit none

	integer nlvddi,n
	real dt
	real val(0:nlvddi,n)

	if( dt .le. 0. ) return
	val = val / dt

	end

c******************************************************************
c******************************************************************
c******************************************************************
c open boundary handling routines
c******************************************************************
c******************************************************************
c******************************************************************

	subroutine box_ob_init

c sets up open boundary inputs

	use basin
	use box

	implicit none

	integer nbc,ibc,itype
	integer nk,i,k,ib,ibk
	integer ierr
	integer iauxv(nkn)

	integer nbnds,itybnd,nkbnds,kbnds,ndim
	integer ipext

	ierr = 0
	kfluxm = kfluxm + 1
	kflux(kfluxm) = 0
	ndim = nkn

	nbc = nbnds()

	do ibc = 1,nbc
	  !write(6,*) ibc,nbc
	  itype = itybnd(ibc)
	  nk = nkbnds(ibc)
	  ib = 0
	  do i=1,nk
	    k = kbnds(ibc,i)
	    ibk = ikboxes(k)
	    if( ibk .le. 0 ) goto 99	!node not in unique box
	    if( ib .eq. 0 ) ib = ibk
	    if( ibk .ne. ib ) goto 98	!boundary nodes in different boxes
	  end do
	  iscbnd(1,ibc) = nk		!total number of boundary nodes
	  iscbnd(2,ibc) = itype		!type of boundary
	  iscbnd(3,ibc) = -ibc		!from box
	  iscbnd(4,ibc) = ib		!to box

	  if( itype .eq. 1 ) then	!insert sections of z boundary
	    nsect = nsect + 1
	    if( nsect .gt. nscboxdim ) goto 95
	    call irbnds(ibc,ndim,nk,iauxv)
	    call box_insert_section(nsect,nk,-ibc,ib,iauxv,ierr)
	    if( kfluxm .gt. nfxboxdim ) goto 96
	    if( ierr /= 0 ) goto 93
	  end if
	end do

	kfluxm = kfluxm - 1

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

c******************************************************************

	subroutine box_ob_compute(nbc_ob,fluxes_ob)

c computes open boundary inputs
c
c these fluxes are barotropic, so we insert them in 0 (baro) and 1 (1 layer)

	implicit none

	integer nbc_ob
	real fluxes_ob(0:1,3,nbc_ob)		!discharges at open bound

	integer ibc,nbc
	real val

	real get_discharge

	nbc = nbc_ob

	do ibc = 1,nbc
	  val = get_discharge(ibc)
	  if( val .ge. 0 ) then
	    fluxes_ob(0,1,ibc) = val
	    fluxes_ob(0,2,ibc) = val
	    fluxes_ob(0,3,ibc) = 0.
	    fluxes_ob(1,1,ibc) = val
	    fluxes_ob(1,2,ibc) = val
	    fluxes_ob(1,3,ibc) = 0.
	  else
	    fluxes_ob(0,1,ibc) = val
	    fluxes_ob(0,2,ibc) = 0.
	    fluxes_ob(0,3,ibc) = -val
	    fluxes_ob(1,1,ibc) = val
	    fluxes_ob(1,2,ibc) = 0.
	    fluxes_ob(1,3,ibc) = -val
	  end if
	end do

	end

c******************************************************************
c******************************************************************
c******************************************************************
c meteo routines
c******************************************************************
c******************************************************************
c******************************************************************

        subroutine boxes_multi_init(n,nj,val)

c as boxes_init, but for multiple values (only 2d)

        implicit none

        integer n,nj
        real val(nj,n)

	val = 0.

        end

c******************************************************************

	subroutine boxes_multi_aver(n,nj,dt,val)

c as boxes_aver, but for multiple values (only 2d)

	implicit none

	integer n,nj
	real dt
	real val(nj,n)

	if( dt .le. 0. ) return
	val = val / dt

	end

c******************************************************************

	subroutine boxes_meteo_accum(nbox,dt,valm,val,ievap)

! meteo parameter accumulation
!
! rain and evaporation are in [m/s]
! before writing they are converted to [mm/day]

	use mod_meteo

	implicit none

	integer nbox
	real dt
	real valm(7,nbox)
	real val(nbox)
	integer ievap

        real zconv
        parameter( zconv = 1.e-3 / 86400. )	!convert from [mm/day] to [m/s]

	integer i,ip
	real econv,rconv

	rconv = 1. / zconv
	econv = 1. / zconv
	if( ievap .le. 0 ) econv = 0.

	ip = 1
	call box_2d_aver_scalar(metrad,val)
	valm(ip,:) = valm(ip,:) + dt * val(:)

	ip = 2
	call box_2d_aver_scalar(mettair,val)
	valm(ip,:) = valm(ip,:) + dt * val(:)

	ip = 3
	call box_2d_aver_scalar(methum,val)
	valm(ip,:) = valm(ip,:) + dt * val(:)

	ip = 4
	call box_2d_aver_scalar(metws,val)
	valm(ip,:) = valm(ip,:) + dt * val(:)

	ip = 5
	call box_2d_aver_scalar(metcc,val)
	valm(ip,:) = valm(ip,:) + dt * val(:)

	ip = 6
	call box_2d_aver_scalar(metrain,val)
	valm(ip,:) = valm(ip,:) + dt * val(:)*rconv		![mm/day]

	ip = 7
	call box_2d_aver_scalar(evapv,val)
	valm(ip,:) = valm(ip,:) + dt * val(:)*econv		![mm/day]

	end

c******************************************************************

	subroutine box_write_meteo(dtime,aline,valm)

! parameters are:
!	1	srad
!	2	tair
!	3	rhum
!	4	wspeed
!	5	cc
!	6	rain
!	7	evap

	use basin
	use box

	implicit none

	double precision dtime
	character*(*) aline
	real valm(7,nbox)			!meteo variables

	integer ib,j

	integer ifileo
	character*80 file,dline
	integer, save :: iu = 0

	if( iu == 0 ) then
	  file = 'boxes_meteo.txt'
	  iu = ifileo(0,file,'formatted','new')
	  if( iu <= 0 ) stop 'error stop boxes: opening file'
	end if

	write(iu,*) dtime,trim(aline)

	write(iu,*) nbox
	do ib=1,nbox
	  write(iu,1001) ib,(valm(j,ib),j=1,7)
	end do

 1000	format(i14,7g12.4)
 1001	format(i4,1x,5f8.3,2g12.4)
	end

c******************************************************************
c******************************************************************
c******************************************************************
c mass balance routines
c******************************************************************
c******************************************************************
c******************************************************************

	subroutine boxes_mass_balance_2d(dtbox,bvolume
     +					,fluxes,fluxes_ob,vale)

c computes mass balance

	use basin
	use levels
	use box

	implicit none

	real dtbox
	real bvolume(nbox)
	real fluxes(0:nlvdi,3,nsect)
	real fluxes_ob(0:1,3,nbc_ob)
	real vale(nbox)		!volume difference due to eta

	integer ib,is,ib1,ib2
	real volf,vole,volb,vola,vold,volr
	real errmax
	double precision vol(nbox)
	double precision flux
	logical bdebug

	real, parameter :: eps = 1.e-3

	bdebug = .true.
	bdebug = .false.

	vol = 0.

	do is=1,nsect
	  ib1 = isects(3,is)
	  ib2 = isects(4,is)
	  flux = dtbox * fluxes(0,1,is)
	  if( ib1 > 0 ) vol(ib1) = vol(ib1) - flux
	  if( ib2 > 0 ) vol(ib2) = vol(ib2) + flux
	end do

	do is=1,nbc_ob
	  ib1 = iscbnd(3,is)
	  ib2 = iscbnd(4,is)
	  flux = dtbox * fluxes_ob(0,1,is)
	  if( ib1 > 0 ) vol(ib1) = vol(ib1) - flux
	  if( ib2 > 0 ) vol(ib2) = vol(ib2) + flux
	end do

	errmax = 0.
	do ib=1,nbox
	  volf = vol(ib)
	  vole = vale(ib)
	  volb = bvolume(ib)
	  vola = max(abs(volf),abs(vole))
	  vold = abs(volf-vole)
	  volr = 0.
	  !volr = vold/vola
	  volr = vold/volb	!use total volume
	  if( volr .gt. eps ) then
	    write(6,*) '*** mass balance error in box ',ib
	    write(6,*) '***',ib,volf,vole,vold,volr
	  end if
	  errmax = max(errmax,abs(volr))
	  !write(115,*) ib,volf,vole,vold,volr
	end do

	if( bdebug ) write(6,*) 'errmax 2d: ',errmax
	if( errmax .gt. eps ) then
	  write(6,*) 'errmax 2d: ',errmax
	end if

	end

c******************************************************************

	subroutine boxes_mass_balance_3d(dtbox,bvolume
     +					,nblayers,nslayers
     +					,fluxes,fluxes_ob
     +					,vale,wflux)

c computes mass balance

	use basin
	use levels
	use box

	implicit none

	real dtbox
	real bvolume(nbox)
	integer nblayers(nbox)
	integer nslayers(nsect)		!number of layers in section
	real fluxes(0:nlvdi,3,nsect)
	real fluxes_ob(0:1,3,nbc_ob)
	real vale(nbox)			!volume difference due to eta
	real wflux(0:nlvdi,nbox)		!vertical fluxes [m**3/s]

	logical bdebug
	integer ib,is,ib1,ib2
	integer l,lmax,ltot
	real volf,vole,volb,vola,vold,volr,volv
	real wtop,errmax
	double precision vol(nlvdi,nbox)
	double precision flux

	real, parameter :: eps = 1.e-4

	bdebug = .true.
	bdebug = .false.

	ltot = 0
	do ib=1,nbox
	  ltot = max(ltot,nblayers(ib))
	end do

	vol = 0.

	do is=1,nsect
	  ib1 = isects(3,is)
	  ib2 = isects(4,is)
	  lmax = nslayers(is)
	  do l=1,lmax
	    flux = dtbox * fluxes(l,1,is)
	    if( ib1 > 0 ) vol(l,ib1) = vol(l,ib1) - flux
	    if( ib2 > 0 ) vol(l,ib2) = vol(l,ib2) + flux
	    if( bdebug .and. abs(flux) > 0. ) then
	      if( nblayers(ib1) < l .or. nblayers(ib2) < l ) then
		write(6,'(a,6i5,g14.4)') 'non existing flux: '
     +				,is,ib1,ib2
     +				,nblayers(ib1),nblayers(ib2),l,flux
	      end if
	    end if
	  end do
	end do

	do is=1,nbc_ob
	  ib1 = iscbnd(3,is)
	  ib2 = iscbnd(4,is)
	  flux = dtbox * fluxes_ob(1,1,is)
	  if( ib1 > 0 ) vol(1,ib1) = vol(1,ib1) - flux
	  if( ib2 > 0 ) vol(1,ib2) = vol(1,ib2) + flux
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
	      write(6,'(a,3i5,g14.4)') 'non existing box: '
     +				,ib,lmax,l,vol(l,ib)
	    end if
	  end do
	  
	  volb = bvolume(ib)
	  do l=1,ltot
	    volf = vol(l,ib)
	    volv = dtbox * (wflux(l-1,ib) - wflux(l,ib))
	    vold = abs(volf + volv)
	    volr = vold/volb
	    !write(116,'(2i5,4g14.4)') ib,l,volf,volv,vold,volr
	    errmax = max(errmax,abs(volr))
	  end do

	  wtop = -dtbox * wflux(0,ib)
	  vole = vale(ib)
	  vold = abs(wtop - vole)
	  volr = vold/volb
	  !write(116,'(a,i5,4g14.4)') 'first layer: ',ib,wtop,vole,vold,volr
	  errmax = max(errmax,abs(volr))

	end do

	if( bdebug ) write(6,*) 'errmax 3d: ',errmax
	if( errmax .gt. eps ) then
	  write(6,*) 'errmax 3d: ',errmax
	end if

	end

c******************************************************************
c******************************************************************
c******************************************************************

