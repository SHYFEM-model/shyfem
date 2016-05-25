c
c $Id: noselab.f,v 1.8 2008-11-20 10:51:34 georg Exp $
c
c revision log :
c
c 18.11.1998    ggu     check dimensions with dimnos
c 06.04.1999    ggu     some cosmetic changes
c 03.12.2001    ggu     some extra output -> place of min/max
c 09.12.2003    ggu     check for NaN introduced
c 07.03.2007    ggu     easier call
c 08.11.2008    ggu     do not compute min/max in non-existing layers
c 07.12.2010    ggu     write statistics on depth distribution (depth_stats)
c 06.05.2015    ggu     noselab started
c 05.06.2015    ggu     many more features added
c 11.09.2015    ggu     split feature added
c 17.03.2016    ggu     outformat git added
c
c**************************************************************

	subroutine ouselab

	use clo
	use elabutil
	use elabtime

	use basin
        use mod_depth
        use mod_hydro
        use mod_hydro_baro
        use evgeom
        use levels
        use shyutil

c elaborates ous file

	implicit none

	real, allocatable :: u2v(:)
	real, allocatable :: v2v(:)
	real, allocatable :: zv(:)
	real, allocatable :: uprv(:,:)
	real, allocatable :: vprv(:,:)
	real, allocatable :: sv(:,:)
	real, allocatable :: dv(:,:)

	real, allocatable :: vars(:,:,:)
	real, allocatable :: varsaux(:,:,:)
	real, allocatable :: vars2d(:,:,:)
	integer, allocatable :: idims(:,:)

	real, allocatable :: hl(:)

	integer nread,nelab,nrec,nin,nwrite
	integer nndim,nvar,iv
	integer nvers
	integer nknous,nelous
	integer invar
	integer ierr
	integer it,ivar,itvar,itnew,itold,iaux
	integer i,l,k,lmax
	integer iano,ks
	integer ip,nb,naccum
	integer date,time
	character*80 title,name
	character*20 dline
	character*80 basnam,simnam
	real rnull
	real zmin,zmax
	real umin,umax,vmin,vmax
	real volume,area
	double precision dtime

	integer iapini
	integer ifem_open_file

c--------------------------------------------------------------

        nwrite = 0
	nread=0
	nelab=0
	nrec=0
	rnull=0.
	rnull=-999.
	bopen = .false.

	ks = -1			!write special node
	iano = -1		!no computation for this area code

c--------------------------------------------------------------
c set command line parameters
c--------------------------------------------------------------

	call elabutil_init('OUS')

	bneedbasin = bneedbasin .or. bsplit
	if( bneedbasin ) modeb = 3

	!--------------------------------------------------------------
	! open input files
	!--------------------------------------------------------------

	call ap_init(bask,modeb,0,0)

        call open_ous_type('.ous','old',nin)

        call ous_is_ous_file(nin,nvers)
        if( nvers .le. 0 ) then
          write(6,*) 'nvers: ',nvers
          stop 'error stop ouselab: not a valid ous file'
        end if

	call peek_ous_header(nin,nknous,nelous,nlv)

	if( bneedbasin ) then
	  if( nkn /= nknous .or. nel /= nelous ) goto 92
	else
	  nkn = nknous
	  nel = nelous
	end if

	call mod_hydro_baro_init(nel)
	call mod_depth_init(nkn,nel)
	call levels_init(nkn,nel,nlv)
	call mod_hydro_init(nkn,nel,nlv)

	if( nlv /= nlvdi ) stop 'error stop: problem.... nlv,nlvdi'

	allocate(u2v(nel),v2v(nel))
	allocate(hl(nlvdi))
	allocate(zv(nkn))
	allocate(uprv(nlvdi,nkn))
	allocate(vprv(nlvdi,nkn))
	allocate(sv(nlvdi,nkn))
	allocate(dv(nlvdi,nkn))

        call read_ous_header(nin,nkn,nel,nlvdi,ilhv,hlv,hev)
        call ous_get_params(nin,nkn,nel,nlv)

	call init_sigma_info(nlv,hlv)

	if( bneedbasin ) then
	  call outfile_make_hkv(nkn,nel,nen3v,hev,hkv)
	  call ev_init(nelous)
          call set_ev
	  call ilhe2k(nkn,nel,nen3v,ilhv,ilhkv)
	end if

	if( bverb ) call depth_stats(nkn,nlvdi,ilhkv)

	call handle_nodes

	!--------------------------------------------------------------
	! time management
	!--------------------------------------------------------------

	call ous_get_date(nin,date,time)
	call elabtime_date_and_time(date,time)

	!--------------------------------------------------------------
	! averaging
	!--------------------------------------------------------------

	call elabutil_set_averaging(nvar)
	!btrans = .false.			!FIXME
	nndim = max(3*nel,nkn)
	nvar = 4
	allocate(vars(nlvdi,nndim,nvar))
	allocate(varsaux(nlvdi,nndim,nvar))
	allocate(vars2d(1,nndim,nvar))
	allocate(idims(4,nvar))			!n,m,lmax,nvar

	vars = 0.
	vars2d = 0.
	idims = 0

	idims(1,:) = nel
	idims(1,1) = nkn
	idims(2,:) = 1
	idims(2,2) = 3
	idims(3,:) = 1
	idims(3,3) = nlvdi
	idims(3,4) = nlvdi
	idims(4,:) = 0		!not used

	if( btrans ) then
	  call shyutil_init_accum(nlvdi,nndim,nvar,istep)
	else
	  call shyutil_init_accum(1,1,1,1)
	end if

	!write(6,*) 'mode: ',mode,ifreq,istep

	!--------------------------------------------------------------
	! open output file
	!--------------------------------------------------------------

	boutput = boutput .or. btrans
	bopen = boutput .and. .not. bsplit

	if( bopen ) then
          call open_ous_file('out','new',nb)
          call ous_init(nb,0)
          call ous_clone_params(nin,nb)
	  if( b2d ) then
	    call ous_set_params(nb,0,0,1)
	  end if
          call write_ous_header(nb,ilhv,hlv,hev)
	end if

	if( outformat == 'gis' ) call gis_write_connect

c--------------------------------------------------------------
c loop on data
c--------------------------------------------------------------

	area = 0.
	volume = 0.
	it = 0
	if( .not. bquiet ) write(6,*)

	do

	  itold = it
	  vars = 0.

	  call new_read_record(nin,it,nlvdi,nndim,nvar,vars,ierr)

          if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
          if(ierr.ne.0) exit

	  dtime = it
	  nread=nread+1
	  nrec = nrec + 1

	  if( nrec == 1 ) itold = it
	  call ous_peek_record(nin,itnew,ierr)
	  if( ierr .ne. 0 ) itnew = it

	  if( .not. elabtime_check_time(it,itnew,itold) ) cycle

	  nelab=nelab+1

	  if( .not. bquiet ) then
	    dline = ' '
	    if( bdate ) call dtsgf(it,dline)
	    write(6,*) 'time : ',it,'  ',dline
	  end if

	  if( bwrite ) then

	    !call convert_2d(nlvdi,nndim,nvar,vars,vars2d)
	    !call transfer_uvz(1,nndim,nvar,vars2d,znv,zenv,unv,vnv)
	    !call mimar(znv,nkn,zmin,zmax,rnull)
            !call comp_vel2d(nel,hev,zenv,unv,vnv,u2v,v2v
!     +                          ,umin,vmin,umax,vmax)
	    !if( bneedbasin ) then
            !  call compute_volume_ia(iano,zenv,volume,area)
	    !end if

            !write(6,*) 'zmin/zmax : ',zmin,zmax
            !write(6,*) 'umin/umax : ',umin,umax
            !write(6,*) 'vmin/vmax : ',vmin,vmax
            !write(6,*) 'volume    : ',volume,area

	    call shy_write_min_max2(nlvdi,nkn,1,vars(1,1,1))
	    call shy_write_min_max2(nlvdi,3*nel,1,vars(1,1,2))
	    call shy_write_min_max2(nlvdi,nel,nlvdi,vars(1,1,3))
	    call shy_write_min_max2(nlvdi,nel,nlvdi,vars(1,1,4))

	    if( ks > 0 ) write(666,*) it,znv(ks),volume,area

	  end if

	  if( btrans ) then
	    do iv=1,nvar
	      call shy_time_aver(mode,iv,nread,ifreq,istep,nndim
     +                   ,idims(:,iv),threshold,vars(:,:,iv),boutput)
	    end do
	  end if

	  if( baverbas ) then
!	    call make_aver(nlvdi,nkn,ilhkv,cv3,vol3
!     +                          ,cmin,cmax,cmed,vtot)
!	    call write_aver(it,ivar,cmin,cmax,cmed,vtot)
	  end if

	  if( b2d ) then
	    call convert_2d(nlvdi,nndim,nvar,vars,vars2d)
	  end if

	  if( bsplit ) then
	    call transfer_uvz(nlvdi,nndim,nvar,vars
     +				,znv,zenv,utlnv,vtlnv)
            call transp2vel(nel,nkn,nlv,nlvdi,hev,zenv,nen3v
     +                          ,ilhv,hlv,utlnv,vtlnv
     +                          ,uprv,vprv)
	    call velzeta2scal(nel,nkn,nlv,nlvdi,nen3v,ilhkv
     +				,zenv,uprv,vprv
     +				,zv,sv,dv)
            call write_split(nin,nlvdi,ilhkv,hlv,hev,it
     +				,zv,sv,dv,uprv,vprv)
	    cycle
	  end if

	  if( boutput ) then
	    if( bverb ) write(6,*) 'writing to output: ',ivar
	    if( bsumvar ) ivar = 30
            nwrite = nwrite + 1
	    ierr = 0
	    if( outformat == 'gis' ) then
	      call transfer_uvz(nlvdi,nndim,nvar,vars
     +				,znv,zenv,utlnv,vtlnv)
              call transp2vel(nel,nkn,nlv,nlvdi,hev,zenv,nen3v
     +                          ,ilhv,hlv,utlnv,vtlnv
     +                          ,uprv,vprv)
	      call gis_write_hydro(it,nlvdi,ilhkv,znv,uprv,vprv)
	    else
	     if( b2d ) then
	      call new_write_record(nb,it,1,nndim,nvar,vars2d,ierr)
	     else
	      call new_write_record(nb,it,nlvdi,nndim,nvar,vars,ierr)
	     end if
	    end if
            if( ierr .ne. 0 ) goto 99
	  end if

	  if( bnodes ) then
	    call transfer_uvz(nlvdi,nndim,nvar,vars
     +				,znv,zenv,utlnv,vtlnv)
            call transp2vel(nel,nkn,nlv,nlvdi,hev,zenv,nen3v
     +                          ,ilhv,hlv,utlnv,vtlnv
     +                          ,uprv,vprv)
	    call write_nodes_vel(dtime,znv,uprv,vprv)
	  end if

	end do		!time do loop

c--------------------------------------------------------------
c end of loop on data
c--------------------------------------------------------------

c--------------------------------------------------------------
c final write of variables
c--------------------------------------------------------------

	if( btrans ) then
	  !write(6,*) 'istep,naccu: ',istep,naccu
	  do ip=1,istep
	   boutput = .false.
           ierr = 0
	   do iv=1,nvar
	    naccum = naccu(iv,ip)
	    !write(6,*) 'naccum: ',naccum
	    if( naccum > 0 ) then
	      !write(6,*) 'final aver: ',ip,naccum
	      call shy_time_aver(-mode,iv,ip,ifreq,istep,nndim
     +                   ,idims(:,iv),threshold,vars(:,:,iv),boutput)
              if( ierr .ne. 0 ) goto 99
	    end if
	   end do
	   if( boutput ) then
            nwrite = nwrite + 1
	    call new_write_record(nb,it,nlvdi,nndim,nvar,vars,ierr)
	   end if
	  end do
	end if

c--------------------------------------------------------------
c write final message
c--------------------------------------------------------------

	write(6,*)
	write(6,*) nread,' total records read'
	!write(6,*) nrec ,' unique time records read'
	write(6,*) nelab,' records elaborated'
	write(6,*) nwrite,' records written'
	write(6,*)

	if( boutput ) then
	  write(6,*) 'output written to file out.ous'
	end if

	call ap_get_names(basnam,simnam)
	write(6,*) 'names used: '
	write(6,*) 'basin: ',trim(basnam)
	write(6,*) 'simul: ',trim(simnam)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	stop
   92	continue
	write(6,*) 'incompatible basin: '
	write(6,*) 'nkn,nknous: ',nkn,nknous
	write(6,*) 'nel,nelous: ',nel,nelous
	stop 'error stop ouselab: parameter mismatch'
   99	continue
	write(6,*) 'error writing to file unit: ',nb
	stop 'error stop ouselab: write error'
	end

c***************************************************************
c***************************************************************
c***************************************************************

c***************************************************************

        subroutine write_split(nin,nlvddi,ilhkv,hlv,hev,it
     +				,zv,sv,dv,uv,vv)

        implicit none

        integer nin
        integer nlvddi
        integer ilhkv(1)
        real hlv(1)
        real hev(1)
	integer it
	real zv(*)
	real sv(nlvddi,*)
	real dv(nlvddi,*)
	real uv(nlvddi,*)
	real vv(nlvddi,*)

        integer nkn,nel,nlv,nvar
        integer ierr
	integer date,time
        character*80 name,title,femver

	integer, save :: icall = 0
	integer, save :: nz,ns,nd,nu,nv

        if( icall == 0 ) then      !open file
          call ous_get_params(nin,nkn,nel,nlv)
	  call ous_get_date(nin,date,time)
	  call ous_get_title(nin,title)
	  call ous_get_femver(nin,femver)

          call open_nos_file('zeta','new',nz)
          call nos_init(nz,0)
          call nos_set_params(nz,nkn,nel,1,1)
	  call nos_set_date(nz,date,time)
	  call nos_set_title(nz,title)
	  call nos_set_femver(nz,femver)

          call open_nos_file('speed','new',ns)
          call nos_init(ns,0)
	  call nos_clone_params(nz,ns)
          call nos_set_params(ns,nkn,nel,nlv,1)

          call open_nos_file('dir','new',nd)
          call nos_init(nd,0)
	  call nos_clone_params(nz,nd)
          call nos_set_params(nd,nkn,nel,nlv,1)

          call open_nos_file('uvel','new',nu)
          call nos_init(nu,0)
	  call nos_clone_params(nz,nu)
          call nos_set_params(nu,nkn,nel,nlv,1)

          call open_nos_file('vvel','new',nv)
          call nos_init(nv,0)
	  call nos_clone_params(nz,nv)
          call nos_set_params(nv,nkn,nel,nlv,1)

          call write_nos_header(nz,ilhkv,hlv,hev)
          call write_nos_header(ns,ilhkv,hlv,hev)
          call write_nos_header(nd,ilhkv,hlv,hev)
          call write_nos_header(nu,ilhkv,hlv,hev)
          call write_nos_header(nv,ilhkv,hlv,hev)
        end if

	icall = icall + 1

	call nos_write_record(nz,it,1,1,ilhkv,zv,ierr)
	call nos_write_record(ns,it,6,nlvddi,ilhkv,sv,ierr)
	call nos_write_record(nd,it,7,nlvddi,ilhkv,dv,ierr)
	call nos_write_record(nu,it,2,nlvddi,ilhkv,uv,ierr)
	call nos_write_record(nv,it,2,nlvddi,ilhkv,vv,ierr)

        end

c***************************************************************
c***************************************************************
c***************************************************************

	subroutine new_read_record(nb,it,nlvddi,nndim,nvar,vars,ierr)

	use basin
	use levels

	implicit none

	integer nb
	integer it
	integer nlvddi,nndim,nvar
	real vars(nlvddi,nndim,nvar)
	integer ierr

	real znv(nkn)
	real zenv(3*nel)
	real utlnv(nlvddi,nel)
	real vtlnv(nlvddi,nel)

	utlnv = 0.
	vtlnv = 0.
        call ous_read_record(nb,it,nlvddi,ilhv,znv,zenv
     +                          ,utlnv,vtlnv,ierr)

	vars(1,1:nkn,1)   = znv(1:nkn)
	vars(1,1:3*nel,2) = zenv(1:3*nel)
	vars(:,1:nel,3)   = utlnv(:,1:nel)
	vars(:,1:nel,4)   = vtlnv(:,1:nel)

	end

c***************************************************************

	subroutine new_write_record(nb,it,nlvddi,nndim,nvar,vars,ierr)

	use basin
	use levels

	implicit none

	integer nb
	integer it
	integer nlvddi,nndim,nvar
	real vars(nlvddi,nndim,nvar)
	integer ierr

	real znv(nkn)
	real zenv(3*nel)
	real utlnv(nlvddi,nel)
	real vtlnv(nlvddi,nel)

	ierr = 0

	znv(1:nkn)     = vars(1,1:nkn,1)
	zenv(1:3*nel)  = vars(1,1:3*nel,2)
	utlnv(:,1:nel) = vars(:,1:nel,3)
	vtlnv(:,1:nel) = vars(:,1:nel,4)

        call ous_write_record(nb,it,nlvddi,ilhv
     +			,znv,zenv,utlnv,vtlnv,ierr)

	end

c***************************************************************

	subroutine convert_2d(nlvddi,nndim,nvar,vars,vars2d)

	use basin

	implicit none

	integer nlvddi,nndim,nvar
	real vars(nlvddi,nndim,nvar)
	real vars2d(1,nndim,nvar)

	integer i

	vars2d(1,:,1) = vars(1,:,1)
	vars2d(1,:,2) = vars(1,:,2)
	do i=1,nel
	  vars2d(1,i,3) = sum(vars(:,i,3),1)
	  vars2d(1,i,4) = sum(vars(:,i,4),1)
	end do

	end

c***************************************************************

	subroutine transfer_uvz(nlvddi,nndim,nvar,vars,znv,zenv,uv,vv)

	use basin

	implicit none

	integer nlvddi,nndim,nvar
	real vars(nlvddi,nndim,nvar)
	real znv(nkn)
	real zenv(3*nel)
	real uv(nlvddi,nel)
	real vv(nlvddi,nel)

	znv(1:nkn)     = vars(1,1:nkn,1)
	zenv(1:3*nel)  = vars(1,1:3*nel,2)
	uv(:,1:nel)    = vars(:,1:nel,3)
	vv(:,1:nel)    = vars(:,1:nel,4)

	end

c***************************************************************

        subroutine shy_write_min_max2(nlvdi,nn,lmax,cv3)

        implicit none

        integer nlvdi,nn,lmax
        real cv3(nlvdi,nn)

        integer l
        real rnull
        real cmin,cmax,cmed
        real cv2(nn)

        do l=1,lmax
          cv2=cv3(l,:)
          call mimar(cv2,nn,cmin,cmax,rnull)
          call aver(cv2,nn,cmed,rnull)
          call check1Dr(nn,cv2,0.,-1.,"NaN check","cv2")
          write(6,1000) 'l,min,max,aver : ',l,cmin,cmax,cmed
        end do

 1000   format(a,i5,3g16.6)
        end

c***************************************************************
c***************************************************************
c***************************************************************

	subroutine make_hydro_basin_aver(dtime,nndim,vars
     +			,znv,uprv,vprv,sv)

c not finished...

	use basin
	use levels
	use mod_depth

	implicit none

	integer, parameter :: nvar = 4
	double precision dtime
	integer nndim
	real vars(nlvdi,nndim,nvar)
	real znv(nkn)
	real uprv(nlvdi,nkn)
	real vprv(nlvdi,nkn)
	real sv(nlvdi,nkn)

	integer it
	real zv(nkn)
	real zenv(3*nel)
	real utlnv(nlvdi,nel)
	real vtlnv(nlvdi,nel)
	real dv(nlvdi,nkn)

	call transfer_uvz(nlvdi,nndim,nvar,vars
     +				,znv,zenv,utlnv,vtlnv)
        call transp2vel(nel,nkn,nlv,nlvdi,hev,zenv,nen3v
     +                          ,ilhv,hlv,utlnv,vtlnv
     +                          ,uprv,vprv)
	call velzeta2scal(nel,nkn,nlv,nlvdi,nen3v,ilhkv
     +				,zenv,uprv,vprv
     +				,zv,sv,dv)

!	    call make_aver(nlvdi,nkn,ilhkv,cv3,vol3
!     +                          ,cmin,cmax,cmed,vtot)
!	    call write_aver(it,ivar,cmin,cmax,cmed,vtot)

!	    call make_aver(nlvdi,nkn,ilhkv,cv3,vol3
!     +                          ,cmin,cmax,cmed,vtot)
!	    call write_aver(it,ivar,cmin,cmax,cmed,vtot)

	end

c***************************************************************


