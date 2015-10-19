!
! $Id: noselab.f,v 1.8 2008-11-20 10:51:34 georg Exp $
!
! revision log :
!
! 18.11.1998    ggu     check dimensions with dimnos
! 06.04.1999    ggu     some cosmetic changes
! 03.12.2001    ggu     some extra output -> place of min/max
! 09.12.2003    ggu     check for NaN introduced
! 07.03.2007    ggu     easier call
! 08.11.2008    ggu     do not compute min/max in non-existing layers
! 07.12.2010    ggu     write statistics on depth distribution (depth_stats)
! 06.05.2015    ggu     noselab started
! 05.06.2015    ggu     many more features added
! 10.09.2015    ggu     std and rms for averaging implemented
! 11.09.2015    ggu     write in gis format
! 23.09.2015    ggu     handle more than one file (look for itstart)
! 16.10.2015    ggu     started shyelab
!
!**************************************************************

	subroutine shyelab1

	use clo
	use elabutil
	use shyfile

        use basin
        use mod_depth
        use evgeom
        use levels

	implicit none

	include 'param.h'

	integer, parameter :: ndim = 1000
	integer iusplit(ndim)

	real, allocatable :: cv2(:)
	real, allocatable :: cv3(:,:)
	real, allocatable :: cv3all(:,:,:)
	real, allocatable :: vol3(:,:)

	real, allocatable :: znv(:)
	real, allocatable :: zenv(:,:)

	integer, allocatable :: ivars(:)

	integer, allocatable :: naccu(:)
	double precision, allocatable :: accum(:,:,:)
	double precision, allocatable :: std(:,:,:,:)

	real, allocatable :: hl(:)

	logical bnextfile
	integer nwrite,nread,nelab,nrec,nin,nold
	integer nvers
	integer nknnos,nelnos,nvar,npr
	integer ierr
	!integer it,itvar,itnew,itold,itstart
	integer it
	integer ivar,iaux
	integer i,j,l,k,lmax,node
	integer ip,nb,naccum
	integer ifile,ftype
	integer id,idout,idold
	integer n,m,nndim
	character*80 title,name,file
	character*20 dline
	character*80 basnam,simnam
	real rnull
	real cmin,cmax,cmed,vtot
	double precision dtime,dtstart,dtnew,dtvar
	double precision atime,atstart,atold,atnew

	integer iapini
	integer ifem_open_file

!--------------------------------------------------------------
! initialize everything
!--------------------------------------------------------------

	nread=0
	nwrite=0
	nelab=0
	nrec=0
	rnull=0.
	rnull=-1.
	bopen = .false.

	!--------------------------------------------------------------
	! set command line parameters
	!--------------------------------------------------------------

	call elabutil_init('SHY')

	!--------------------------------------------------------------
	! open input files
	!--------------------------------------------------------------

	ifile = 1
	idold = 0
	call open_next_file(ifile,idold,id)
	call get_start_of_next_file(ifile+1,atstart,bnextfile)

	!--------------------------------------------------------------
	! set up params and arrays
	!--------------------------------------------------------------

	call shy_get_params(id,nkn,nel,npr,nlv,nvar)
	call shy_get_ftype(id,ftype)

	if( ftype < 1 .or. ftype > 2 ) goto 76	!relax later
	nndim = nkn
	if( ftype == 1 ) nndim = 3*nel	!needed for zenv

	call shy_info(id)

        call basin_init(nkn,nel)
        call levels_init(nkn,nel,nlv)
        call mod_depth_init(nkn,nel)
	call shy_copy_basin_from_shy(id)
	call shy_copy_levels_from_shy(id)

	allocate(cv2(nndim))
	allocate(cv3(nlv,nndim))
	allocate(vol3(nlv,nndim))
	allocate(cv3all(nlv,nndim,nvar))
        allocate(hl(nlv))
	allocate(ivars(nvar))

	call init_sigma_info(nlv,hlv)

	call outfile_make_depth(nkn,nel,nen3v,hm3v,hev,hkv)
	call init_volume(nlvdi,nkn,nel,nlv,nen3v,ilhkv
     +                          ,hlv,hev,hl,vol3)

	if( bverb ) call depth_stats(nkn,nlvdi,ilhkv)

	call handle_nodes	!single node output

	!--------------------------------------------------------------
	! time management
	!--------------------------------------------------------------

	call shy_get_date(id,date,time)
	call elabutil_date_and_time	!this also sets datetime

	!--------------------------------------------------------------
	! averaging
	!--------------------------------------------------------------

	call elabutil_set_averaging(nvar)

	if( btrans ) then
	  allocate(naccu(istep))
	  allocate(accum(nlvdi,nndim,istep))
	  allocate(std(nlvdi,nndim,16,istep))	!also used for directions
	else
	  allocate(naccu(1))
	  allocate(accum(1,1,1))
	  allocate(std(1,1,1,1))	!also used for directions
	end if
	naccum = 0
	naccu = 0
	accum = 0.
	std = 0.

	!write(6,*) 'mode: ',mode,ifreq,istep

	!--------------------------------------------------------------
	! open output file
	!--------------------------------------------------------------

	iusplit = 0

	boutput = boutput .or. btrans
	bopen = boutput .and. .not. bsplit

	if( bopen ) then
	  idout = shy_init('out.shy')
          call shy_clone(id,idout)
	  if( b2d ) then
	    call shy_set_params(idout,nkn,nel,npr,1,nvar)
	  end if
          call shy_write_header(idout,ierr)
	  if( ierr /= 0 ) goto 75
	end if

	if( outformat == 'gis' ) call gis_write_connect

!--------------------------------------------------------------
! loop on data
!--------------------------------------------------------------

	dtvar = 0.
	dtime = 0.
	call shy_peek_record(id,dtime,iaux,iaux,iaux,iaux,ierr)
	call fem_file_convert_time(datetime,dtime,atime)

	cv3 = 0.

	do

	 atold = atime

	 call read_records(id,dtime,nvar,nndim,ivars,nlvdi
     +				,cv3,cv3all,ierr)

         if(ierr.ne.0) then	!EOF - see if we have to read another file
	   if( ierr > 0 ) exit
	   if( .not. bnextfile ) exit
	   idold = id
	   ifile = ifile + 1
	   call open_next_file(ifile,idold,id)
	   call get_start_of_next_file(ifile+1,atstart,bnextfile)
	   call nos_get_date(nin,date,time)
	   call elabutil_date_and_time
	   atime = atold	!reset time of last successfully read record
	   cycle
	 end if

	 nread=nread+nvar
	 nrec = nrec + 1
	 call fem_file_convert_time(datetime,dtime,atime)

	 call shy_peek_record(id,dtnew,iaux,iaux,iaux,iaux,ierr)
	 call fem_file_convert_time(datetime,dtnew,atnew)
	 if( ierr .ne. 0 ) atnew = atime

	 if( .not. elabutil_check_time_a(atime,atnew,atold) ) cycle

	 do i=1,nvar

	  ivar = ivars(i)
	  cv3(:,:) = cv3all(:,:,i)

	  nelab=nelab+1

	  if( .not. bquiet ) then
	    dline = ' '
	    !if( bdate ) call dtsgf(it,dline)
	    if( bdate ) call dts_format_abs_time(atime,dline)
	    write(6,*) 'time : ',dtime,' ',dline,'  ivar : ',ivar
	  end if

	  if( bwrite ) then
	    do l=1,nlv
	      do k=1,nkn
	        cv2(k)=cv3(l,k)
	        if( l .gt. ilhkv(k) ) cv2(k) = rnull
	      end do
	      call mimar(cv2,nkn,cmin,cmax,rnull)
              call aver(cv2,nkn,cmed,rnull)
              call check1Dr(nkn,cv2,0.,-1.,"NaN check","cv2")
	      write(6,*) 'l,min,max,aver : ',l,cmin,cmax,cmed
	    end do
	  end if

	  if( btrans ) then
	    call nos_time_aver(mode,i,ifreq,istep,nkn,nlvdi
     +				,naccu,accum,std,threshold,cv3,boutput)
	  end if

	  if( baverbas ) then
	    call make_aver(nlvdi,nkn,ilhkv,cv3,vol3
     +                          ,cmin,cmax,cmed,vtot)
	    !call write_aver(it,ivar,cmin,cmax,cmed,vtot)
	  end if

	  if( b2d ) then
	    call make_vert_aver(nlvdi,nkn,ilhkv,cv3,vol3,cv2)
	  end if

	  if( bsplit ) then
            call get_split_iu(ndim,iusplit,ivar,nin,ilhkv,hlv,hev,nb)
	  end if

	  if( boutput ) then
	    nwrite = nwrite + 1
	    if( bverb ) write(6,*) 'writing to output: ',ivar
	    if( bsumvar ) ivar = 30
	    if( b2d ) then
	      call shy_write_scalar_record2d(idout,dtime,ivar,cv2)
	    else
	      call shy_write_scalar_record(idout,dtime,ivar,nlvdi,cv3)
	    end if
            if( ierr .ne. 0 ) goto 99
	  end if

	  if( bnodes ) then
	    call write_nodes(dtime,ivar,nlvdi,cv3)
	  end if

	 end do		!loop on ivar
	end do		!time do loop

!--------------------------------------------------------------
! end of loop on data
!--------------------------------------------------------------

!--------------------------------------------------------------
! final write of variables
!--------------------------------------------------------------

! next not working for b2d ... FIXME

	if( btrans ) then
	  !write(6,*) 'istep,naccu: ',istep,naccu
	  do ip=1,istep
	    naccum = naccu(ip)
	    !write(6,*) 'naccum: ',naccum
	    if( naccum > 0 ) then
	      nwrite = nwrite + 1
	      write(6,*) 'final aver: ',ip,naccum
	      call nos_time_aver(-mode,ip,ifreq,istep,nkn,nlvdi
     +				,naccu,accum,std,threshold,cv3,boutput)
	      if( bsumvar ) ivar = 30
	      call shy_write_scalar_record(idout,dtime,ivar,nlvdi,cv3)
	    end if
	  end do
	end if

!--------------------------------------------------------------
! write final message
!--------------------------------------------------------------

	write(6,*)
	write(6,*) nread, ' records read'
	write(6,*) nrec , ' unique time records read'
	write(6,*) nelab, ' records elaborated'
	write(6,*) ifile, ' files read'
	write(6,*) nwrite,' records written'
	write(6,*)

	if( bsplit ) then
	  write(6,*) 'output written to following files: '
	  do ivar=1,ndim
	    nb = iusplit(ivar)
	    if( nb .gt. 0 ) then
              write(name,'(i4)') ivar
	      write(6,*) trim(adjustl(name))//'.nos'
	      close(nb)
	    end if
	  end do
	else if( boutput ) then
	  write(6,*) 'output written to file out.shy'
	  close(nb)
	end if

	!call ap_get_names(basnam,simnam)
	!write(6,*) 'names used: '
	!write(6,*) 'basin: ',trim(basnam)
	!write(6,*) 'simul: ',trim(simnam)

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	stop
   74	continue
	stop 'error stop shyelab: general error...'
   75	continue
	write(6,*) 'error writing header, ierr = ',ierr
	write(6,*) 'file = ',trim(file)
	stop 'error stop shyelab: writing header'
   76	continue
	write(6,*) 'ftype = ',ftype,'  expecting ',2
	stop 'error stop shyelab: ftype'
   77	continue
	write(6,*) 'error reading header, ierr = ',ierr
	write(6,*) 'file = ',trim(file)
	stop 'error stop shyelab: reading header'
   85	continue
	write(6,*) 'dtime,dtvar,i,ivar,nvar: ',dtime,dtvar,i,ivar,nvar
	stop 'error stop shyelab: time mismatch'
   92	continue
	write(6,*) 'incompatible basin: '
	write(6,*) 'nkn,nknnos: ',nkn,nknnos
	write(6,*) 'nel,nelnos: ',nel,nelnos
	stop 'error stop shyelab: parameter mismatch'
   99	continue
	write(6,*) 'error writing to file unit: ',nb
	stop 'error stop shyelab: write error'
	end

!***************************************************************
!***************************************************************
!***************************************************************

	subroutine open_next_file(ifile,idold,id)

	use shyfile
	use clo

	implicit none

	integer ifile
	integer idold,id

	integer ierr
	character*80 file

	call clo_get_file(ifile,file)
	if( .not. shy_is_shy_file(file) ) then
	  write(6,*) 'file: ',trim(file)
	  stop 'error stop shyelab: not a valid shy file'
	end if

	id = shy_init(file)
	if( id == 0 ) then
	  stop 'error stop open_next_file: internal error (1)'
	end if

	write(6,*) '================================'
	write(6,*) 'reading file: ',trim(file)
	write(6,*) '================================'

	call shy_read_header(id,ierr)
	if( ierr /= 0 ) then
	  write(6,*) 'error reading header, ierr = ',ierr
	  write(6,*) 'file = ',trim(file)
	  stop 'error stop open_next_file: reading header'
	end if

	if( idold > 0 ) then
	  if( .not. shy_are_compatible(idold,id) ) then
	    write(6,*) 'files are not compatible'
	    stop 'error stop open_next_file: not compatible'
	  end if
	  call shy_close(idold)
	end if

	!call shy_info(id)

	end

!***************************************************************

	subroutine get_start_of_next_file(ifile,atstart,bok)

	use clo

	implicit none

	integer ifile
	double precision atstart
	logical bok

	integer datetime(2)
	double precision dtime
	character*80 file

	atstart = -1.

	call clo_get_file(ifile,file)
	if( file == ' ' ) return

	call shy_get_tstart(file,datetime,dtime,bok)
	if( .not. bok ) return

	call fem_file_convert_time(datetime,dtime,atstart)

	end 

!***************************************************************

	subroutine read_records(id,dtime,nvar,nndim,ivars,nlvddi
     +				,cv3,cv3all,ierr)

	use elabutil
	use shyfile

	implicit none

	integer id
	double precision dtime
	integer nvar,nndim
	integer ivars(nvar)
	integer nlvddi
	real cv3(nlvddi,nndim)
	real cv3all(nlvddi,nndim,nvar)
	integer ierr

	integer idims(3,nvar)
	integer nexp

	integer i
	integer ivar,n,m,lmax
	double precision dtvar

	do i=1,nvar
	  call shy_read_record(id,dtime,ivar,n,m,lmax,nlvddi,cv3,ierr)
	  nexp = n * m
	  idims(1,i) = n
	  idims(2,i) = m
	  idims(3,i) = lmax
	  if( nexp > nndim ) goto 74
	  if( lmax > nlvddi ) goto 74
          if( ierr .gt. 0 ) goto 75
          if( ierr .ne. 0 ) exit
	  if( i == 1 ) dtvar = dtime
	  if( dtvar /= dtime ) goto 85
	  ivars(i) = ivar
	  cv3all(:,:,i) = cv3(:,:) * fact
	end do

	if( ierr /= 0 .and. i .ne. 1 ) goto 76

	return
   74	continue
	write(6,*) n,m,nexp,nndim,lmax,nlvddi
	stop 'error stop shyelab: dimension error...'
   75	continue
        write(6,*) 'error in reading file : ',ierr
	stop 'error stop shyelab: reading file'
   76	continue
        write(6,*) 'end of file between variables'
	stop 'error stop shyelab: EOF unexpected'
   85	continue
	write(6,*) 'dtime,dtvar,i,ivar,nvar: ',dtime,dtvar,i,ivar,nvar
	stop 'error stop shyelab: time mismatch'
	end

!***************************************************************
!***************************************************************
!***************************************************************

	subroutine nos_time_aver0(mode,nread,ifreq,istep,nkn,nlvddi
     +				,naccu,accum,std,threshold,cv3,bout)

! mode:  1:aver  2:sum  3:min  4:max  5:std  6:rms  7:thres  8:averdir
!
! mode negative: only transform, do not accumulate

	implicit none

	integer mode
	integer nread,ifreq,istep
	integer nkn,nlvddi
	integer naccu(istep)
	double precision accum(nlvddi,nkn,istep)
	double precision std(nlvddi,nkn,16,istep)
	double precision threshold
	real cv3(nlvddi,nkn)
	logical bout

	integer ip,naccum,mmode
	integer k,l,id,idmax
	double precision dmax

	if( mode .eq. 0 ) return

	bout = .false.
	ip = mod(nread,istep)
	if( ip .eq. 0 ) ip = istep

	!write(6,*) 'ip: ',ip,istep,nread,mode

	if( mode == 1 .or. mode == 2 ) then
	  accum(:,:,ip) = accum(:,:,ip) + cv3(:,:)
	else if( mode == 3 ) then
	  do k=1,nkn
	    do l=1,nlvddi
	      accum(l,k,ip) = min(accum(l,k,ip),cv3(l,k))
	    end do
	  end do
	else if( mode == 4 ) then
	  do k=1,nkn
	    do l=1,nlvddi
	      accum(l,k,ip) = max(accum(l,k,ip),cv3(l,k))
	    end do
	  end do
	else if( mode == 5 ) then
	  accum(:,:,ip) = accum(:,:,ip) + cv3(:,:)
	  std(:,:,1,ip) = std(:,:,1,ip) + cv3(:,:)**2
	else if( mode == 6 ) then
	  accum(:,:,ip) = accum(:,:,ip) + cv3(:,:)**2
	else if( mode == 7 ) then
	  where( cv3(:,:) >= threshold )
	    accum(:,:,ip) = accum(:,:,ip) + 1.
	  end where
	else if( mode == 8 ) then
	  do k=1,nkn
	    do l=1,nlvddi
	      id = nint( cv3(l,k)/22.5 )
	      if( id == 0 ) id = 16
	      if( id < 0 .or. id > 16 ) stop 'error stop: direction'
	      std(l,k,id,ip) = std(l,k,id,ip) + 1.
	    end do
	  end do
	end if

	if( mode > 0 ) naccu(ip) = naccu(ip) + 1
	!write(6,*) '... ',ifreq,mode,ip,istep,naccu(ip)

	if( naccu(ip) == ifreq .or. mode < 0 ) then	!here ip == 1
	  naccum = max(1,naccu(ip))
	  mmode = abs(mode)
	  if( mmode == 3 ) naccum = 1			!min
	  if( mmode == 4 ) naccum = 1			!max
	  if( mmode == 7 ) naccum = 1			!threshold
	  if( naccum > 0 ) cv3(:,:) = accum(:,:,ip)/naccum
	  if( mmode == 5 ) then
	    cv3(:,:) = sqrt( std(:,:,1,ip)/naccum - cv3(:,:)**2 )
	  else if( mmode == 6 ) then
	    cv3(:,:) = sqrt( cv3(:,:) )
	  else if( mmode == 8 ) then
	    do k=1,nkn
	      do l=1,nlvddi
		dmax = 0.
		idmax = 0
	        do id=1,16
		  if( std(l,k,id,ip) > dmax ) then
		    idmax = id
		    dmax = std(l,k,id,ip)
		  end if
		end do
		if( idmax == 16 ) idmax = 0
		cv3(l,k) = idmax * 22.5
	      end do
	    end do
	  end if
	  write(6,*) 'averaging: ',ip,naccum,naccu(ip)
	  bout = .true.
	  naccu(ip) = 0
	  accum(:,:,ip) = 0.
	  std(:,:,:,ip) = 0.
	end if

	end

!***************************************************************

        subroutine get_split_iu0(ndim,iu,ivar,nin,ilhkv,hlv,hev,nb)

        implicit none

        integer ndim
        integer iu(ndim)
        integer ivar
        integer nin
        integer ilhkv(1)
        real hlv(1)
        real hev(1)
	integer nb		!unit to use for writing (return)

        integer nkn,nel,nlv,nvar
        integer ierr
        character*80 name

        if( ivar > ndim ) then
          write(6,*) 'ndim,ivar: ',ndim,ivar
          stop 'error stop: ndim'
        end if

        if( iu(ivar) .le. 0 ) then      !open file
          write(name,'(i4)') ivar
          call open_nos_file(name,'new',nb)
          call nos_init(nb,0)
          call nos_clone_params(nin,nb)
          call nos_get_params(nb,nkn,nel,nlv,nvar)
          call nos_set_params(nb,nkn,nel,nlv,1)
          call write_nos_header(nb,ilhkv,hlv,hev)
          iu(ivar) = nb
        end if

        nb = iu(ivar)

        end

!***************************************************************

	subroutine noselab_write_record0(nb,it,ivar,nlvdi,ilhkv,cv,ierr)

        use elabutil

	implicit none

	integer nb,it,ivar,nlvdi
	integer ilhkv(nlvdi)
	real cv(nlvdi,*)
	integer ierr

	ierr = 0

	if( outformat == 'nos' .or. outformat == 'native') then
          call nos_write_record(nb,it,ivar,nlvdi,ilhkv,cv,ierr)
	else if( outformat == 'gis' ) then
          call gis_write_record(nb,it,ivar,nlvdi,ilhkv,cv)
	else
	  write(6,*) 'output format unknown: ',outformat
	  stop 'error stop noselab_write_record: output format'
	end if

        end

!***************************************************************

        subroutine gis_write_record0(nb,it,ivar,nlvdi,ilhkv,cv)

! writes one record to file nb (3D)

        use basin

        implicit none

        integer nb,it,ivar,nlvdi
        integer ilhkv(nlvdi)
        real cv(nlvdi,*)

        integer k,l,lmax
	integer nout
        real x,y
	character*80 format,name
	character*20 line
	character*3 var

	integer ifileo

	call dtsgf(it,line)
	call i2s0(ivar,var)

	name = 'extract_'//var//'_'//line//'.gis'
        nout = ifileo(60,name,'form','new')
	!write(6,*) 'writing: ',trim(name)

        write(nout,*) it,nkn,ivar,line

        do k=1,nkn
          lmax = ilhkv(k)
          x = xgv(k)
          y = ygv(k)

	  write(format,'(a,i5,a)') '(i10,2g14.6,i5,',lmax,'g14.6)'
          write(nout,format) k,x,y,lmax,(cv(l,k),l=1,lmax)
        end do

	close(nout)

        end

!***************************************************************

        subroutine gis_write_connect0

! writes connectivity

        use basin

        implicit none

	integer ie,ii

	open(1,file='connectivity.gis',form='formatted',status='unknown')

	write(1,*) nel
	do ie=1,nel
	  write(1,*) ie,(nen3v(ii,ie),ii=1,3)
	end do

	close(1)

	end

!***************************************************************

