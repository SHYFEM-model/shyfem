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
	use shyutil

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

	integer, allocatable :: ivars(:,:)
	integer, allocatable :: il(:)

	real, allocatable :: hl(:)

	logical bnextfile
	integer nwrite,nread,nelab,nrec,nin,nold
	integer nvers
	integer nknnos,nelnos,nvar,npr
	integer ierr
	!integer it,itvar,itnew,itold,itstart
	integer it
	integer ivar,iaux
	integer iv,j,l,k,lmax,node
	integer ip,nb
	integer ifile,ftype
	integer id,idout,idold
	integer n,m,nndim,nn
	integer naccum
	character*80 title,name,file
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
	bzeta = .false.		!file has zeta information

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

	if( ftype == 1 ) then		!OUS
	  if( nvar /= 4 ) goto 71
	  nndim = 3*nel
	  allocate(il(nel))
	  il = ilhv
	else if( ftype == 2 ) then	!NOS
	  nndim = nkn
	  allocate(il(nkn))
	  il = ilhkv
	else
	  goto 76	!relax later
	end if

	call shy_info(id)

        call basin_init(nkn,nel)
        call levels_init(nkn,nel,nlv)
        call mod_depth_init(nkn,nel)
	call set_ev
	call shy_copy_basin_from_shy(id)
	call shy_copy_levels_from_shy(id)

	allocate(cv2(nndim))
	allocate(cv3(nlv,nndim))
	allocate(cv3all(nlv,nndim,0:nvar))

        allocate(hl(nlv))
	allocate(ivars(4,nvar))

	call shyutil_init(nkn,nel,nlv)

	call init_sigma_info(nlv,hlv)

	call shy_make_area
	call outfile_make_depth(nkn,nel,nen3v,hm3v,hev,hkv)

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

	call elabutil_set_averaging(nvar)	!sets btrans

	if( btrans ) then
	  call shyutil_init_accum(nlvdi,nndim,nvar,istep)
	else
	  call shyutil_init_accum(1,1,1,1)
	end if

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
	cv3all = 0.

	do

	 atold = atime

	 call read_records(id,dtime,nvar,nndim,nlvdi,ivars
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

	 nread = nread + nvar
	 nrec = nrec + 1
	 call fem_file_convert_time(datetime,dtime,atime)

	 call shy_peek_record(id,dtnew,iaux,iaux,iaux,iaux,ierr)
	 call fem_file_convert_time(datetime,dtnew,atnew)
	 if( ierr .ne. 0 ) atnew = atime

	 if( .not. elabutil_check_time_a(atime,atnew,atold) ) cycle

	 call shy_make_zeta(ftype)
	 call shy_make_volume

	 do iv=1,nvar

	  ivar = ivars(4,iv)
	  lmax = ivars(3,iv)
	  nn = ivars(1,iv) * ivars(2,iv)

	  cv3(:,:) = cv3all(:,:,iv)

	  nelab=nelab+1

	  if( .not. bquiet ) then
	    call shy_write_time(bdate,dtime,atime,ivar)
	  end if

	  if( bwrite ) then
	    call shy_write_min_max(nlvdi,nn,lmax,cv3)
	  end if

	  if( btrans ) then
	    call shy_time_aver(mode,iv,nrec,ifreq,istep,nndim
     +			,ivars,threshold,cv3,boutput)
	  end if

	  if( baverbas ) then
	    call shy_make_aver(ivars(:,iv),nndim,cv3
     +                          ,cmin,cmax,cmed,vtot)
	    !call write_aver(it,ivar,cmin,cmax,cmed,vtot)
	  end if

	  if( b2d ) then
	    call shy_make_vert_aver(ivars(:,iv),nndim,cv3,cv2)
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
	    naccum = naccu(iv,ip)
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
   71	continue
	write(6,*) 'ftype = ',ftype,'  nvar = ',nvar
	write(6,*) 'nvar should be 4'
	stop 'error stop shyelab: ftype,nvar'
   74	continue
	stop 'error stop shyelab: general error...'
   75	continue
	write(6,*) 'error writing header, ierr = ',ierr
	write(6,*) 'file = ',trim(file)
	stop 'error stop shyelab: writing header'
   76	continue
	write(6,*) 'ftype = ',ftype,'  expecting 1 or 2'
	stop 'error stop shyelab: ftype'
   77	continue
	write(6,*) 'error reading header, ierr = ',ierr
	write(6,*) 'file = ',trim(file)
	stop 'error stop shyelab: reading header'
   85	continue
	write(6,*) 'dtime,dtvar,iv,ivar,nvar: ',dtime,dtvar,iv,ivar,nvar
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

	subroutine shy_write_time(bdate,dtime,atime,ivar)

	implicit none

	logical bdate
	double precision dtime,atime
	integer ivar

	character*20 dline

	dline = ' '
	!if( bdate ) call dtsgf(it,dline)
	if( bdate ) call dts_format_abs_time(atime,dline)
	write(6,*) 'time : ',dtime,' ',dline,'  ivar : ',ivar

	end

!***************************************************************

	subroutine shy_write_min_max(nlvdi,nn,lmax,cv3)

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
	  write(6,*) 'l,min,max,aver : ',l,cmin,cmax,cmed
	end do

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
!***************************************************************
!***************************************************************

	subroutine read_records(id,dtime,nvar,nndim,nlvddi
     +				,ivars,cv3,cv3all,ierr)

	use elabutil
	use shyfile
	use shyutil

	implicit none

	integer id
	double precision dtime
	integer nvar,nndim
	integer nlvddi
	integer nkn,nel
	integer ivars(4,nvar)
	real cv3(nlvddi,nndim)
	real cv3all(nlvddi,nndim,0:nvar)
	integer ierr

	integer nexp

	integer iv
	integer ivar,n,m,lmax
	logical bfirst
	double precision dtvar

	shy_znv = 0.
	shy_zenv = 0.

	iv = 0
	bfirst = .true.
	bzeta = .false.

	do
	  iv = iv + 1
	  if( iv > nvar ) exit
	  call shy_read_record(id,dtime,ivar,n,m,lmax,nlvddi,cv3,ierr)
	  nexp = n * m
	  bzeta = ivar == -1
	  if( bzeta ) iv = iv - 1
	  ivars(1,iv) = n
	  ivars(2,iv) = m
	  ivars(3,iv) = lmax
	  if( nexp > nndim ) goto 74
	  if( lmax > nlvddi ) goto 74
          if( ierr .gt. 0 ) goto 75
          if( ierr .ne. 0 ) exit
	  if( bfirst ) dtvar = dtime
	  if( dtvar /= dtime ) goto 85
	  ivars(4,iv) = ivar
	  cv3all(:,:,iv) = cv3(:,:) * fact
	  if( abs(ivar) == 1 ) then		! water level
	    if( ivar == -1 .or. iv == 1 ) then
	      shy_znv = cv3(1,1:n)
	    else
	      !zenv = cv3(1,1:3*n)
	      shy_zenv = reshape(cv3(1,1:3*n),(/3,n/))	!FIXME
	    end if
	  end if
	end do

	if( ierr /= 0 .and. iv .ne. 1 ) goto 76

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
	write(6,*) 'dtime,dtvar,iv,ivar,nvar: ',dtime,dtvar,iv,ivar,nvar
	stop 'error stop shyelab: time mismatch'
	end

!***************************************************************
!***************************************************************
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

