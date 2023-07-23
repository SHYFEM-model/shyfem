
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998-1999,2001,2003,2007-2008,2010  Georg Umgiesser
!    Copyright (C) 2015-2020  Georg Umgiesser
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

c revision log :
c
c 18.11.1998	ggu	check dimensions with dimnos
c 06.04.1999	ggu	some cosmetic changes
c 03.12.2001	ggu	some extra output -> place of min/max
c 09.12.2003	ggu	check for NaN introduced
c 07.03.2007	ggu	easier call
c 08.11.2008	ggu	do not compute min/max in non-existing layers
c 07.12.2010	ggu	write statistics on depth distribution (depth_stats)
c 06.05.2015	ggu	noselab started
c 05.06.2015	ggu	many more features added
c 14.09.2015	ggu	changed VERS_7_2_2
c 05.10.2015	ggu	variables in xv were used in the wromg order - fixed
c 12.10.2015	ggu	changed VERS_7_3_3
c 19.10.2015	ggu	changed VERS_7_3_6
c 25.05.2016	ggu	changed VERS_7_5_10
c 05.10.2017	ggu	implement quiet, silent option, write dir
c 09.10.2017	ggu	consistent treatment of output files
c 20.10.2017	ggu	write time in string format
c 26.10.2017	ggu	various user related improvements
c 04.11.2017	ggu	changed VERS_7_5_34
c 14.11.2017	ggu	changed VERS_7_5_36
c 05.12.2017	ggu	changed VERS_7_5_39
c 22.02.2018	ggu	changed VERS_7_5_42
c 16.10.2018	ggu	changed VERS_7_5_50
c 16.02.2019	ggu	changed VERS_7_5_60
c 03.02.2020	ggu	cleaned, new headers for split
c 06.03.2020	ggu	check for time step
c
c**************************************************************

	subroutine extelab

	use clo
	use elabutil
	use elabtime
	use shyfem_strings

        use basin
        use mod_depth
        use evgeom
        use levels

c elaborates ext file

	implicit none

	integer, allocatable :: knaus(:)
	integer, allocatable :: il(:)
	real, allocatable :: hdep(:)
	real, allocatable :: x(:)
	real, allocatable :: y(:)
	real, allocatable :: hl(:)
	character*80, allocatable :: strings(:)
	integer, allocatable :: ivars(:)
	real, allocatable :: xv(:,:)
	real, allocatable :: vals(:,:,:)
	real, allocatable :: val(:,:)
	real, allocatable :: val0(:)
	real, allocatable :: zeta(:)
	real, allocatable :: uu(:)
	real, allocatable :: vv(:)
	real, allocatable :: speed(:)
	real, allocatable :: dir(:)

	integer, allocatable :: naccu(:)
	double precision, allocatable :: accum(:,:,:)

	logical btskip
	integer nread,nelab,nrec,nout,nin,nn
	integer nvers
	integer nknnos,nelnos,nvar
	integer ierr
	integer it,ivar,ivarn,iv,itvar,iaux
	integer ii,i,j,l,k,lmax,node,m
	integer ip,nb,naccum
	integer knausm
	integer date,time
	character*80 title,name,file,femver,format,range
	character*20 dline
	character*80 basnam,simnam
	real rnull
	real cmin,cmax,cmed,vtot
	real vmin,vmax,vmed
	real href,hzmin,nzadapt,rzmov
	real s,d
	double precision atime,atfirst,atlast,atold,atnew,atwrite,atime0
	character*10 :: short
	character*40 :: full
	integer, parameter :: niu = 6
	character*5 :: what(niu) = (/'velx ','vely ','zeta '
     +				,'speed','dir  ','all  '/)
	character*26 :: descrp(niu) = (/
     +		 'velocity in x-direction   '
     +		,'velocity in y-direction   '
     +		,'water level               '
     +		,'current speed             '
     +		,'current direction         '
     +		,'all hydrodynamic variables'
     +			/)

	integer iapini,ifileo
	integer ifem_open_file

c--------------------------------------------------------------

	nread=0
	nelab=0
	nrec=0
	nout=0
	rnull=0.
	rnull=-1.
	bopen = .false.

c--------------------------------------------------------------
c set command line parameters
c--------------------------------------------------------------

	call elabutil_init('EXT','extelab')

	!--------------------------------------------------------------
	! open input files
	!--------------------------------------------------------------

	call clo_reset_files
	call clo_get_next_file(file)
	nin = ifileo(0,file,'unform','old')

        call populate_strings

        call ext_is_ext_file(nin,nvers)
        if( nvers .le. 0 ) then
          write(6,*) 'nvers: ',nvers
          stop 'error stop extelab: not a valid ext file'
        end if

	call ext_read_header(nin,nvers,knausm,lmax,nvar,ierr)
	if( ierr /= 0 ) goto 93
	nlv = 1

	allocate(knaus(knausm))
	allocate(hdep(knausm))
	allocate(il(knausm))
	allocate(x(knausm))
	allocate(y(knausm))
	allocate(strings(knausm))
	allocate(hl(lmax))
	allocate(xv(knausm,3))
	allocate(vals(lmax,knausm,3))
	allocate(val(lmax,knausm))
	allocate(val0(knausm))
	allocate(zeta(knausm))
	allocate(uu(knausm))
	allocate(vv(knausm))
	allocate(speed(knausm))
	allocate(dir(knausm))
	allocate(ivars(nvar))
	femver = ' '
	ivars = 0

	call ext_read_header2(nin,nvers,knausm,lmax
     +                          ,atime0
     +                          ,href,hzmin,nzadapt,title,femver
     +                          ,knaus,hdep,il,x,y,strings,hl
     +				,ierr)
	if( ierr /= 0 ) goto 93

	if( .not. bquiet ) then
          write(6,*) 'file       : ',trim(file)
          write(6,*) 'nvers      : ',nvers
          write(6,*) 'knausm     : ',knausm
          write(6,*) 'lmax       : ',lmax
          write(6,*) 'nvar       : ',nvar
          write(6,*) 'href,hzmin : ',href,hzmin
	  write(6,*) 'nzadapt    : ',nint(nzadapt)
          write(6,*) 'title      : ',trim(title)
	  write(6,*) 'Nodes contained in file:'
          write(6,*) ' i    node  il      hdep' //
     +			'           x           y  description'
	  do i=1,knausm
            write(6,1000) i,knaus(i),il(i),hdep(i),x(i),y(i)
     +					,'  ',trim(strings(i))
 1000	    format(i3,i8,i4,f10.2,2f12.2,a,a)
	  end do
	end if

	!--------------------------------------------------------------
	! time management
	!--------------------------------------------------------------

	date = 0
	time = 0
	call elabtime_date_and_time(date,time)	!we work with absolute time
	call elabtime_set_minmax(stmin,stmax)

	call ext_peek_record(nin,nvers,atime,ivar,ierr)
	if( ierr /= 0 ) goto 91
	atfirst = atime
	atold = atime

	!--------------------------------------------------------------
	! see what is in the file
	!--------------------------------------------------------------

	if( .not. bquiet ) then
	  write(6,*) 'Variables contained in file:'
	  write(6,*) ' i ivar  short     full'
	  do iv=1,nvar
	    call ext_read_record(nin,nvers,atime,knausm,lmax
     +				,ivar,m,il,vals
     +				,ierr)
	    if( ierr /= 0 ) goto 91
	    if( ivar == 0 ) ivar = 1
	    ivars(iv) = ivar
	    call strings_get_short_name(ivar,short)
	    call strings_get_full_name(ivar,full)
	    write(6,'(i3,i5,a,a,a)') iv,ivar,'  ',short,full
	  end do

	  do iv=1,nvar
	    backspace(nin)
	  end do
	end if

	if( binfo ) return

	!--------------------------------------------------------------
	! averaging
	!--------------------------------------------------------------

	nlv = lmax
	call init_sigma_info(nlv,hl)
	call init_rzmov_info(nlv,nint(nzadapt),hl,rzmov)
	!call elabutil_set_averaging(nvar)

	btrans = .false.
	if( btrans ) then
	  allocate(naccu(istep))
	  allocate(accum(nlvdi,nkn,istep))
	else
          allocate(naccu(1))
          allocate(accum(1,1,1))
        end if
        naccum = 0
        naccu = 0
        accum = 0.

	!write(6,*) 'mode: ',mode,ifreq,istep

	!--------------------------------------------------------------
	! open output file
	!--------------------------------------------------------------

	boutput = boutput .or. btrans
	bopen = boutput
	if( bsplit ) then
	  boutput = .false.
	  bopen = .false.
	end if

	if( bopen ) then
	  nb = ifileo(0,'out.ext','unform','new')
	  call ext_write_header(nb,0,knausm,lmax,nvar,ierr)
	  if( ierr /= 0 ) goto 99
          call ext_write_header2(nb,0,knausm,lmax
     +				,atime0
     +                          ,href,hzmin,nzadapt,title,femver
     +                          ,knaus,hdep,il,x,y,strings,hl
     +                          ,ierr)
	  if( ierr /= 0 ) goto 99
	end if

c--------------------------------------------------------------
c loop on data
c--------------------------------------------------------------

	atime = atold
	atwrite = -1
	if( .not. bquiet ) write(6,*)

	do

	  atold = atime
	  call ext_read_record(nin,nvers,atime,knausm,lmax
     +				,ivar,m,il,vals
     +				,ierr)
          if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
          if(ierr.ne.0) exit
	  nread = nread + 1
	  if( ivar == 1 ) nrec = nrec + 1

	  atlast = atime
	  call ext_peek_record(nin,nvers,atnew,ivarn,ierr)
	  if( ierr .ne. 0 ) atnew = atime

	  if( ivar == 1 ) then
            call handle_timestep(atime,bcheckdt,btskip)
            if( btskip ) cycle
	  end if

          if( elabtime_over_time(atime,atnew,atold) ) exit
          if( .not. elabtime_in_time(atime,atnew,atold) ) cycle
	  !if( .not. elabtime_check_time(atime,atnew,atold) ) cycle

	  if( ivar == 1 ) nelab=nelab+1

	  if( bverb .and. atwrite /= atime ) then
	    call dts_format_abs_time(atime,dline)
	    write(6,*) 'time : ',atime,'  ',dline
	    atwrite = atime
	  end if

	  if( bwrite .or. bsplit ) then
	    if( ivar == 1 ) then
	      uu(:)   = vals(1,:,1)
	      vv(:)   = vals(1,:,2)
	      zeta(:) = vals(1,:,3)
	      do j=1,knausm
	        call c2p_ocean(uu(j),vv(j),speed(j),dir(j))
	      end do
	    else
	      val(:,:) = vals(:,:,1)
	      call average_val(knausm,lmax,il,hl,hdep,zeta,val,val0)
	    end if
	  end if

	  if( bwrite ) then
	    if( ivar == 0 ) then
	      call minmaxmed(zeta,knausm,vmin,vmax,vmed)
	      write(6,*) 'zeta  min,max,aver : ',vmin,vmax,vmed
	      call minmaxmed(uu,knausm,vmin,vmax,vmed)
	      write(6,*) 'velx  min,max,aver : ',vmin,vmax,vmed
	      call minmaxmed(vv,knausm,vmin,vmax,vmed)
	      write(6,*) 'vely  min,max,aver : ',vmin,vmax,vmed
	      call minmaxmed(speed,knausm,vmin,vmax,vmed)
	      write(6,*) 'speed min,max,aver : ',vmin,vmax,vmed
	    else if( ivar == 2 ) then	!velocities... already handled
	      !nothing
	    else
	      call strings_get_short_name(ivar,short)
	      call minmaxmed(val0,knausm,vmin,vmax,vmed)
	      write(6,*) trim(short)//' min,max,aver : ',vmin,vmax,vmed
	    end if
	  end if

	  if( btrans ) then
!	    call nos_time_aver(mode,i,ifreq,istep,nkn,nlvdi
!     +					,naccu,accum,cv3,boutput)
	  end if

	  if( bsplit ) then
	    if( ivar == 1 ) then	!this is always the first record
	      iv = 1
	      call split_var0d(atime,knausm,nvar,what,zeta,uu,vv)
	    else if( ivar == 230 ) then
              call split_var_wave(atime,knausm,lmax,nvar,ivar,iv,vals)
	    else
	      iv = iv + 1
	      if( lmax > 1 ) then
                call split_var3d(atime,knausm,lmax,il
     +				,m,nvar,ivar,iv,vals)
	      end if
	      if( ivar /= 2 ) then	!only if not velocity
                call split_var2d(atime,knausm
     +				,m,nvar,ivar,iv,val0)
	      end if
	    end if
	  end if

	  if( boutput ) then
	    if( ivar == 0 ) nout = nout + 1
	    if( bverb ) write(6,*) 'writing to output: ',ivar,atime
            call ext_write_record(nb,0,atime,knausm,lmax
     +                                  ,ivar,m,il,vals,ierr)
	    if( ierr /= 0 ) goto 99
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
	    naccum = naccu(ip)
	    !write(6,*) 'naccum: ',naccum
	    if( naccum > 0 ) then
	      write(6,*) 'final aver: ',ip,naccum
!	      call nos_time_aver(-mode,ip,ifreq,istep,nkn,nlvdi
!     +					,naccu,accum,cv3,boutput)
!	      if( bsumvar ) ivar = 30
!              call nos_write_record(nb,it,ivar,nlvdi,ilhkv,cv3,ierr)
              if( ierr .ne. 0 ) goto 99
	    end if
	  end do
	end if

c--------------------------------------------------------------
c write final message
c--------------------------------------------------------------

	if( .not. bsilent ) then
          write(6,*)
	  call dts_format_abs_time(atfirst,dline)
          write(6,*) 'first time record: ',dline
	  call dts_format_abs_time(atlast,dline)
          write(6,*) 'last time record:  ',dline

	  call handle_timestep_last(bcheckdt)

	  write(6,*)
	  write(6,*) nread,' data records read'
	  write(6,*) nrec ,' time records read'
	  write(6,*) nelab,' time records elaborated'
	  write(6,*) nout ,' time records written to file'
	  write(6,*)
	end if

	if( .not. bquiet ) then
	 if( bsplit ) then
	  write(6,*) 'output written to following files: '
	  write(6,*) '  what.dim.node'
	  write(6,*) 'what is one of the following:'
	  call write_special_vars(niu,what,descrp)	!write hydro variables
	  if( nvar > 2 ) then
	    call write_vars(nvar-2,ivars(3:))		!write rest of variables
	  end if
	  write(6,*) 'dim is 2d or 3d'
	  write(6,*) '  2d for depth averaged variables'
	  write(6,*) '  3d for output at each layer'
	  call compute_range(knausm,range)
	  write(6,1123) ' node is consecutive node numbering: '
     +					,trim(range)
	 else if( boutput ) then
	  write(6,*) 'output written to file out.ext'
	 end if
	end if

 1123	format(a,a)

	!if( .not. bquiet ) then
	! call ap_get_names(basnam,simnam)
	! write(6,*) 'names used: '
	! write(6,*) '  simul: ',trim(simnam)
	!end if

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	return
   85	continue
	write(6,*) 'it,itvar,i,ivar,nvar: ',it,itvar,i,ivar,nvar
	stop 'error stop extelab: time mismatch'
   91	continue
	write(6,*) 'error reading first data record'
	write(6,*) 'maybe the file is empty'
	stop 'error stop extelab: empty record'
   93	continue
	write(6,*) 'error reading header of file'
	write(6,*) 'maybe the file is empty'
	stop 'error stop extelab: empty header'
   99	continue
	write(6,*) 'error writing to file unit: ',nb
	stop 'error stop extelab: write error'
	end

c***************************************************************
c***************************************************************
c***************************************************************

	subroutine split_var0d(atime,knausm,nvar,what,zeta,uu,vv)

! splits first (old format) record which contains zeta, velx, vely

        implicit none

	integer, parameter :: niu = 6

	double precision atime
        integer knausm,nvar
        character*5 what(niu)
        real zeta(knausm)
        real uu(knausm)
        real vv(knausm)

	integer j,ii,iu
	real s,d
        character*80 name
        character*20 dline
	integer, save :: icall = 0
	integer, save, allocatable :: iusplit(:,:)
	character*20, save :: time_date = '#          date_time'
	character*80, save :: title(niu) = (
     +			/'      velx'
     +			,'      vely'
     +			,'      zeta'
     +			,'     speed'
     +			,'       dir'
     +			,'       all'
     +			/)
	character*80, save :: all =	
     +  '      zeta        velx        vely       speed         dir'

	if( icall == 0 ) then
	  allocate(iusplit(niu,knausm))
	  title(6) = all
	  do j=1,knausm
	    do ii=1,niu
	      call make_iunit_name(what(ii),'','2d',j,iu)
	      iusplit(ii,j) = iu
	      write(iu,'(3a)') time_date,'  ',trim(title(ii))
	    end do
	  end do
	end if

	icall = icall + 1

	call dts_format_abs_time(atime,dline)

	do j=1,knausm
	  iu = iusplit(1,j)
	  write(iu,1000) dline,uu(j)
	  iu = iusplit(2,j)
	  write(iu,1000) dline,vv(j)
	  iu = iusplit(3,j)
	  write(iu,1000) dline,zeta(j)
	  call c2p_ocean(uu(j),vv(j),s,d)
	  iu = iusplit(4,j)
	  write(iu,1000) dline,s
	  iu = iusplit(5,j)
	  write(iu,1000) dline,d
	  iu = iusplit(6,j)
	  write(iu,'(a20,5f12.4)') dline,zeta(j),uu(j),vv(j),s,d
	end do

	return
 1000	format(a20,f12.4)
        end

c***************************************************************

        subroutine split_var2d(atime,knausm
     +				,m,nvar,ivar,iv,val0)

! splits 2d scalar records

	use shyfem_strings

        implicit none

	double precision atime
        integer knausm,m,nvar,ivar,iv
	real val0(knausm)

	integer j,ii,iu,it
	integer l,lm
        character*80 name,format
        character*20 filename
        character*14 title
        character*20 dline
	character*10 short
	integer, save :: icall = 0
	integer, save, allocatable :: iusplit(:,:)
	character*20, save :: time_date = '#          date_time'

	if( icall == 0 ) then
	  allocate(iusplit(knausm,nvar))
	  iusplit = 0
	end if

	if( iusplit(1,iv) == 0 ) then
	  call ivar2filename(ivar,filename)
	  title = adjustr(filename(1:14))
	  do j=1,knausm
	    call make_iunit_name(filename,'','2d',j,iu)
	    write(iu,'(2a)') time_date,trim(title)
	    iusplit(j,iv) = iu
	  end do
	end if

	icall = icall + 1

	call dts_format_abs_time(atime,dline)

	if( m == 1 ) then	!regular variable
	  do j=1,knausm
	    iu = iusplit(j,iv)
	    write(iu,1000) dline,val0(j)
	  end do
	else
	  write(6,*) 'm,ivar: ',m,ivar
	  write(6,*) 'no multi-variable only for 2d'
	  stop 'error stop split_var2d: no multi-variable'
	end if

	return
 1000	format(a20,f14.4)
	end

c***************************************************************

        subroutine split_var_wave(atime,knausm,lmax,nvar,ivar,iv,vals)

! splits wave data
!
! to be integrated into split_var2d

	use shyfem_strings

        implicit none

        double precision atime
        integer knausm,lmax,nvar,ivar,iv
        integer il(knausm)
        real vals(lmax,knausm,3)

	integer j,ii,iu,it
	integer l,lm
        character*80 name,format
        character*20 filename
        character*20 dline
	character*10 short
	integer, save :: icall = 0
	integer, save, allocatable :: iusplit(:)

	if( icall == 0 ) then
	  allocate(iusplit(knausm))
	  iusplit = 0
	end if

	if( iusplit(1) == 0 ) then
	  call ivar2filename(ivar,filename)
	  do j=1,knausm
	    call make_iunit_name(filename,'','2d',j,iu)
	    iusplit(j) = iu
	  end do
	end if

	icall = icall + 1

	call dts_format_abs_time(atime,dline)

	do j=1,knausm
	  iu = iusplit(j)
	  write(iu,*) dline,vals(1,j,1),vals(1,j,2),vals(1,j,3)
	end do

	end

c***************************************************************

        subroutine split_var3d(atime,knausm,lmax,il
     +				,m,nvar,ivar,iv,vals)

! splits 3d values
	use shyfem_strings

        implicit none

	integer, parameter :: niu = 4

	double precision atime
        integer knausm,lmax,m,nvar,ivar,iv
	integer il(knausm)
	real vals(lmax,knausm,3)

	integer j,ii,iu,it
	integer l,lm
	real u(lmax),v(lmax),s(lmax),d(lmax)
        character*80 name,format
        character*20 dline
        character*20 filename
	character*10 short
	character*5 :: what(niu) = (/'velx ','vely ','speed','dir  '/)
	character*20, save :: time_date = '#          date_time'
	integer, save, allocatable :: iusplit(:,:,:)
	integer, save :: icall = 0

	if( icall == 0 ) then
	  allocate(iusplit(niu,knausm,nvar))
	  iusplit = 0
	end if

	if( iusplit(1,1,iv) == 0 ) then
	  call ivar2filename(ivar,filename)
	  do j=1,knausm
	    if( ivar == 2 ) then	!velocities
	      do ii=1,niu
	        call make_iunit_name(what(ii),'','3d',j,iu)
	        iusplit(ii,j,iv) = iu
	      end do
	    else
	      call make_iunit_name(filename,'','3d',j,iu)
	      iusplit(1,j,iv) = iu
	    end if
	  end do
	end if

	icall = icall + 1

	call dts_format_abs_time(atime,dline)
	write(format,'(a,i3,a)') '(a20,',lmax,'f8.3)'

	if( m == 1 ) then	!regular variable
	  do j=1,knausm
	    lm = min(il(j),lmax)
	    iu = iusplit(1,j,iv)
	    write(iu,format) dline,(vals(l,j,1),l=1,lm)
	  end do
	else if( ivar == 2 ) then
	  do j=1,knausm
	    lm = min(il(j),lmax)
	    u(1:lm) = vals(1:lm,j,1)
	    v(1:lm) = vals(1:lm,j,2)
	    do l=1,lm
	      call c2p_ocean(u(l),v(l),s(l),d(l))
	    end do
	    iu = iusplit(1,j,iv)
	    write(iu,format) dline,(u(l),l=1,lm)
	    iu = iusplit(2,j,iv)
	    write(iu,format) dline,(v(l),l=1,lm)
	    iu = iusplit(3,j,iv)
	    write(iu,format) dline,(s(l),l=1,lm)
	    iu = iusplit(4,j,iv)
	    write(iu,format) dline,(d(l),l=1,lm)
	  end do
	else
	  write(6,*) 'm,ivar: ',m,ivar
	  write(6,*) 'multi-variable only for velocity'
	  stop 'error stop split_var: no multi-variable'
	end if

        end

c***************************************************************
c***************************************************************
c***************************************************************

	subroutine minmaxmed(val,n,vmin,vmax,vmed)

! computes min/max/med of val array
	
	implicit none

	integer n
	real val(n)
	real vmin,vmax,vmed

	real, parameter :: rflag = -999.

	call mimar(val,n,vmin,vmax,rflag)
        call aver(val,n,vmed,rflag)

	end

c***************************************************************

	subroutine average_val(knausm,lmax,il,hl,hdep,zeta,vals,val0)

	implicit none

	integer knausm,lmax
	integer il(knausm)
	real hl(lmax)
	real hdep(knausm)
	real zeta(knausm)
	real vals(lmax,knausm)
	real val0(knausm)

	integer j,lm
	real z,h,v0
	real v(lmax)

	do j=1,knausm
	  lm = il(j)
	  z = zeta(j)
	  h = hdep(j)
	  v = vals(:,j)
	  call average_vertical_node(lm,hl,z,h,v,v0)
	  val0(j) = v0
	end do

	end

c***************************************************************
