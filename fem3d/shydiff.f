
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

! revision log :
!
! 18.11.1998	ggu	check dimensions with dimnos
! 06.04.1999	ggu	some cosmetic changes
! 03.12.2001	ggu	some extra output -> place of min/max
! 09.12.2003	ggu	check for NaN introduced
! 07.03.2007	ggu	easier call
! 08.11.2008	ggu	do not compute min/max in non-existing layers
! 07.12.2010	ggu	write statistics on depth distribution (depth_stats)
! 06.05.2015	ggu	noselab started
! 05.06.2015	ggu	many more features added
! 10.09.2015	ggu	std and rms for averaging implemented
! 11.09.2015	ggu	write in gis format
! 23.09.2015	ggu	handle more than one file (look for itstart)
! 16.10.2015	ggu	started shyelab
! 28.04.2016	ggu	changed VERS_7_5_9
! 25.05.2016	ggu	changed VERS_7_5_10
! 30.05.2016	ggu	changed VERS_7_5_11
! 05.10.2016	ggu	changed VERS_7_5_19
! 04.11.2017	ggu	changed VERS_7_5_34
! 24.01.2018	ggu	changed VERS_7_5_41
! 16.02.2019	ggu	changed VERS_7_5_60
!
!**************************************************************

	program shydiff

	use clo
	use elabutil
	use elabtime
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

	real, allocatable :: cv3(:,:)
	real, allocatable :: cv4(:,:)
	real, allocatable :: cv3all(:,:,:)
	real, allocatable :: cv4all(:,:,:)

	integer, allocatable :: ivars(:,:)
	integer, allocatable :: il(:)

	real, allocatable :: hl(:)

	logical bnextfile
	integer nwrite,nread,nelab,nrec,nin,nold,ndiff
	integer nvers
	integer nexit
	integer nknnos,nelnos,nvar,npr
	integer ierr
	!integer it,itvar,itnew,itold,itstart
	integer it
	integer ivar,iaux
	integer iv,j,l,k,lmax,node
	integer ip,nb
	integer ifile,ftype
	integer id,idout,idold,id2
	integer n,m,nndim,nn
	integer naccum
	integer date,time
	character*80 title,name,file
	character*80 basnam,simnam
	real rnull
	real cmin,cmax,cmed,vtot
	double precision dtime,dtime2,dtstart,dtnew,dtvar
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
	ndiff = 0
	rnull=0.
	rnull=-1.
	bopen = .false.
	bzeta = .false.		!file has zeta information

	!--------------------------------------------------------------
	! set command line parameters
	!--------------------------------------------------------------

	call elabutil_init('SHY','shydiff')

	!--------------------------------------------------------------
	! open input files
	!--------------------------------------------------------------

	ifile = 1
	call open_next_file(ifile,0,id)
	ifile = 2
        call clo_get_file(ifile,file)
        if( file == ' ' ) then
	  write(6,*) 'need two files for shydiff'
	  stop 'error stop shydiff: two files needed'
	end if
	call open_next_file(ifile,id,id2)

	!--------------------------------------------------------------
	! set up params and arrays
	!--------------------------------------------------------------

	call shy_get_params(id,nkn,nel,npr,nlv,nvar)
	call shy_get_ftype(id,ftype)

	call shy_info(id)

        call basin_init(nkn,nel)
        call levels_init(nkn,nel,nlv)
        call mod_depth_init(nkn,nel)
	call shy_copy_basin_from_shy(id)
	call shy_copy_levels_from_shy(id)
        call ev_init(nel)
	call set_ev

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

	allocate(cv3(nlv,nndim))
	allocate(cv4(nlv,nndim))
	allocate(cv3all(nlv,nndim,0:nvar))
	allocate(cv4all(nlv,nndim,0:nvar))

        allocate(hl(nlv))
	allocate(ivars(4,nvar))

	call shyutil_init(nkn,nel,nlv)

	call init_sigma_info(nlv,hlv)

	call shy_make_area
	call outfile_make_depth(nkn,nel,nen3v,hm3v,hev,hkv)

	if( bverb ) call depth_stats(nkn,nlvdi,ilhkv)

	!--------------------------------------------------------------
	! time management
	!--------------------------------------------------------------

	call shy_get_date(id,date,time)
	call elabtime_date_and_time(date,time)
        call elabtime_set_minmax(stmin,stmax)
        call elabtime_set_inclusive(.false.)

!--------------------------------------------------------------
! loop on data
!--------------------------------------------------------------

	dtvar = 0.
	dtime = 0.
	call shy_peek_record(id,dtime,iaux,iaux,iaux,iaux,ierr)
	call dts_convert_to_atime(datetime_elab,dtime,atime)

	cv3 = 0.
	cv4 = 0.
	cv3all = 0.
	cv4all = 0.

	do

	 atold = atime

	 call read_records(id,dtime,nvar,nndim,nlvdi,ivars
     +				,cv3,cv3all,ierr)
         if(ierr.ne.0) exit

	 call read_records(id2,dtime2,nvar,nndim,nlvdi,ivars
     +				,cv4,cv4all,ierr)
         if(ierr.ne.0) exit

	 !if( dtime /= dtime2 ) then
	 !  ndiff = ndiff + 1
	 !  write(6,*) '*** time is different: ',dtime,dtime2
	 !  exit
	 !end if

	 nread = nread + nvar
	 nrec = nrec + 1
	 call dts_convert_to_atime(datetime_elab,dtime,atime)

	 call shy_peek_record(id,dtnew,iaux,iaux,iaux,iaux,ierr)
	 call dts_convert_to_atime(datetime_elab,dtnew,atnew)
	 if( ierr .ne. 0 ) atnew = atime

	 if( elabtime_over_time(atime,atnew,atold) ) exit
	 if( .not. elabtime_in_time(atime,atnew,atold) ) cycle

	 do iv=1,nvar

	  ivar = ivars(4,iv)
	  lmax = ivars(3,iv)
	  nn = ivars(1,iv) * ivars(2,iv)

	  nelab=nelab+1

	  if( .not. bquiet ) then
	    call shy_write_time(.true.,dtime,atime,ivar)
	  end if

	  if( bwrite ) then
	    call shy_write_min_max(nlvdi,nn,il,lmax,cv3)
	  end if

	  cv3(:,:) = cv3all(:,:,iv)
	  cv4(:,:) = cv4all(:,:,iv)

	  call diff_shy(nlv,nn,cv3,cv4,ndiff)

	 end do		!loop on ivar
	end do		!time do loop

!--------------------------------------------------------------
! end of loop on data
!--------------------------------------------------------------

!--------------------------------------------------------------
! final write of variables
!--------------------------------------------------------------

!--------------------------------------------------------------
! write final message
!--------------------------------------------------------------

	write(6,*)
	write(6,*) nread, ' records read'
	write(6,*) nrec , ' unique time records read'
	write(6,*) nelab, ' records elaborated'
	write(6,*) ifile, ' files read'
	write(6,*) ndiff, ' differing records found'
	write(6,*)
	if( ndiff == 0 ) then
	  nexit = 99
	  write(6,*) 'files are equal'
	else
	  nexit = 0
	  write(6,*) '*** files are different'
	end if
	write(6,*)
	call exit(nexit)

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

	subroutine diff_shy(nlvddi,nn,cv3,cv4,ndiff)

	implicit none

	integer nlvddi,nn
	real cv3(nlvddi,nn)
	real cv4(nlvddi,nn)
	integer ndiff

	real, parameter :: eps = 1.e-6
	!real, parameter :: eps = 0.
	real diff
	integer iloc(2)

	diff = maxval( abs(cv3-cv4) )

	if( diff > eps ) then
	  ndiff = ndiff + 1
	  iloc = maxloc( abs(cv3-cv4) )
	  write(6,*) '*** record differing: ',diff,iloc
	end if

	end

!***************************************************************

