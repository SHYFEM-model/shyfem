
!--------------------------------------------------------------------------
!
!    Copyright (C) 2018-2020  Georg Umgiesser
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

! routines to run fem_check for fem files

! revision log :
!
! 22.02.2018	ggu	changed VERS_7_5_42
! 03.04.2018	ggu	changed VERS_7_5_43
! 06.07.2018	ggu	changed VERS_7_5_48
! 16.10.2018	ggu	changed VERS_7_5_50
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 23.09.2019	ggu	in fem_check handle case when atime==-1
! 17.04.2020	ggu	better handling of directional variables
! 18.05.2020	ggu	period week implemented

!**************************************************************************

c*****************************************************************
c*****************************************************************
c*****************************************************************
c below files for -check
c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine fem_check(atime,np,lmax,nvar,data,flag
     +				,strings,scheck,bquiet)

	use iso8601

	implicit none

	double precision atime
	integer np,lmax,nvar
	real data(lmax,np,nvar)
	real flag
	character*(*) strings(nvar)
	character*(*) scheck		!period on which to average
	logical bquiet

	logical bwrite,bw,bfile
	integer date,time
	integer iv,nacu,i,l,ivar,ivarss,ivardd,isub,iw,idw
	real data_profile(lmax)
	double precision aver(nvar)
	double precision dtime,dtot
	real :: val,vals,vald
	logical bvel,bwind
	character*20 dline
	character*80 varnum,filename,file,aux,range,string
	character*80 dir,unit
	integer, save :: dt(8),dt0(8)
	integer, save :: iu = 0
	integer, allocatable, save :: ius(:)
	integer, allocatable, save :: ivars(:)
	character*80, allocatable, save :: filenames(:)
	character*80, allocatable, save :: strings_out(:)
	integer, save :: idt
	integer, save :: naccum,ivect
	logical, save :: bfirst = .true.
	logical, save :: bmeteo = .false.
	double precision, save :: atime0,aatime,atimelast
	double precision, allocatable, save :: accum(:)
	double precision, allocatable, save :: astd(:)
	double precision, allocatable, save :: amin(:)
	double precision, allocatable, save :: amax(:)
	double precision, allocatable, save :: facts(:)
	double precision, parameter :: high = 1.e+30
	double precision, parameter :: dlim = 0.9  !fraction of period needed

	real, parameter :: fact(8) = (/365.25,30.5,1.,0.,0.,0.,7.,0./)

	logical string_is_this_short
	integer ifileo

	if( scheck == ' ' ) return

!	-------------------------------
!	initialize
!	-------------------------------

	bw = .not. bquiet

	if( iu == 0 ) then
	  !iu = ifileo(88,'out.txt','form','new')
	  atime0 = atime
	  atimelast = atime
	  call dts_from_abs_time(date,time,atime)
	  call datetime2dt((/date,time/),dt0)
	  call week_of_year(atime,dt0(7))	!for week - memorize at idt=7
	  allocate(accum(nvar))
	  allocate(astd(nvar))
	  allocate(amin(nvar))
	  allocate(amax(nvar))
	  allocate(facts(nvar))
	  allocate(ivars(nvar))
	  allocate(ius(nvar))
	  allocate(filenames(nvar))
	  allocate(strings_out(nvar))
	  naccum = 0
	  accum = 0.
	  astd = 0.
	  amin = high
	  amax = -high
	  facts = 0.
	  idt = 0				!compute total
	  ivect = 0
	  aatime = 0.
	  string = "all year month day week none"
	  i = index(string,trim(scheck))
	  if( i == 0 ) then
	    write(6,*) 'period for check not recognized: ',trim(scheck)
	    write(6,*) 'possible periods: ',trim(string)
	    stop 'error stop fem_check: no such period'
	  end if
	  if( scheck == 'all' ) idt = 0		!compute for whole period
	  if( scheck == 'year' ) idt = 1	!compute on year
	  if( scheck == 'month' ) idt = 2	!compute on month
	  if( scheck == 'day' ) idt = 3		!compute on day
	  if( scheck == 'week' ) idt = 7	!compute on day
	  if( scheck == 'none' ) idt = -1	!output every time step
	  do iv=1,nvar
	    string = strings(iv)

	    if( string_is_this_short('rain',string) ) then
	      if( idt > 0 ) then
	        if( bw ) write(6,*) 'setting facts: ',idt,fact(idt)
	        facts(iv) = fact(idt)
	      end if
	    end if

            call string2ivar(string,ivar)

	    call get_direction_ivars(ivar,ivarss,ivardd)
	    if( ivarss > 0 ) then	!directional
	      call strings_meteo_convention(ivar,bmeteo)
              call string_direction_and_unit(string,dir,unit)
              if( dir == 'x' ) then
  	        ivar = ivarss
		ivect = iv
              else if( dir == 'y' ) then
  	        ivar = ivardd
		if( iv /= ivect+1 ) then
		  stop 'error stop fem_check: internal error (1)'
		end if
              else
                write(6,*) 'unknown direction: ',trim(string),'  ',dir
                stop 'error stop fem_check: unknown direction'
	      end if
            end if

            call ivar2filename(ivar,filename)
	    call ivar2string(ivar,string,isub)
            call strings_pop_direction(filename)
            file = 'aver.' // trim(filename) // '.txt'
	    call get_new_unit(iu)
            open(iu,file=file,form='formatted',status='unknown')
	    write(iu,'(a)') '#      date_and_time    minimum'//
     +			'       average       maximum       std'
	    filenames(iv) = file
	    ius(iv) = iu
	    ivars(iv) = ivar
	    strings_out(iv) = string
	  end do

	  if( ivect > 0 .and. .not. bquiet ) then
	    write(6,*) 'file contains directional data...'
	    write(6,*) 'average done on the following variables:'
            write(6,*) '   varnum     varid    varname'
	    do iv=1,nvar
	      ivar = ivars(iv)
	      string = strings_out(iv)
              write(6,'(2i10,4x,a)') iv,ivar,trim(string)
	    end do
	  end if
	end if

!	-------------------------------
!	elaborate time
!	-------------------------------

	call dts_from_abs_time(date,time,atime)
	call datetime2dt((/date,time/),dt)

	if( idt == 7 ) then	!handle week
	  call week_of_year(atime,dt(idt))
	  if( dt(idt) < dt0(idt) ) then		!new year
	    call weekday(atime,idw)
	    if( idw /= 1 ) dt0(idt) = dt(idt)	!still same week, not Monday
	  end if
	end if

!	-------------------------------
!	average spatially
!	-------------------------------

	do iv=1,nvar
	  if( iv == ivect ) then
	    call aver_vect_data(np,lmax,data(:,:,iv:iv+1)
     +				,flag,bmeteo,vals,vald)
	    aver(iv) = vals
	  else if( ivect > 0 .and. iv == ivect+1 ) then
	    aver(iv) = vald
	  else
	    call aver_data(np,lmax,data(:,:,iv),flag,val)
	    aver(iv) = val
	  end if
	end do

!	-------------------------------
!	check if we have to write
!	-------------------------------

	bwrite = .true.
	if( atime == -2 ) then		!last message
	  bwrite = .false.
	else if( idt == -1 ) then	!output every time step
	  accum = aver
	  amin = aver
	  amax = aver
	  astd = 0
	  aatime = atime
	else if( idt > 0 .and. dt(idt) /= dt0(idt) ) then !over period (FIXME)
	  accum = accum / naccum
	  astd = sqrt( astd/naccum - accum*accum )
	  aatime = atime0 + aatime / naccum
	else if( atime == -1 .and. naccum > 0 ) then	!last time step
	  accum = accum / naccum
	  astd = sqrt( astd/naccum - accum*accum )
	  aatime = atime0 + aatime / naccum
	else
	  bwrite = .false.
	end if

!	-------------------------------
!	write results
!	-------------------------------

	if( bwrite ) then
	  bfile = .true.
	  if( idt > 0 ) then
	    dtime = atimelast - atime0
	    dtime = atime - atime0
	    dtot = fact(idt) * 86400
	    if( dtime/dtot < dlim ) bfile = .false.
	    !write(6,*) idt,naccum,bfile,dtime/dtot,dlim
	  end if

	  where( facts /= 0. )
	    accum = accum * facts
	  end where
	  call dts_format_abs_time(aatime,dline)
	  if( aatime == -1. ) bfile = .false.

	  if( bfirst  .and. bw .and. bfile ) then
	    write(6,*)
	    write(6,'(a)') 'varid naccum         date_and_time'//
     +			'    minimum       average       maximum'
	    bfirst = .false.
	  end if

	  do iv=1,nvar
	    if( bw .and. bfile ) then
	      ivar = ivars(iv)
	      write(6,1000) ivar,naccum,dline,amin(iv),accum(iv),amax(iv)
	    end if
 1000	    format(i5,i7,2x,a20,2x,3e14.6)
	    iu = ius(iv)
	    if( bfile ) then
	      write(iu,1010) dline,amin(iv),accum(iv),amax(iv),astd(iv)
	    end if
 1010	    format(a20,2x,4e14.6)
	  end do
	  dt0 = dt
	  naccum = 0.
	  accum = 0.
	  astd = 0.
	  amin = high
	  amax = -high
	  aatime = 0.
	  atime0 = atime
	end if

!	-------------------------------
!	accumulate in time
!	-------------------------------

	naccum = naccum + 1
	accum = accum + aver
	astd = astd + aver**2
	aatime = aatime + (atime-atime0)
	atimelast = atime
	do iv=1,nvar
	  amin(iv) = min(amin(iv),aver(iv))
	  amax(iv) = max(amax(iv),aver(iv))
	end do

!	-------------------------------
!	final message
!	-------------------------------

	if( bw .and. atime == -2. ) then
	  call compute_range(nvar,range)
	  write(6,*) 'output written to following files:'
	  do iv=1,nvar
	    write(6,*) '  ',trim(filenames(iv))
	  end do
	  write(6,*) 'the four colums are min/aver/max/std'
	  write(6,*) 'the averaging has been done over period: '
     +				,trim(scheck)
	end if

!	-------------------------------
!	end of routine
!	-------------------------------

	end 

c*****************************************************************

	subroutine aver_data(np,lmax,data,flag,aver)

	implicit none

	integer np,lmax
	real data(lmax,np)
	real flag
	real aver

	logical, parameter :: bmax = .false.
	!logical, parameter :: bmax = .true.
	double precision, parameter :: high = 1.e+30
	integer nacu,l,i
	real val,rmax
	double precision acu

	nacu = 0
	acu = 0.
	rmax = -high

	do i=1,np
	  do l=1,lmax
	    val = data(l,i)
	    if( val /= flag ) then
	      nacu = nacu + 1
	      acu = acu + val
	      rmax = max(rmax,val)
	    end if
	  end do
	end do

	if( nacu == 0 ) then
	  aver = flag
	else if( bmax ) then
	  aver = rmax
	else
	  aver = acu / nacu
	end if

	end

c*****************************************************************

	subroutine aver_vect_data(np,lmax,data,flag,bmeteo,avers,averd)

	implicit none

	integer np,lmax
	real data(lmax,np,2)
	real flag
	logical bmeteo
	real avers,averd

	logical, parameter :: bmax = .false.
	!logical, parameter :: bmax = .true.
	double precision, parameter :: high = 1.e+30
	integer nacu,l,i
	real valx,valy,vals,vald,val,rmaxs,rmaxd
	double precision acus,acud

	nacu = 0
	acus = 0.
	acud = 0.
	rmaxs = -high
	rmaxd = -high

	do i=1,np
	  do l=1,lmax
	    valx = data(l,i,1)
	    valy = data(l,i,2)
	    if( valx /= flag .and. valy /= flag ) then
	      nacu = nacu + 1
	      call convert_uv_sd(valx,valy,vals,vald,bmeteo)
	      acus = acus + vals
	      acud = acud + vald
	      rmaxs = max(rmaxs,vals)
	      rmaxd = max(rmaxd,vald)
	    end if
	  end do
	end do

	if( nacu == 0 ) then
	  avers = flag
	  averd = flag
	else if( bmax ) then
	  avers = rmaxs
	  averd = rmaxd
	else
	  avers = acus / nacu
	  averd = acud / nacu
	end if

	end

c*****************************************************************

