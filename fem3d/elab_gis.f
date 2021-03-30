
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015-2016,2018-2019  Georg Umgiesser
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

! utility routines for shyelab: elabutil
!
! revision log :
!
! 15.07.2015	ggu	written from scratch
! 22.09.2015	ggu	new routine open_shy_file()
! 10.10.2015	ggu	code added to handle FLX routines
! 22.02.2016	ggu	handle catmode
! 15.04.2016	ggu	handle gis files with substitution of colon
! 30.05.2016	ggu	changed VERS_7_5_11
! 22.02.2018	ggu	changed VERS_7_5_42
! 16.02.2019	ggu	changed VERS_7_5_60
! 30.03.2021	ggu	better output documentation
!
!************************************************************


c***************************************************************
c***************************************************************
c***************************************************************

        subroutine gis_write_record(dtime,ivar,np,nlvddi,il,cv,xv,yv)

c writes one record to file (3D)

        !use basin
        use shyfem_strings

        implicit none

	double precision dtime
        integer ivar,np,nlvddi
        integer il(np)
        real cv(nlvddi,np)
	real xv(np),yv(np)

	logical bold
	integer it
        integer i,l,lmax
	integer nout
        real x,y
	character*80 format,name,short,full
	character*20 line,dateline
	character*3 var

	integer ifileo

	bold = .false.		!old or new format

	it = nint(dtime)
	call dtsgf(it,dateline)
	call gis_subst_colon(dateline,line)
	call strings_get_short_name(ivar,short)
	call strings_get_full_name(ivar,full)

	name = 'extract_'//trim(short)//'_'//line//'.gis'
        nout = ifileo(60,name,'form','new')
	!write(6,*) 'writing: ',trim(name)

	if( bold ) then
          write(nout,*) it,np,ivar,dateline
	else
          write(nout,1000) '# data extracted to gis format by shyelab'
          write(nout,1001) '# date:    ',dateline
          write(nout,1002) '# ivar:    ',ivar
          write(nout,1001) '# varname: ',trim(full)
          write(nout,1002) '# nodes:   ',np
          write(nout,1000) '#    inode'//
     +			'             x             y'//
     +			'  layers      data (all layers)'
	end if

	lmax = 1

        do i=1,np
          if( nlvddi > 1 ) lmax = il(i)
          x = xv(i)
          y = yv(i)
	  write(format,'(a,i5,a)') '(i10,2e14.6,i8,',lmax,'g14.6)'
          write(nout,format) i,x,y,lmax,(cv(l,i),l=1,lmax)
        end do

	close(nout)

	return
 1000	format(a)
 1001	format(a,a)
 1002	format(a,i10)
        end

c***************************************************************

        subroutine gis_write_hydro(dtime,np,nlvddi,il,zv,uv,vv,xv,yv)

c writes one record to file (3D)

        !use basin

        implicit none

	double precision dtime
	integer np
        integer nlvddi
        integer il(np)
        real zv(np)
        real uv(nlvddi,np)
        real vv(nlvddi,np)
	real xv(np)
	real yv(np)

	logical bold
        integer l,lmax,nn,i,it
	integer nout
	integer ivar
        real x,y
	character*80 format,name,short,full
	character*20 line,dateline
	character*3 var

	integer ifileo

	bold = .false.		!old or new format

	it = nint(dtime)
	call dtsgf(it,dateline)
	call gis_subst_colon(dateline,line)

	name = 'extract_hydro_'//line//'.gis'
        nout = ifileo(60,name,'form','new')
	!write(6,*) 'writing: ',trim(name)

	if( bold ) then
          write(nout,*) it,np,0,dateline
	else
	  write(nout,1000) '# data extracted to gis format by shyelab'
	  write(nout,1001) '# date:    ',dateline
	  write(nout,1002) '# ivar:    ',0
	  write(nout,1001) '# varname: '
     +			,'water level (zeta) and velocities'
	  write(nout,1002) '# nodes:   ',np
	  write(nout,1000) '#   inode              x             y'//
     +				'  layers      zeta'//
     +				'     velocities x/y (all layers)'
	end if

	lmax = 1

        do i=1,np
          if( nlvddi > 1 ) lmax = il(i)
          x = xv(i)
          y = yv(i)
	  nn = 1 + 2*lmax
	  write(format,'(a,i5,a)') '(i10,2e14.6,i8,',nn,'g14.6)'
          write(nout,format) i,x,y,lmax,zv(i)
     +			,(uv(l,i),vv(l,i),l=1,lmax)
        end do

	close(nout)

	return
 1000	format(a)
 1001	format(a,a)
 1002	format(a,i10)
        end

c***************************************************************

	subroutine gis_subst_colon(line_old,line_new)

	implicit none

	character*(*) line_old,line_new

	integer n,i

	n = min(len(line_old),len(line_new))
	line_new = line_old

	do i=1,n
	  if( line_new(i:i) == ':' ) line_new(i:i) = '_'
	  if( line_new(i:i) == ' ' ) line_new(i:i) = '_'
	end do

	end

c***************************************************************

        subroutine gis_write_connect

c writes connectivity

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

c***************************************************************
