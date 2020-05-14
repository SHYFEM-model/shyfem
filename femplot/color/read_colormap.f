
!--------------------------------------------------------------------------
!
!    Copyright (C) 2016,2019  Georg Umgiesser
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

!******************************************************************

! revision log :
!
! 17.06.2016	ggu	changed VERS_7_5_15
! 14.02.2019	ggu	changed VERS_7_5_56

!******************************************************************


	subroutine read_maps(file,nmaps,ncdim,cdata)

	implicit none

	character*(*) file
	integer nmaps,ncdim
	real cdata(3,ncdim,nmaps)

	end

!******************************************************************

	subroutine read_colormap(file,cname,ncdim,cdata,berr)

	implicit none

	character*(*) file
	character*(*) cname	!check if empty, read if set
	integer ncdim
	real cdata(3,ncdim)
	logical berr

	logical bend,bread,bfound
	logical bincolormap
	integer iend,nc,nmap,ios
	real col(3)
	character*80 line,name

	nmap = 0
	nc = 0
	bincolormap = .false.		!inside colormap
	bread = .false.			!have to read data?
	bfound = .false.		!have found right colormap?
	berr = .true.			!have encountered error?
	iend = 0
	bend = .false.

	cdata = -1.

	open(1,file=file,status='old',form='formatted',iostat=ios)
	if( ios /= 0 ) then
	  write(6,*) 'cannot read colormap ',trim(file)
	  berr = .true.
	  return
	end if

	do 
	  read(1,'(a)',iostat=ios) line
	  if( ios /= 0 ) exit
	  if( line == ' ' ) cycle
	  line = adjustl(line)
	  if( line(1:1) == '#' ) cycle
	  if( line(1:1) == '_' ) then
	    bincolormap = .true.
	    nmap = nmap + 1
	    nc = 0
	    iend = scan(line(2:),'_')
	    name = line(2:iend)
	    if( name == cname ) then
	      bread = .true.
	      bfound = .true.
	    end if
	    iend = scan(line,'[') + 1
	    line = adjustl(line(iend:))
	    write(6,*) 'start of colormap: ',trim(name),'  ',trim(line)
	    if( bread ) write(6,*) '...inserting'
	  end if
	  if( bincolormap ) then
	    bend = .false.
	    iend = len_trim(line)
	    if( line(iend-1:iend) == ']]' ) then
	      bend = .true.
	    else if( line(iend-1:iend) == '],' ) then
	      !still in colormap
	    else
	      write(6,*) 'error reading colormap ',trim(line)
	    end if
	  end if
	  line = line(2:iend-2)
	  !write(6,*) 'scanning line: ',trim(line)
	  read(line,*) col
	  !write(6,*) 'scanned line: ',col
	  nc = nc + 1
	  if( bread ) then
	    if( nc > ncdim ) then
	      write(6,*) 'colortable too big: ',ncdim
	      berr = .true.
	      return
	    end if
	    cdata(:,nc) = col
	  end if

	  if( bend ) then
	    bincolormap = .false.
	    bread = .false.
	    write(6,*) 'end of colormap: ',trim(name),'  ',nc
	  end if
	end do

	if( bincolormap ) then
	  write(6,*) 'no end of colormap ',trim(file)
	  berr = .true.
	  return
	end if

	if( ios > 0 ) then
	  write(6,*) 'error reading colormap ',trim(file)
	  berr = .true.
	  return
	end if

	close(1)

	berr = .not. bfound

	end

!******************************************************************

	subroutine test_colormaps

	implicit none

	logical berr
	integer, parameter :: ncdim = 256
	real cdata(3,ncdim)
	character*80 file,cname

	file='colormap.dat'
	cname='inferno'
	cname=' '

	call read_colormap(file,cname,ncdim,cdata,berr)

	if( berr ) then
	  stop 'error stop: reading color table'
	end if

	end

!******************************************************************

	program test_colormaps_main
	call test_colormaps
	end

!******************************************************************

