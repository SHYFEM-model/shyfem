
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

c utilities for ETS files
c
c revision log :
c
c 24.01.2014    ggu     copied from nosutil.f
c 30.05.2015    ggu     new code for reading header 1 and 2 indipendently
c
c***************************************************************

c***************************************************************
c***************************************************************
c***************************************************************

	subroutine write_ets_header(iu,ilhkv,hlv,hev)

c other variables are stored internally
c
c must have been initialized with ets_init
c all other variables must have already been stored internally (title,date..)

	implicit none

	integer iu
	integer ilhkv(1)
	real hlv(1)
	real hev(1)

	integer nkn,nel,nlv,nvar
	integer ierr

	stop 'error stop write_ets_header: not ready'

	call ets_get_params(iu,nkn,nlv,nvar)
	call ets_write_header(iu,nkn,nlv,nvar,ierr)
	if( ierr .ne. 0 ) goto 99
	!call ets_write_header2(iu,ilhkv,hlv,hev,ierr)
	if( ierr .ne. 0 ) goto 99

	return
   99	continue
	write(6,*) 'error in writing header of ETS file'
	stop 'error stop write_ets_header: writing header'
	end

c***************************************************************

	subroutine read_ets_header(iu,nknddi,nlvddi,ilets,hlv,hets
     +					,nodes,xg,yg,desc)

c other variables are stored internally

	implicit none

	integer iu
	integer nknddi,nelddi,nlvddi
	integer ilets(nknddi)
	real hlv(nlvddi)
	real hets(nknddi)
	integer nodes(nknddi)
	real xg(nknddi)
	real yg(nknddi)
	character*(*) desc(nknddi)

	integer nkn,nlv,nvar

	call read_ets_header1(iu,nkn,nlv,nvar)

	call dimets(iu,nknddi,nlvddi)

	call read_ets_header2(iu,nkn,nlv,ilets,hlv,hets
     +				,nodes,xg,yg,desc)

	end

c***************************************************************

	subroutine read_ets_header1(iu,nkn,nlv,nvar)

c reads first header of ets file

	implicit none

	integer iu
	integer nkn,nlv,nvar

	integer nvers
	integer ierr
	integer date,time
	character*50 title,femver

	nvers = 1

	call ets_init(iu,nvers)

	call ets_read_header(iu,nkn,nlv,nvar,ierr)
	if( ierr .ne. 0 ) goto 99

	call getets(iu,nvers,nkn,nlv,nvar)
	call ets_get_date(iu,date,time)
	call ets_get_title(iu,title)
	call ets_get_femver(iu,femver)

        write(6,*) 'nvers     : ',nvers
        write(6,*) 'nkn       : ',nkn
        write(6,*) 'nlv,nvar  : ',nlv,nvar
        write(6,*) 'title     : ',title
        write(6,*) 'femver    : ',femver
        write(6,*) 'date,time : ',date,time

	return
   99	continue
	write(6,*) 'error in reading header of ETS file'
	stop 'error stop read_ets_header1: reading header'
	end

c***************************************************************

	subroutine read_ets_header2(iu,nkn,nlv,ilets,hlv,hets
     +					,nodes,xg,yg,desc)

c other variables are stored internally

	implicit none

	integer iu
	integer nkn,nlv
	integer ilets(nkn)
	real hlv(nlv)
	real hets(nkn)
	integer nodes(nkn)
	real xg(nkn)
	real yg(nkn)
	character*(*) desc(nkn)

	integer ierr
	integer l,lmax,i,k
	real x,y,h
	character*60 s

	call ets_read_header2(iu,ilets,hlv,hets
     +				,nodes,xg,yg,desc,ierr)
	if( ierr .ne. 0 ) goto 99

        write(6,*) 'Available levels: ',nlv
        write(6,*) (hlv(l),l=1,nlv)

        write(6,*) 'Available nodes: ',nkn
        do i=1,nkn
          k = nodes(i)
	  lmax = ilets(i)
          x = xg(i)
          y = xg(i)
	  h = hets(i)
          write(6,1009) i,k,lmax,x,y,h
        end do
        write(6,*) 'Description of nodes: ',nkn
        do i=1,nkn
          s = desc(i)
          write(6,1008) i,s
	end do
 
	return
   99	continue
	write(6,*) 'error in reading header of ETS file'
	stop 'error stop read_ets_header2: reading header'
 1009   format(i3,i10,i5,3e14.6)
 1008   format(i3,1x,a)
	end

c***************************************************************
c***************************************************************
c***************************************************************

	subroutine open_ets_type(type,status,nunit)

c open ETS file with default simulation name and given extension

	implicit none

	character*(*) type,status
	integer nunit

	integer nb
	character*80 file

        integer ifileo

	call def_make(type,file)
	nb = ifileo(0,file,'unform',status)

	if( nb .le. 0 ) then
	  write(6,*) 'file: ',file
	  stop 'error stop open_ets_type: opening file'
	end if

	nunit = nb

	end

c***************************************************************

	subroutine open_ets_file(name,status,nunit)

	implicit none

	character*(*) name,status
	integer nunit

	integer nb
	character*80 file

        integer ifileo

	call mkname(' ',name,'.ets',file)
	nb = ifileo(0,file,'unform',status)

	if( nb .le. 0 ) then
	  write(6,*) 'file: ',file
	  stop 'error stop open_ets_file: opening file'
	end if

	nunit = nb

	end

c***************************************************************

        subroutine qopen_ets_file(text,status,nunit)

c asks for name and opens ets file

        implicit none

        character*(*) text,status
        integer nunit

        character*80 name

        write(6,*) text
        read(5,'(a)') name
        write(6,*) name

        call open_ets_file(name,status,nunit)

        end

c***************************************************************
c***************************************************************
c***************************************************************

	subroutine ets_get_it_start(file,itstart)

c gets it of first record

	implicit none

	character*(*) file
	integer itstart

	integer nunit,nvers,nvar
	integer it,ivar,ierr
	character*80 title

	nvers = 1
	itstart = -1

	call open_ets_file(file,'old',nunit)
	call ets_init(nunit,nvers)
	call ets_skip_header(nunit,nvar,ierr)
	if( ierr .ne. 0 ) return
	call ets_skip_record(nunit,it,ivar,ierr)
	if( ierr .ne. 0 ) return
	itstart = it

	end

c***************************************************************

	subroutine ets_get_it_end(file,itend)

c gets it of last record

	implicit none

	character*(*) file
	integer itend

	integer nunit,nvers,nvar
	integer it,itlast,ivar,ierr
	character*80 title

	nvers = 1
	itend = -1
	itlast = -1

	call open_ets_file(file,'old',nunit)
	call ets_init(nunit,nvers)
	call ets_skip_header(nunit,nvar,ierr)
	if( ierr .ne. 0 ) return

    1	continue
	call ets_skip_record(nunit,it,ivar,ierr)
	if( ierr .gt. 0 ) return
	if( ierr .lt. 0 ) goto 2
	itlast = it
	goto 1
    2	continue
	itend = itlast

	end

c***************************************************************

	subroutine ets_get_vars(nin,nvar,ivars)

	implicit none

	integer nin
	integer nvar
	integer ivars(nvar)

	integer i,ivar,it,ierr

	do i=1,nvar
	  call ets_next_record(nin,it,ivar,ierr)
	  if( ierr .ne. 0 ) goto 99
	  ivars(i) = ivar
	end do

	do i=1,nvar
	  call ets_back_record(nin)
	end do

	return
   99	continue
	write(6,*) 'not enough variables: ',nvar,i,ierr
	stop 'error stop ets_get_vars: nvar'
	end

c***************************************************************

