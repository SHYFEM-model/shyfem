
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

c revision log :
c
c 18.11.1998    ggu     check dimensions with dimnos
c 06.04.1999    ggu     some cosmetic changes
c 03.12.2001    ggu     some extra output -> place of min/max
c 09.12.2003    ggu     check for NaN introduced
c 07.03.2007    ggu     easier call
c 23.04.2015    ggu     for new version 5
c
c**************************************************************

	program lgrinf

c reads lgr file

	implicit none

	integer nread,nin,i,it,nc
	integer mtype,nvers
	integer nbdy,nn,nout
	integer id,ie,ies,lmax,lb
	real x,y,z,xs,ys
	integer id_special
	character*80 file

	integer iapini
	integer ifem_open_file

c--------------------------------------------------------------

	nread=0
	id_special = 4
	id_special = 1
	id_special = 0

c--------------------------------------------------------------
c open simulation
c--------------------------------------------------------------

	nin = 1
	nc = command_argument_count()
        if( nc == 0 ) stop 'no file given'
        call get_command_argument(1,file)
	open(nin,file=file,status='old',form='unformatted')

	read(nin) mtype,nvers
	write(6,*) 'mtype,nvers: ',mtype,nvers
	if( nvers > 4 ) read(nin) lmax
	if( mtype .ne. 367265 ) stop 'error stop: mtype'

c--------------------------------------------------------------
c loop on data
c--------------------------------------------------------------

	do while(.true.)

	   read(nin,end=100) it,nbdy,nn,nout
	   write(6,*) it,nbdy,nn,nout

	   nread = nread + 1

	   do i=1,nn
	     if( nvers < 5 ) then
	       read(nin) id,x,y,z,ie,xs,ys,ies
	     else
	       read(nin) id,x,y,z,ie,lb,xs,ys,ies
	     end if
	     !if( ie .lt. 0 ) write(6,*) -ie,x,y
	     if( id_special < 0 .or. id == id_special ) then
	       write(6,*) id,x,y,z,lb
	     end if
	   end do

	end do	!do while

c--------------------------------------------------------------
c end of loop on data
c--------------------------------------------------------------

  100	continue

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	end

c***************************************************************

