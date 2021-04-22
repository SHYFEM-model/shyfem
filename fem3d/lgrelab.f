
!--------------------------------------------------------------------------
!
!    Copyright (C) 2008,2010-2011,2014-2015,2019  Georg Umgiesser
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
c 06.06.2008	ggu	new from scratch
c 23.03.2010	ggu	changed v6.1.1
c 22.11.2011	ggu	changed VERS_6_1_37
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_3
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 16.02.2019	ggu	changed VERS_7_5_60
c
c**************************************************************

	program lgrelab

c reads nos file

	use mod_depth
	use basin

	implicit none

	include 'param.h'

c--------------------------------------------------



c--------------------------------------------------

	integer nread,nin,i,it
	integer mtype,nvers
	integer nbdy,nn,nout
	integer ip,ie,ies
	integer nlast,nfirst
	real x,y,xs,ys

	integer ibocche(0:3)
	integer, allocatable :: ielems(:,:)

	integer iapini
	integer ifem_open_file

c--------------------------------------------------------------

	nread=0
	nfirst = 0

c--------------------------------------------------------------
c open simulation
c--------------------------------------------------------------

	if(iapini(3,0,0,0).eq.0) then
		stop 'error stop : iapini'
	end if

	nin = ifem_open_file('.lgr','old')

	read(nin) mtype,nvers
	write(6,*) 'mtype,nvers: ',mtype,nvers
	if( mtype .ne. 367265 ) stop 'error stop: mtype'

	call populate_strings

	allocate(ielems(0:3,nel))

	do i=0,3
	  ibocche(i) = 0
	  do ie=1,nel
	    ielems(i,ie) = 0
	  end do
	end do

	do ie=1,nel
	  hev(ie) = 1.
	end do

c--------------------------------------------------------------
c loop on data
c--------------------------------------------------------------

	do while(.true.)

	   read(nin,end=100) it,nbdy,nn,nout
	   write(6,*) it,nbdy,nn,nout

	   nread = nread + 1

	   do i=1,nn
	     read(nin) ip,x,y,ie,xs,ys,ies
	     if( ie .lt. 0 ) call elab_out(x,y,ies,ielems,ibocche)
	   end do

	   nlast = nn
	   if( nfirst .eq. 0 ) nfirst = nn

	end do	!do while

c--------------------------------------------------------------
c end of loop on data
c--------------------------------------------------------------

  100	continue

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)

	write(6,*) 'still in system: ',nlast,' out of ',nfirst
	write(6,*) 'ibocche: ',ibocche

	call elab_elems(ielems)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	end

c***************************************************************

	subroutine elab_out(x,y,ies,ielems,ibocche)

	implicit none

	real x,y
	integer ies
	integer ielems(0:3,1)
	integer ibocche(0:3)

	logical is_in_rect

	integer ib

	if( is_in_rect(x,y,38500.,32500.,39500.,33500.) ) then		!lido
	  ib = 1
	else if( is_in_rect(x,y,31000.,23000.,31500.,24000.) ) then	!mala
	  ib = 2
	else if( is_in_rect(x,y,29000.,11500.,29500.,12500.) ) then	!chio
	  ib = 3
	else
	  ib = 0
	end if

	ielems(ib,ies) = ielems(ib,ies) + 1
	ibocche(ib) = ibocche(ib) + 1

	end


c***************************************************************

	function is_in_rect(x,y,xmin,ymin,xmax,ymax)

	implicit none

	logical is_in_rect
	real x,y
	real xmin,ymin,xmax,ymax

	if( xmin .le. x .and. x .le. xmax ) then
	  if( ymin .le. y .and. y .le. ymax ) then
	    is_in_rect = .true.
	    return
	  end if
	end if

	is_in_rect = .false.

	end

c***************************************************************

	subroutine elab_elems(ielems)

	use basin

	implicit none

	integer ielems(0:3,1)

	include 'param.h'



	integer nodes(0:3,nkn)
	integer icount(nkn)
	real vals(nkn)

	integer ic,ietot,ievar,ies,i3var
	integer ie,i,k,ii,iu
	real tot,val,valtot

	do i=0,3
	  do k=1,nkn
	    nodes(i,k) = 0
	  end do
	end do

	ietot = 0
	ievar = 0
	i3var = 0

	do ie=1,nel
	  ies = 0
	  do i=1,3
	    ic = ielems(i,ie)
	    if( ic .gt. 0 ) then
	      ies = ies + 1
	      do ii=1,3
	        k = nen3v(ii,ie)
	        nodes(i,k) = nodes(i,k) + 1
	        nodes(0,k) = nodes(0,k) + 1
	      end do
	    end if
	  end do
	  if( ies .gt. 0 ) then
	    ietot = ietot + 1
	    if( ies .gt. 1 ) ievar = ievar + 1
	    if( ies .gt. 2 ) i3var = i3var + 1
	  end if
	end do

	write(6,*) 'elab_elems: ',i3var,ievar,ietot,nel

	iu = 0

	do i=1,3
	  do k=1,nkn
	    tot = nodes(0,k)
	    val = nodes(i,k)
	    if( tot .gt. 0. ) vals(k) = 100*val/tot
	    !write(77,*) i,k,tot,val,vals(k)
	  end do
	  call conwrite2d(iu,'bocche.nos',i,1,543,vals)
	end do

	do k=1,nkn
	  icount(k) = 0
	end do

	do k=1,nkn
	  ies = 0
	  tot = nodes(0,k)
	  valtot = 1.
	  do i=1,3
	    val = nodes(i,k)
	    if( val .gt. 0. ) then
		valtot = valtot * val / tot
		ies = ies + 1
	    end if
	  end do
	  icount(k) = ies
	  if( ies .le. 1 .or. ies .eq. 3 ) valtot = 0.
	  vals(k) = 4.*100.*valtot
	end do

	call conwrite2d(iu,'bocche.nos',4,1,543,vals)

	do k=1,nkn
	  vals(k) = icount(k)
	end do

	call conwrite2d(iu,'bocche.nos',5,1,543,vals)

	end

c***************************************************************

        subroutine conwrite2d(iu,name,it,nvar,ivar,c)

c shell for writing file unconditionally to disk

	use mod_depth
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        integer iu              !unit (0 for first call, set on return)
        character*(*) name      !name of file
	integer it
        integer nvar            !total number of variables
        integer ivar            !id of variable to be written
        real c(1)         !concentration to write

	include 'param.h'


	integer ierr,nvers,nlv
	integer ilhkv(1)
	real hlv(1)
	character*80 title

	integer ifileo

        if( iu .eq. 0 ) then
	  nvers = 3
	  nlv = 1
	  title = 'created file'
	  iu = ifileo(iu,name,'unform','new')
	  call whnos        (iu,nvers
     +                          ,nkn,nel,nlv,nvar
     +                          ,ilhkv,hlv,hev
     +                          ,title
     +                          )
        end if

	call wrnos(iu,it,ivar,1,ilhkv,c,ierr)

        end

c*************************************************************

