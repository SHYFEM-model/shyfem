
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017-2019  Georg Umgiesser
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
! 19.05.2020	ccf	started from scratch
! 12.04.2022	ggu	adapted
!
!****************************************************************

        subroutine do_partition(nkn,nel,nen3v,nparts,npart,epart)

! shyparts dummy routine

        implicit none

        integer nkn,nel
        integer nen3v(3,nel)
        integer nparts
        integer npart(nkn)
        integer epart(nel)

        write(6,*)' For automatic partitioning of a grid install' 
        write(6,*)' one of the following libraries and set the'
        write(6,*)' parameters PARTS and PARTSDIR in the'
        write(6,*)' Rules.make configuration file' 
        write(6,*)'   - METIS'
        write(6,*)' Then recompile: "make fem"'

	!stop 'error stop do_partition: no metis available'
	npart = 0
	epart = 0
	write(6,*) 'running do_custom with np = ',nparts
        call do_custom(nkn,nel,nen3v,nparts,npart)

	end

!*******************************************************************

	subroutine check_partition(npart,epart,ierr1,ierr2)

	use basin

	implicit none

        integer npart(nkn)
        integer epart(nel)
	integer ierr1,ierr2

	ierr1 = 0
	ierr2 = 0

	end

!*******************************************************************

        subroutine do_custom(nkn,nel,nen3v,nparts,npart)

! shyparts custom routine

	use mod_color

        implicit none

        integer nkn,nel
        integer nen3v(3,nel)
        integer nparts
        integer npart(nkn)

	integer np,ngr,ngrmin,k,ncol,nmax
	integer kroots(nparts+1)
	integer, allocatable :: ngrade(:),egrade(:)
	integer, allocatable :: epart(:)
	integer, allocatable :: ncolor(:)
	integer, allocatable :: matrix(:,:)

	allocate(ngrade(nkn),egrade(nkn))
	allocate(epart(nel))
	allocate(ncolor(nkn))

!------------------------------------------------------------
! initialize connectivity routines
!------------------------------------------------------------

	call connect_init(nkn,nel,nen3v)
	call connect_get_grades(nkn,ngrade,egrade,ngr)

!------------------------------------------------------------
! find one node with lowest grade
!------------------------------------------------------------

	ngrmin = minval(ngrade)
	do k=1,nkn
	  if( ngrade(k) == ngrmin ) exit
	end do
	if( k > nkn ) stop 'error stop: k>nkn'

	kroots(1) = k
	npart(k) = 1
	epart = 0

!------------------------------------------------------------
! loop over partitioning
!------------------------------------------------------------

	do np=1,nparts
	  write(6,*) 'calling ffill ',np
	  call ffill(nkn,nel,nen3v,ngrade,np,kroots,npart)
	  call make_epart(nkn,nel,nen3v,npart,epart)
	  call check_connected_domains(nkn,nel,npart,epart)
	  call make_node_matrix(nel,nen3v,nkn,npart,nmax,matrix)
	  call color_graph(nkn,npart,nmax,matrix,ncol,ncolor)
	  !call color_nodes(nel,nen3v,nkn,npart,ncol,ncolor)
	  call release_color_matrix(matrix)
	  write(6,*) 'ncol = ',ncol
	  write(6,*) 'npart min/max: ',minval(npart),maxval(npart)
	  write(6,*) 'epart min/max: ',minval(epart),maxval(epart)
	  call write_partition_to_grd('domain_custom',.false.
     +		,np,npart,epart)
	  call write_partition_to_grd('domain_color',.false.
     +		,np,ncolor,epart)
	end do

!------------------------------------------------------------
! end of routine
!------------------------------------------------------------

	end

!*******************************************************************

	subroutine ffill(nkn,nel,nen3v,ngrade,np,kroots,npart)

! shyparts custom routine

        implicit none

        integer nkn,nel
        integer nen3v(3,nel)
        integer ngrade(nkn)
        integer np
	integer kroots(np+1)
        integer npart(nkn)

	integer ie,ii,k,is,ncol,ic
	integer i,nchange,ia,dmax,id,iloop
	integer kn,ng
	integer kk(3)

	integer color1(nkn)
	integer color2(nkn)
	integer idist(nkn)

	color1 = 0
	idist = 0

	do i=1,np
	  k = kroots(i)
	  color1(k) = i
	  idist(k) = 1
	end do
	color2 = color1

	iloop = 0

	do
	  iloop = iloop + 1
	  id = iloop + 1
	  nchange = 0
	  color1 = color2
	  do ie=1,nel
	    ncol = 0
	    is = 0
	    do ii=1,3
	      k = nen3v(ii,ie)
	      kk(ii) = k
	      if( color1(k) /= 0 ) then
		ncol = ncol + 1
		is = is + ii
	      end if
	    end do
	    if( ncol > 0 .and. ncol < 3 ) then
	      nchange = nchange + 1
	      if( ncol == 1 ) then
	        ic = color1(kk(is))
	      else
		is = 6 - is
		ia = 1+mod(is,3)	!use any color of the two available
	        ic = color1(kk(ia))
	      end if
	      if( ic <= 0 ) stop 'error stop (1)'
	      do ii=1,3
		k = kk(ii)
		if( color1(k) == 0 ) then
		  color2(k) = ic
		  idist(k) = id
		end if
	      end do
	    end if
	  end do
	  if( nchange == 0 ) exit
	end do

	npart = color1
	dmax = maxval(idist)

	kn = 0
	ng = nkn
	do k=1,nkn
	  if( idist(k) == dmax ) then
	    if( ngrade(k) < ng ) then
	      ng = ngrade(k)
	      kn = k
	    end if
	  end if
	end do
	if( kn == 0 ) stop 'error stop: kn==0'

	kroots(np+1) = kn

	end

!*******************************************************************

	subroutine make_epart(nkn,nel,nen3v,npart,epart)

! make epart from npart - only elements fully in domain are colored, other -1

	implicit none

	integer nkn,nel
	integer nen3v(3,nel)
	integer npart(nkn)		!nodal partition [1-np]
	integer epart(nel)		!elemental partition [1-np] (return)

	integer ie,ii,k,ic,ics(3)

	epart = -1

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    ics(ii) = npart(k)
	  end do
	  ic = ics(1)
	  if( all( ics == ic ) ) epart(ie) = ic
	end do

	end

!*******************************************************************

	subroutine check_connected_domains(nkn,nel,npart,epart)

	use mod_flood

	implicit none

	integer nkn,nel
	integer npart(nkn)
	integer epart(nel)
	
	integer nmax,i,iae

	nmax = maxval(npart)
	write(6,*) 'stats: ',nmax
	call domain_stats(nmax,nel,epart)

	do i=0,nmax
	  call floodfill_elem(nel,epart,i,iae)
	  write(6,*) 'domain by elems: ',i,iae
	  if( iae > 1 ) goto 99
	  call floodfill_elem(nkn,npart,i,iae)
	  write(6,*) 'domain by nodes: ',i,iae
	  if( iae > 1 ) goto 99
	end do

	return
   99	continue
	write(6,*) 'domain is not connected: ',i,iae
	stop 'error stop check_connected_domains: not connected'
	end

!*******************************************************************

	subroutine domain_stats(nmax,n,part)

	implicit none

	integer nmax,n
	integer part(n)

	integer i,ic,nmin
	integer, allocatable :: count(:)

	nmin = -1
	allocate(count(nmin:nmax))

	count = 0

	do i=1,n
	  ic = part(i)
	  if( ic < nmin .or. ic > nmax ) goto 99
	  count(ic) = count(ic) + 1
	end do

	write(6,*) '----------------------'
	do i=nmin,nmax
	  write(6,*) i,count(i)
	end do
	write(6,*) '----------------------'

	return
   99	continue
	write(6,*) i,ic
	stop 'error stop domain_stats: value not possible'
	end

!*******************************************************************

