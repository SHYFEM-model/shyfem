
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
! 24.05.2020	ggu	debug option added
! 28.05.2020	ggu	some more informational messages
! 15.07.2020	ggu	bug fix for counting elements
! 22.04.2021	ggu	resolve bound check error (not yet finished)
! 12.04.2022	ggu	preapred for online partitioning
!
!****************************************************************

	subroutine info_partition(nparts,area_node)

! write partition information to terminal

	use basin

	implicit none

	integer nparts
	integer area_node(nkn)

	integer ic,k,ia
	integer netot,neint
	integer min,max

        integer, allocatable  :: nc(:)                  !array for check
        integer, allocatable  :: ne(:)                  !array for check
        integer, allocatable  :: ni(:)                  !array for check

	write(6,*) 'writing information on partion to terminal...'

        allocate(nc(0:nparts))
        allocate(ne(0:nparts))
        allocate(ni(0:nparts))

        nc = 0
	netot = 0
	neint = 0
	min = minval(area_node)
	max = maxval(area_node)
	write(6,*) 'min/max: ',min,max
        if( min < 0 .or. max > nparts ) then
          write(6,*) 'nparts: ',nparts
          write(6,*) 'min/max: ',min,max
          stop 'error stop info_partition: internal error (1)'
        end if
	do k=1,nkn
          ic = area_node(k)
          nc(ic) = nc(ic) + 1
	end do
	do ic=min,max
	  call count_elements(nkn,nel,nen3v,ic,area_node,netot,neint)
	  !write(6,*) nel,netot,neint,(100.*neint)/netot
	  !write(6,*) 'ic,netot,neint:',ic,netot,neint
	  ne(ic) = netot
	  ni(ic) = neint
        end do
        write(6,*) 
        write(6,*) 'total number of nodes: ',nkn
        write(6,*) 'total number of elems: ',nel
        write(6,*) 
        write(6,*) 'Information on domains: ',nparts
        write(6,*) 
        write(6,*) '   domain      area     nodes   percent'
     +				// '  elements     ghost'
     +				// '   percent'
        do ic=min,max
	  ia = ic
	  if( min == 0 ) ia = ic + 1
          write(6,'(3i10,f10.2,2i10,f10.2)') 
     +		 ia-1,ia,nc(ic),(100.*nc(ic))/nkn
     +		,ne(ic),ne(ic)-ni(ic),(100.*(ne(ic)-ni(ic)))/ne(ic)
        end do
        write(6,*) 

	end

!****************************************************************

	subroutine write_partition_to_grd(grdfile,bdebug
     +			,nparts,npart,epart)

! write grd files

	use basin
	use grd

	implicit none

	character*(*) grdfile
	logical bdebug
	integer nparts
	integer npart(nkn)
	integer epart(nel)

	character*3 numb
	character*80 basnam,name

	write(6,*) 'writing grd-file...'
	call basin_to_grd
	call grd_set_write(.false.)

	write(numb,'(i3)') nparts
        numb = adjustl(numb)
        basnam = grdfile
        call delete_extension(basnam,'.grd')

        ianv = npart
        name = trim(basnam)//'.'//trim(numb)//'.'//'node.grd'
        call grd_write(name)
	write(6,*) ''
	write(6,*) 'Grid with partition on nodes in file: ',trim(name)

        iaev = epart
        name = trim(basnam)//'.'//trim(numb)//'.'//'elem.grd'
        call grd_write(name)
	write(6,*) 'Grid with partition on elements in file: ',trim(name) 

	if( bdebug ) then
	  call grd_write_debug(basnam,nparts,npart)
	end if

        end

!****************************************************************

	subroutine grd_write_debug(basnam,nparts,npart)

	use basin
        use grd

	implicit none

	character*(*) basnam
	integer nparts
	integer npart(nkn)

	logical bhasnode
	integer i,k,ie,ii
	character*80 name,pre,post,numb
	integer epart(nel)

	write(numb,'(i3)') nparts
        numb = adjustl(numb)
	pre = trim(basnam)//'.'//trim(numb)//'.'
	post = '.node.grd'

	write(6,*) 'writing debug grd files...'

	do i=1,nparts
	  write(numb,'(i3)') i
          numb = adjustl(numb)
          name = trim(pre)//trim(numb)//trim(post)
	  write(6,*) 'writing debug file: ',trim(name)
          ianv = npart
	  where( ianv /= i ) ianv = 0
	  iaev = 0
	  do ie=1,nel
	    bhasnode = .false.
	    do ii=1,3
	      k = nen3v(ii,ie)
	      if( ianv(k) == i ) bhasnode = .true.
	    end do
	    if( bhasnode ) iaev(ie) = i
	  end do
          call grd_write(name)
	end do

	end

c*******************************************************************

