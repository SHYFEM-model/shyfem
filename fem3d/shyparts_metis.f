
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

c revision log :
c
c 19.05.2020	ccf	started from scratch
c 24.05.2020	ggu	debug option added
c 28.05.2020	ggu	some more informational messages
c 15.07.2020	ggu	bug fix for counting elements
c
c****************************************************************

        program shyparts

c partitions grd file using the METIS library
c requires metis-5.1.0.tar.gz

	use mod_geom
	use evgeom
	use basin
	use clo
	use basutil
	use shympi
        use grd

	implicit none

	integer		      :: k,ie,ii,l,ic
	integer		      :: min,max
        integer               :: nparts			!number of parts to partition the mesh
	integer, allocatable  :: eptr(:) 		!index for eind
	integer, allocatable  :: eind(:) 		!list of nodes in elements
	integer, pointer      :: vwgt(: )=>null() 	!weights of the vertices
        integer, pointer      :: vsize=>null()		!size of the nodes
        real(kind=8), pointer :: tpwgts(: )=>null()	!desired weight for each partition

        integer               :: objval			!edge-cut or the total comm vol
	integer, allocatable  :: epart(:)		!partition vector for the elements
	integer, allocatable  :: npart(:)		!partition vector for the nodes
	integer               :: options(40)		!metis options
        integer, allocatable  :: nc(:)			!array for check
        integer, allocatable  :: ne(:)			!array for check
        integer, allocatable  :: ni(:)			!array for check

	logical bdebug
	integer kerr,ierr1,ierr2
	integer netot,neint
        character*3 numb
	character*80 grdfile,basnam,name

c-----------------------------------------------------------------
c read in basin
c-----------------------------------------------------------------

        call shyfem_copyright('shyparts - partitioning a SHYFEM grid')
	write(6,*) 'Partition with METIS'
	write(6,*) ''

        call shyparts_init(grdfile,nparts,bdebug)

        if( grdfile == ' ' ) call clo_usage
	call grd_set_write(.false.)
        call read_command_line_file(grdfile)

        call shympi_init(.false.)

c-----------------------------------------------------------------
c initialize modules
c-----------------------------------------------------------------

        call ev_init(nel)
        call set_ev

        call mod_geom_init(nkn,nel,ngr)
        call set_geom

c-----------------------------------------------------------------
c initialiaze arrays
c-----------------------------------------------------------------

        allocate(eptr(nel+1))
        allocate(eind(nel*3))
        allocate(epart(nel))
        allocate(npart(nkn))
        epart = 0
        npart = 0

c-----------------------------------------------------------------
c set up METIS eptr and eind arrays structures
c-----------------------------------------------------------------

        eptr(1)=1
        do ie = 1,nel
          do ii = 1,3
            l = eptr(ie) + ii -1
            eind(l) = nen3v(ii,ie)
          end do
          eptr(ie+1) = eptr(ie) + 3
        enddo 

c-----------------------------------------------------------------
c Set METIS options
c For the full list and order of options see metis-5.1.0/include/metis.h
c and the documentation in manual/manual.pdf
c-----------------------------------------------------------------

        call METIS_SetDefaultOptions(options)
        options(1) = 1		!PTYPE (0=rb,1=kway)
        options(2) = 1		!OBJTYPE (0=cut,1=vol)
        options(12) = 1		!CONTIG (0=defoult,1=force contiguous)
        options(18) = 1		!NUMBERING (0 C-style, 1 Fortran-style)

c-----------------------------------------------------------------
c Call METIS for patitioning on nodes
c-----------------------------------------------------------------

	write(6,*) 'partitioning with METIS...'
        call METIS_PartMeshNodal(nel, nkn, eptr, eind, vwgt, vsize, 
     +       nparts, tpwgts, options, objval, epart, npart)

c-----------------------------------------------------------------
c check partition
c-----------------------------------------------------------------

	write(6,*) 'checking connectivity and connections...'

	call link_set_stop(.false.)	!do not stop after error
	call link_set_write(.false.)	!do not write error

        iarnv = npart
        iarv = epart
        call check_connectivity(ierr1)
        call check_connections(ierr2)
	npart = iarnv
	epart = iarv

	call link_set_stop(.true.)

!-----------------------------------------------------------------
! write partition information to terminal
!-----------------------------------------------------------------

	write(6,*) 'writing information on partion to terminal...'
        allocate(nc(0:nparts))
        allocate(ne(0:nparts))
        allocate(ni(0:nparts))
        nc = 0
	netot = 0
	neint = 0
	min = minval(iarnv)
	max = maxval(iarnv)
        if( min < 1 .or. max > nparts ) then
          write(6,*) 'ic,nparts: ',ic,nparts
          stop 'error stop bas_partition: internal error (1)'
        end if
	do k=1,nkn
          ic = iarnv(k)
          nc(ic) = nc(ic) + 1
	end do
	do ic=1,nparts
	  call count_elements(nkn,nel,nen3v,ic,iarnv,netot,neint)
	  !write(6,*) nel,netot,neint,(100.*neint)/netot
	  ne(ic) = netot
	  ni(ic) = neint
        end do
        write(6,*) 
        write(6,*) 'total number of nodes: ',nkn
        write(6,*) 'total number of elems: ',nel
        write(6,*) 
        write(6,*) 'Information on domains: ',nparts
        write(6,*) 
        write(6,*) '   domain     nodes   percent  elements     ghost'
     +				//'   percent'
        do ic=1,nparts
          write(6,'(2i10,f10.2,2i10,f10.2)') 
     +		 ic,nc(ic),(100.*nc(ic))/nkn
     +		,ne(ic),ne(ic)-ni(ic),(100.*(ne(ic)-ni(ic)))/ne(ic)
        end do
        write(6,*) 

c-----------------------------------------------------------------
c write grd files
c-----------------------------------------------------------------

	write(6,*) 'writing grd-file...'
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

	if( ierr1 == 0 .and. ierr2 == 0 ) then
	  write(6,*) 'domain successfully partitioned'
	  call exit(0)
	else
	  write(6,*) 'errors in domains: ',ierr1,ierr2
	  call exit(9)
	end if

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

        end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine read_command_line_file(file)

	use basin
	use basutil

	implicit none

	character*(*) file
	logical is_grd_file

	if( basin_is_basin(file) ) then
	  write(6,*) 'reading BAS file ',trim(file)
	  call basin_read(file)
	  breadbas = .true.
	else if( is_grd_file(file) ) then
	  write(6,*) 'reading GRD file ',trim(file)
	  call grd_read(file)
	  call grd_to_basin
	  call estimate_ngr(ngr)
	  breadbas = .false.
	else
	  write(6,*) 'Cannot read this file: ',trim(file)
	  stop 'error stop read_given_file: format not recognized'
	end if

	end

c*******************************************************************

	subroutine shyparts_init(grdfile,nparts,bdebug)

	use clo

	implicit none

	character*(*) grdfile
        integer nparts
	logical bdebug

	call clo_init('shyparts','grd-file','3.0')

        call clo_add_info('partitioning of grd file with METIS')

	call clo_add_sep('options for partitioning')
        call clo_add_option('nparts',-1,'number of partitions')
        call clo_add_option('debug',.false.,'write debug grd files')

	call clo_parse_options

	call clo_check_files(1)
	call clo_get_file(1,grdfile)

        call clo_get_option('nparts',nparts)
        call clo_get_option('debug',bdebug)

        if (nparts < 2 ) then
          write(6,*) 'nparts: ',nparts
	  stop 'error stop shyparts_init: nparts < 2'
        end if

	end

c*******************************************************************

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

	subroutine node_test
	end

c*******************************************************************

