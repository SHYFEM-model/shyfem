
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
c 22.04.2021	ggu	resolve bound check error (not yet finished)
c
c****************************************************************

        subroutine do_partition(nkn,nel,nen3v,nparts,npart,epart)

c partitions grd file using the METIS library
c requires metis-5.1.0.tar.gz

	implicit none

	integer nkn,nel
	integer nen3v(3,nel)
	integer nparts
	integer npart(nkn)
	integer epart(nel)

	integer		      :: ie,ii,l
	integer, allocatable  :: eptr(:) 		!index for eind
	integer, allocatable  :: eind(:) 		!nodelist in elements
	integer, pointer      :: vwgt(:)=>null() 	!weights of vertices
        integer, pointer      :: vsize(:)=>null()	!size of the nodes
        real(kind=8), pointer :: tpwgts(:)=>null()	!desired weight 

        integer               :: objval		!edge-cut or total comm vol
	integer               :: options(40)	!metis options

c-----------------------------------------------------------------
c initialiaze arrays
c-----------------------------------------------------------------

	npart = 0
	epart = 0

        allocate(eptr(nel+1))
        allocate(eind(nel*3))

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
c end of routine
c-----------------------------------------------------------------

        end

c*******************************************************************

