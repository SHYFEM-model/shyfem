
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

! module for boundary information
!
! revision log :
!
! 01.12.2015    ggu     ready to accept no open boundary (bug fix)
! 01.04.2016    ggu     parameter and file names transfered to here
!
!******************************************************************

!==================================================================
        module mod_bnd
!==================================================================

        implicit none

        integer, private, save :: nbc_bnd = 0

	integer, parameter :: nbvdim = 27	!total number of values
	integer, parameter :: nbfdim = 16	!total number of files

        character*6, save :: bnd_par_names(nbvdim) =
     +			(/
     +                   'iqual ','ibtyp ','kranf ','krend ','zval  '
     +                  ,'ampli ','period','phase ','zref  ','ktilt '
     +                  ,'intpol','levmax','igrad0','zfact '
     +                  ,'conz  ','temp  ','salt  ','levmin','kref  '
     +                  ,'ztilt ','tramp ','levflx','nad   ','lgrpps'
     +                  ,'tnudge','kmanf ','kmend '
     +                  /)

	character*6, save :: bnd_file_names(nbfdim) =
     +			(/
     +                   'boundn','conzn ','saltn ','tempn '
     +                  ,'bio2dn','tox3dn','mercn '
     +                  ,'sed2dn','mud2dn','lam2dn','dmf2dn'
     +                  ,'bfm1bc','bfm2bc','bfm3bc','bfmbcn'
     +                  ,'vel3dn'
     +                  /)

	integer, save :: nbc = 0
	integer, save :: nrb = 0

        real, allocatable, save :: bnd(:,:)
        character*80, allocatable, save :: bnd_file(:,:)

!==================================================================
	contains
!==================================================================

        subroutine mod_bnd_init(nb)

        integer nb

	!write(6,*) 'mod_bnd_init: ',nb,nb,nbc_bnd

        if( nb == nbc_bnd ) return

        if( nbc_bnd > 0 ) then
          deallocate(bnd)
          deallocate(bnd_file)
        end if

        nbc_bnd = nb

        if( nb == 0 ) return

        allocate(bnd(nbvdim,nb))
        allocate(bnd_file(nbfdim,nb))

        end subroutine mod_bnd_init

!***************************************************************

	subroutine mod_bnd_adjust(nb)

        integer nb

        integer ndim
        real, allocatable :: bnd_aux(:,:)
        character*80, allocatable :: bnd_file_aux(:,:)

	!write(6,*) 'mod_bnd_adjust: ',nb,nbc_bnd

        ndim = nbc_bnd

        if( ndim == 0 ) then
          ndim = max(10,2*nb)
          allocate(bnd(nbvdim,ndim))
          allocate(bnd_file(nbfdim,ndim))
        else if( nb > ndim ) then
          ndim = ndim*2
          allocate(bnd_aux(nbvdim,ndim))
          allocate(bnd_file_aux(nbfdim,ndim))
          bnd_aux(:,1:ndim/2) = bnd(:,1:ndim/2)
          bnd_file_aux(:,1:ndim/2) = bnd_file(:,1:ndim/2)
          call move_alloc(bnd_aux,bnd)
          call move_alloc(bnd_file_aux,bnd_file)
        end if

	nbc_bnd = ndim

	end subroutine mod_bnd_adjust

!***************************************************************

        subroutine mod_bnd_reinit(nb)

        integer nb
	integer nbb

	real, allocatable :: bnd_aux(:,:)
        character*80, allocatable :: bnd_file_aux(:,:)

	!write(6,*) 'mod_bnd_reinit: ',nb,nbc_bnd

	if( nb > nbc_bnd ) then
	  write(6,*) 'nb,nbc_bnd: ',nb,nbc_bnd
	  stop 'error stop mod_bnd_reinit: nb > nbc_bnd'
	end if

	nbb = nb
	if( nbc_bnd == 0 ) then		!no boundary, also nb == 0
	  nbb = 1
	  nbc_bnd = nbb
          allocate(bnd(nbvdim,nbb))
          allocate(bnd_file(nbfdim,nbb))
	  bnd = 0
	end if

	allocate(bnd_aux(nbvdim,nbb))
	allocate(bnd_file_aux(nbfdim,nbb))
	bnd_aux(:,1:nbb) = bnd(:,1:nbb)
	bnd_file_aux(:,1:nbb) = bnd_file(:,1:nbb)

        call mod_bnd_init(nbb)

	bnd = 0
	bnd_file = ' '
	bnd(:,1:nbb) = bnd_aux(:,1:nbb)
	bnd_file(:,1:nbb) = bnd_file_aux(:,1:nbb)
	deallocate(bnd_aux)
	deallocate(bnd_file_aux)

	nbc_bnd = nbb

        end subroutine mod_bnd_reinit

!==================================================================
        end module mod_bnd
!==================================================================

