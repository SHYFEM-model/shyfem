
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

c utilities for reading/writing fem files
c
c contents :
c
c revision log :
c
c 31.08.2017    ggu     started from scratch
c 30.01.2018    ggu     new routine for combining records
c
c**************************************************************
c**************************************************************
c**************************************************************

!==================================================================
	module fem_util
!==================================================================

	implicit none

	integer, parameter, private  :: nvers0 = 3	!highest version of fem

	type :: femfile_type
	  integer :: iunit
	  integer :: nvers
	  integer :: iformat
	end type femfile_type

	type :: fem_type
	  double precision :: dtime,atime
	  integer :: datetime(2)
	  integer :: np,lmax,nvar,ntype
	  real :: regpar(7)
	  real, allocatable :: hlv(:)
	  character*80, allocatable :: strings(:)
	  integer, allocatable :: ilhkv(:)
	  real, allocatable :: hd(:)
	  real, allocatable :: data(:,:,:)
	end type fem_type

!==================================================================
	contains
!==================================================================

	subroutine femutil_init_record(finfo)

	type(fem_type) :: finfo

	finfo%dtime = 0.
	finfo%atime = 0.
	finfo%datetime = 0
	finfo%regpar = 0.
	finfo%np = 0
	finfo%lmax = 0
	finfo%nvar = 0
	finfo%ntype = 0
	if( allocated(finfo%hlv) ) deallocate(finfo%hlv)
	if( allocated(finfo%strings) ) deallocate(finfo%strings)
	if( allocated(finfo%ilhkv) ) deallocate(finfo%ilhkv)
	if( allocated(finfo%hd) ) deallocate(finfo%hd)
	if( allocated(finfo%data) ) deallocate(finfo%data)

	end subroutine

!******************************************************************

	subroutine femutil_alloc_record(finfo,np,lmax,nvar)

	type(fem_type) :: finfo
	integer np,lmax,nvar

	call femutil_init_record(finfo)

	finfo%np = np
	finfo%lmax = lmax
	finfo%nvar = nvar

	allocate(finfo%hlv(lmax))
	allocate(finfo%strings(nvar))
	allocate(finfo%ilhkv(np))
	allocate(finfo%hd(np))
	allocate(finfo%data(lmax,np,nvar))

	end subroutine

!******************************************************************

	function femutil_is_compatible(finfo1,finfo2)

! checks if the two fem records are compatible
!
! dtime,datetime,atime,data may differ, also string may differ

	logical femutil_is_compatible
	type(fem_type) :: finfo1
	type(fem_type) :: finfo2

	femutil_is_compatible = .false.

	if( finfo1%np /= finfo2%np ) return
	if( finfo1%lmax /= finfo2%lmax ) return
	if( finfo1%nvar /= finfo2%nvar ) return
	if( finfo1%ntype /= finfo2%ntype ) return

	if( any(finfo1%regpar/=finfo2%regpar) ) return
	if( any(finfo1%hlv/=finfo2%hlv) ) return
	if( any(finfo1%ilhkv/=finfo2%ilhkv) ) return
	if( any(finfo1%hd/=finfo2%hd) ) return
	!if( any(finfo1%strings/=finfo2%strings) ) return

	femutil_is_compatible = .true.

	end function

!******************************************************************

	subroutine femutil_get_time(finfo,atime)

! returns atime from datetime/dtime

	type(fem_type) :: finfo
	double precision atime

        call dts_convert_to_atime(finfo%datetime,finfo%dtime,atime)
        finfo%atime = atime

	end subroutine

!******************************************************************

	subroutine femutil_set_time(finfo,atime)

! sets new datetime/dtime from atime

	type(fem_type) :: finfo
	double precision atime

        call dts_convert_from_atime(finfo%datetime,finfo%dtime,atime)
        finfo%atime = atime

	end subroutine

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine femutil_open_for_write(file,iformat,ffinfo)

	character*(*) file
	integer iformat
	type(femfile_type) :: ffinfo

	integer iunit

	ffinfo = femfile_type(0,0,0)

	call fem_file_write_open(file,iformat,iunit)

	ffinfo%iunit = iunit
	ffinfo%nvers = nvers0
	ffinfo%iformat = iformat
	
	end subroutine

!******************************************************************

	subroutine femutil_open_for_read(file,nexp,ffinfo,ierr)

	character*(*) file
	integer nexp
	type(femfile_type) :: ffinfo
	integer ierr

	integer iunit,iformat

	ierr = 1
	ffinfo = femfile_type(0,0,0)

	call fem_file_read_open(file,nexp,iformat,iunit)

	if( iunit <= 0 ) return

	ierr = 0
	ffinfo%iunit = iunit
	ffinfo%nvers = nvers0
	ffinfo%iformat = iformat
	
	end subroutine

!******************************************************************

	subroutine femutil_close(ffinfo)

	type(femfile_type) :: ffinfo

	close(ffinfo%iunit)
	ffinfo = femfile_type(0,0,0)

	end subroutine

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine femutil_write_record(ffinfo,finfo)

	type(femfile_type) :: ffinfo
	type(fem_type) :: finfo

	integer iunit,iformat,nvers
	integer lmax,nvar,np,ntype
	integer iv

	iunit = ffinfo%iunit
	nvers = ffinfo%nvers
	iformat = ffinfo%iformat

	lmax = finfo%lmax
	nvar = finfo%nvar
	np = finfo%np
	ntype = finfo%ntype

	call fem_file_write_header(iformat,iunit
     +				,finfo%dtime
     +                          ,nvers
     +				,np
     +				,lmax
     +                          ,nvar
     +				,ntype
     +                          ,lmax
     +				,finfo%hlv
     +				,finfo%datetime
     +				,finfo%regpar)

	do iv=1,nvar
	  call fem_file_write_data(iformat,iunit
     +                          ,nvers
     +				,np
     +				,lmax
     +                          ,finfo%strings(iv)
     +                          ,finfo%ilhkv
     +				,finfo%hd
     +                          ,lmax
     +				,finfo%data(:,:,iv))
	end do

	end subroutine

!******************************************************************

	subroutine femutil_read_record(ffinfo,finfo,ierr)

	type(femfile_type) :: ffinfo
	type(fem_type) :: finfo
	integer ierr

	logical brealloc
	integer iunit,iformat,nvers
	integer lmax,nvar,np,ntype
	integer iv

	iunit = ffinfo%iunit
	iformat = ffinfo%iformat

	call fem_file_read_params(iformat,iunit
     +				,finfo%dtime
     +                          ,nvers
     +				,np
     +				,lmax
     +                          ,nvar
     +				,ntype
     +				,finfo%datetime
     +				,ierr)
	if( ierr /= 0 ) return

	ffinfo%nvers = nvers
	finfo%ntype = ntype
	brealloc = .false.

	if( lmax /= finfo%lmax ) then
	  if( allocated(finfo%hlv) ) deallocate(finfo%hlv)
	  allocate(finfo%hlv(lmax))
	  finfo%lmax = lmax
	  brealloc = .true.
	end if
	if( nvar /= finfo%nvar ) then
	  if( allocated(finfo%strings) ) deallocate(finfo%strings)
	  allocate(finfo%strings(nvar))
	  finfo%nvar = nvar
	  brealloc = .true.
	end if
	if( nvar /= finfo%np ) then
	  if( allocated(finfo%ilhkv) ) deallocate(finfo%ilhkv)
	  if( allocated(finfo%hd) ) deallocate(finfo%hd)
	  allocate(finfo%ilhkv(np))
	  allocate(finfo%hd(np))
	  finfo%np = np
	  brealloc = .true.
	end if
	if( brealloc ) then
	  if( allocated(finfo%data) ) deallocate(finfo%data)
	  allocate(finfo%data(lmax,np,nvar))
	end if

        call fem_file_read_2header(iformat,iunit,ntype,lmax
     +                  ,finfo%hlv,finfo%regpar,ierr)
	if( ierr /= 0 ) return

	do iv=1,nvar
	  call fem_file_read_data(iformat,iunit
     +                          ,nvers
     +				,np
     +				,lmax
     +                          ,finfo%strings(iv)
     +                          ,finfo%ilhkv
     +				,finfo%hd
     +                          ,lmax
     +				,finfo%data(:,:,iv)
     +                          ,ierr)
	  if( ierr /= 0 ) return
	end do

	end subroutine

!******************************************************************

	subroutine femutil_combine_data(nfile,finfo,fout)

	implicit none

	integer nfile
	type(fem_type) :: finfo(nfile)
	type(fem_type) :: fout

	integer i,nvar,lmax,np,nv,ip

	nvar = 0
	do i=1,nfile
	  nvar = nvar + finfo(i)%nvar
	end do

	fout = finfo(1)
	fout%nvar = nvar
	lmax = fout%lmax
	np = fout%np

	deallocate(fout%strings)
	deallocate(fout%data)
	allocate(fout%strings(nvar))
	allocate(fout%data(lmax,np,nvar))

	ip = 0
	do i=1,nfile
	  nv = finfo(i)%nvar
	  fout%strings(ip+1:ip+nv) = finfo(i)%strings(1:nv)
	  fout%data(:,:,ip+1:ip+nv) = finfo(i)%data(:,:,1:nv)
	  ip = ip + nv
	end do

	end subroutine

!==================================================================
	end module fem_util
!==================================================================
!
!	program fem_util_test_main
!	end program
!
!******************************************************************

