
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017-2020  Georg Umgiesser
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
c 31.08.2017	ggu	started from scratch
c 02.09.2017	ggu	changed VERS_7_5_31
c 30.01.2018	ggu	new routine for combining records
c 22.02.2018	ggu	changed VERS_7_5_42
c 16.02.2019	ggu	changed VERS_7_5_60
c 01.04.2020    ggu     new routines to write regular fem file
c 22.04.2020    ggu     module procedures introduced
c 18.05.2020    ggu     check read/write of files, flag in structure
c 10.07.2020    ggu     compiler warnings resolved (do not init arrays)
c 26.01.2024    ggu     bug fix in femutil_read_record_ff() with optional
c
c**************************************************************
c**************************************************************
c**************************************************************

!==================================================================
	module fem_util
!==================================================================

	implicit none

	integer, parameter, private  :: nvers0 = 3	!highest version of fem

	real, parameter, private  :: flag_fem = -999.

	type :: femfile_type
	  integer :: iunit = 0
	  integer :: nvers = 0
	  integer :: iformat = 0
	  logical :: bread = .true.
	end type femfile_type

	type :: femrec_type
	  logical :: bchanged = .true.
	  real :: flag = flag_fem
	  integer :: np = 0
	  integer :: lmax = 0
	  integer :: nvar = 0
	  integer :: ntype = 0
	  integer :: datetime(2)
	  real :: regpar(7)
	  double precision :: dtime = 0.
	  double precision :: atime = 0.
	  real, allocatable :: hlv(:)
	  integer, allocatable :: ilhkv(:)
	  real, allocatable :: hd(:)
	  real, allocatable :: data(:,:,:)
	  character*80, allocatable :: strings(:)
	end type femrec_type

	type :: fem_type
	  type(femfile_type) :: femfile
	  type(femrec_type)  :: femrec
	end type fem_type

        INTERFACE femutil_get_time
        MODULE PROCEDURE femutil_get_time_rec,femutil_get_time_fem
        END INTERFACE

        INTERFACE femutil_set_time
        MODULE PROCEDURE femutil_set_time_rec,femutil_set_time_fem
        END INTERFACE

        INTERFACE femutil_open_for_write
        MODULE PROCEDURE femutil_open_for_write_fem
     +			,femutil_open_for_write_file
        END INTERFACE

        INTERFACE femutil_open_for_read
        MODULE PROCEDURE femutil_open_for_read_fem
     +			,femutil_open_for_read_file
        END INTERFACE

        INTERFACE femutil_close
        MODULE PROCEDURE femutil_close_fem
     +			,femutil_close_file
        END INTERFACE

        INTERFACE femutil_write_record
        MODULE PROCEDURE femutil_write_record_fem
     +			,femutil_write_record_ff
        END INTERFACE

        INTERFACE femutil_read_record
        MODULE PROCEDURE femutil_read_record_fem
     +			,femutil_read_record_ff
        END INTERFACE

        INTERFACE femutil_peek_time
        MODULE PROCEDURE femutil_peek_time_fem
     +			,femutil_peek_time_rec
        END INTERFACE

        INTERFACE femutil_get_flag
        MODULE PROCEDURE femutil_get_flag_fem
     +			,femutil_get_flag_rec
        END INTERFACE

        INTERFACE femutil_is_regular
        MODULE PROCEDURE femutil_is_regular_fem
     +			,femutil_is_regular_rec
        END INTERFACE

        INTERFACE femutil_copy
        MODULE PROCEDURE femutil_copy_fem
     +			,femutil_copy_rec
        END INTERFACE

!==================================================================
	contains
!==================================================================

	subroutine femutil_init_record(frec)

	type(femrec_type) :: frec

	frec%bchanged = .true.
	frec%flag = flag_fem
	frec%np = 0
	frec%lmax = 0
	frec%nvar = 0
	frec%ntype = 0
	frec%datetime = 0
	frec%regpar = 0.
	frec%dtime = 0.
	frec%atime = 0.
	if( allocated(frec%hlv) ) deallocate(frec%hlv)
	if( allocated(frec%strings) ) deallocate(frec%strings)
	if( allocated(frec%ilhkv) ) deallocate(frec%ilhkv)
	if( allocated(frec%hd) ) deallocate(frec%hd)
	if( allocated(frec%data) ) deallocate(frec%data)

	end subroutine

!******************************************************************

	subroutine femutil_alloc_record(frec,np,lmax,nvar)

	type(femrec_type) :: frec
	integer np,lmax,nvar

	call femutil_init_record(frec)

	frec%np = np
	frec%lmax = lmax
	frec%nvar = nvar

	allocate(frec%hlv(lmax))
	allocate(frec%strings(nvar))
	allocate(frec%ilhkv(np))
	allocate(frec%hd(np))
	allocate(frec%data(lmax,np,nvar))

	end subroutine

!******************************************************************

	subroutine femutil_alloc_record_np(frec,np)

	type(femrec_type) :: frec
	integer np

	integer lmax,nvar

	frec%np = np
	lmax = frec%lmax
	nvar = frec%nvar

	if( allocated(frec%ilhkv) ) deallocate(frec%ilhkv)
	if( allocated(frec%hd) ) deallocate(frec%hd)
	if( allocated(frec%data) ) deallocate(frec%data)
	allocate(frec%ilhkv(np))
	allocate(frec%hd(np))
	allocate(frec%data(lmax,np,nvar))

	end subroutine

!******************************************************************

	function femutil_is_compatible(frec1,frec2)

! checks if the two fem records are compatible
!
! dtime,datetime,atime,data may differ, also string may differ

	logical femutil_is_compatible
	type(femrec_type) :: frec1
	type(femrec_type) :: frec2

	femutil_is_compatible = .false.

	if( frec1%np /= frec2%np ) return
	if( frec1%lmax /= frec2%lmax ) return
	if( frec1%nvar /= frec2%nvar ) return
	if( frec1%ntype /= frec2%ntype ) return

	if( any(frec1%regpar/=frec2%regpar) ) return
	if( any(frec1%hlv/=frec2%hlv) ) return
	if( any(frec1%ilhkv/=frec2%ilhkv) ) return
	if( any(frec1%hd/=frec2%hd) ) return
	!if( any(frec1%strings/=frec2%strings) ) return

	femutil_is_compatible = .true.

	end function

!******************************************************************

	subroutine femutil_get_time_fem(ffem,atime)
	type(fem_type) :: ffem
	double precision atime
	call femutil_get_time_rec(ffem%femrec,atime)
	end subroutine

	subroutine femutil_get_time_rec(frec,atime)
	! returns atime from datetime/dtime
	type(femrec_type) :: frec
	double precision atime
        call dts_convert_to_atime(frec%datetime,frec%dtime,atime)
        frec%atime = atime
	end subroutine

!******************************************************************

	subroutine femutil_set_time_fem(ffem,atime)
	type(fem_type) :: ffem
	double precision atime
	call femutil_set_time_rec(ffem%femrec,atime)
	end subroutine

	subroutine femutil_set_time_rec(frec,atime)
	! sets new datetime/dtime from atime
	type(femrec_type) :: frec
	double precision atime
        call dts_convert_from_atime(frec%datetime,frec%dtime,atime)
        frec%atime = atime
	end subroutine

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine femutil_open_for_write_fem(file,iformat,ffem)
	character*(*) file
	integer iformat
	type(fem_type) :: ffem
	call femutil_open_for_write_file(file,iformat,ffem%femfile)
	call femutil_init_record(ffem%femrec)
	end subroutine

	subroutine femutil_open_for_write_file(file,iformat,ffile)
	character*(*) file
	integer iformat
	type(femfile_type) :: ffile
	integer iunit
	ffile = femfile_type(0,0,0)
	call fem_file_write_open(file,iformat,iunit)
	ffile = femfile_type(iunit,nvers0,iformat)
	ffile%bread = .false.
	end subroutine

!******************************************************************

	subroutine femutil_open_for_read_fem(file,nexp,ffem,ierr)
	character*(*) file
	integer nexp
	type(fem_type) :: ffem
	integer ierr
	call femutil_open_for_read_file(file,nexp,ffem%femfile,ierr)
	call femutil_init_record(ffem%femrec)
	end subroutine

	subroutine femutil_open_for_read_file(file,nexp,ffile,ierr)

	character*(*) file
	integer nexp
	type(femfile_type) :: ffile
	integer ierr

	integer iunit,iformat

	ierr = 1
	ffile = femfile_type(0,0,0)

	call fem_file_read_open(file,nexp,iformat,iunit)

	if( iunit <= 0 ) return

	ierr = 0
	ffile = femfile_type(iunit,nvers0,iformat)
	
	end subroutine

!******************************************************************

	subroutine femutil_close_fem(ffem)
	type(fem_type) :: ffem
	call femutil_close(ffem%femfile)
	end subroutine

	subroutine femutil_close_file(ffile)
	type(femfile_type) :: ffile
	close(ffile%iunit)
	ffile = femfile_type(0,0,0)
	end subroutine

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine femutil_write_record_fem(ffem)
	type(fem_type) :: ffem
	call femutil_write_record(ffem%femfile,ffem%femrec)
	end subroutine

	subroutine femutil_write_record_ff(ffile,frec)

	type(femfile_type) :: ffile
	type(femrec_type) :: frec

	integer iunit,iformat,nvers
	integer lmax,nvar,np,ntype
	integer iv
	character*80 file

	iunit = ffile%iunit
	nvers = ffile%nvers
	iformat = ffile%iformat

	if( ffile%bread ) then
	  inquire(unit=iunit, name=file)
	  write(6,*) 'file name: ',trim(file)
	  write(6,*) 'file has been opened for read... cannot write'
	  stop 'error stop femutil_write_record: inconsistency'
	end if

	lmax = frec%lmax
	nvar = frec%nvar
	np = frec%np
	ntype = frec%ntype

	call fem_file_write_header(iformat,iunit
     +				,frec%dtime
     +                          ,nvers
     +				,np
     +				,lmax
     +                          ,nvar
     +				,ntype
     +                          ,lmax
     +				,frec%hlv
     +				,frec%datetime
     +				,frec%regpar)

	do iv=1,nvar
	  call fem_file_write_data(iformat,iunit
     +                          ,nvers
     +				,np
     +				,lmax
     +                          ,frec%strings(iv)
     +                          ,frec%ilhkv
     +				,frec%hd
     +                          ,lmax
     +				,frec%data(:,:,iv))
	end do

	end subroutine

!******************************************************************

	subroutine femutil_read_record_fem(ffem,ierr,bskip)
	type(fem_type) :: ffem
	integer :: ierr
	logical, optional :: bskip
	call femutil_read_record_ff(ffem%femfile,ffem%femrec,ierr,bskip)
	end subroutine

	subroutine femutil_read_record_ff(ffile,frec,ierr,bskip)

	type(femfile_type) :: ffile
	type(femrec_type) :: frec
	integer :: ierr
	logical, optional :: bskip

	logical brealloc,bbskip
	integer iunit,iformat,nvers
	integer lmax,nvar,np,ntype
	integer iv
	integer llmax(frec%nvar)
	character*80 file

	bbskip = .false.
	if( present(bskip) ) bbskip = bskip

	iunit = ffile%iunit
	iformat = ffile%iformat

	if( .not. ffile%bread ) then
	  inquire(unit=iunit, name=file)
	  write(6,*) 'file name: ',trim(file)
	  write(6,*) 'file has been opened for write... cannot read'
	  stop 'error stop femutil_read_record: inconsistency'
	end if

	call fem_file_read_params(iformat,iunit
     +				,frec%dtime
     +                          ,nvers
     +				,np
     +				,lmax
     +                          ,nvar
     +				,ntype
     +				,frec%datetime
     +				,ierr)
	if( ierr /= 0 ) return

	ffile%nvers = nvers
	frec%ntype = ntype
	brealloc = .false.

	if( lmax /= frec%lmax ) then
	  if( allocated(frec%hlv) ) deallocate(frec%hlv)
	  allocate(frec%hlv(lmax))
	  frec%lmax = lmax
	  brealloc = .true.
	end if
	if( nvar /= frec%nvar ) then
	  if( allocated(frec%strings) ) deallocate(frec%strings)
	  allocate(frec%strings(nvar))
	  frec%nvar = nvar
	  brealloc = .true.
	end if
	if( np /= frec%np ) then
	  if( allocated(frec%ilhkv) ) deallocate(frec%ilhkv)
	  if( allocated(frec%hd) ) deallocate(frec%hd)
	  allocate(frec%ilhkv(np))
	  allocate(frec%hd(np))
	  frec%np = np
	  brealloc = .true.
	end if
	if( brealloc ) then
	  if( allocated(frec%data) ) deallocate(frec%data)
	  allocate(frec%data(lmax,np,nvar))
	  frec%bchanged = .true.
	else
	  frec%bchanged = .false.
	end if

        call fem_file_read_2header(iformat,iunit,ntype,lmax
     +                  ,frec%hlv,frec%regpar,ierr)
	if( ierr /= 0 ) return

	if( femutil_is_regular(frec) ) then
	  frec%flag = frec%regpar(7)
	else
	  frec%flag = flag_fem
	end if

	if( ntype /= frec%ntype ) then
	  frec%ntype = ntype
	  frec%bchanged = .true.
	end if

	do iv=1,nvar
	  if( bbskip ) then
            call fem_file_skip_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,frec%strings(iv),ierr)
	    llmax(iv) = lmax
	  else
	    call fem_file_read_data(iformat,iunit
     +                          ,nvers
     +				,np
     +				,llmax(iv)
     +                          ,frec%strings(iv)
     +                          ,frec%ilhkv
     +				,frec%hd
     +                          ,lmax
     +				,frec%data(:,:,iv)
     +                          ,ierr)
	  end if
	  if( ierr /= 0 ) return
	end do

	if( any( llmax /= lmax ) ) then
	  write(6,*) '*** inconsistency of lmax'
	  write(6,*) 'lmax: ',lmax
	  write(6,*) 'lmax(iv): ',llmax
	  ierr = 5
	end if

	end subroutine

!******************************************************************

	subroutine femutil_peek_time_fem(ffem,atime,ierr)
	type(fem_type) :: ffem
	double precision atime
	integer ierr
	call femutil_peek_time_rec(ffem%femfile,atime,ierr)
	end subroutine

	subroutine femutil_peek_time_rec(ffile,atime,ierr)
	type(femfile_type) :: ffile
	double precision atime
	integer ierr
        double precision dtime
        integer nvers,np,lmax,nvar,ntype
        integer datetime(2)
	atime = -1
	call fem_file_peek_params(ffile%iformat,ffile%iunit,dtime
     +                      ,nvers,np,lmax,nvar,ntype,datetime,ierr)
	if( ierr /= 0 ) return
	call dts_convert_to_atime(datetime,dtime,atime)
	end subroutine

!******************************************************************

	function femutil_get_flag_fem(ffem)
	real femutil_get_flag_fem
	type(fem_type) :: ffem
	femutil_get_flag_fem = femutil_get_flag_rec(ffem%femrec)
	end function

	function femutil_get_flag_rec(frec)
	real femutil_get_flag_rec
	type(femrec_type) :: frec
	if( femutil_is_regular(frec) ) then
	  femutil_get_flag_rec = frec%regpar(7)
	else
	  femutil_get_flag_rec = flag_fem
	end if
	end function

!******************************************************************

	function femutil_is_regular_fem(ffem)
	logical femutil_is_regular_fem
	type(fem_type) :: ffem
	femutil_is_regular_fem = femutil_is_regular_rec(ffem%femrec)
	end function

	function femutil_is_regular_rec(frec)
	logical femutil_is_regular_rec
	type(femrec_type) :: frec
	integer itype(2)
	call fem_file_make_type(frec%ntype,2,itype)
	femutil_is_regular_rec = ( itype(2) .gt. 0 )
	end function

!******************************************************************

	subroutine femutil_copy_fem(ffem_from,ffem_to)
	type(fem_type) :: ffem_from,ffem_to
	ffem_to%femrec = ffem_from%femrec
	end subroutine

	subroutine femutil_copy_rec(frec_from,frec_to)
	type(femrec_type) :: frec_from,frec_to
	frec_to = frec_from
	end subroutine

!******************************************************************

	subroutine femfile_info(ffinfo)

	type(femfile_type) :: ffinfo

	integer iformat,iunit
	character*20 fline
	character*80 file

        iformat = ffinfo%iformat
        iunit = ffinfo%iunit
	inquire(unit=iunit, name=file)
        call fem_file_get_format_description(iformat,fline)

        write(6,*) 'info on femfile'
        write(6,*) 'file:   ',trim(file)
        write(6,*) 'iunit:  ',ffinfo%iunit
        write(6,*) 'nvers:  ',ffinfo%nvers
        write(6,*) 'format: ',iformat,"  (",trim(fline),")"
        write(6,*) 'bread:  ',ffinfo%bread
        write(6,*) 'info on femfile finished'

	end

!******************************************************************

	subroutine femutil_info(ffem)

	type(fem_type) :: ffem

	integer iformat,iunit
	integer itype(2)
	character*20 fline
	character*80 file

        iformat = ffem%femfile%iformat
        iunit = ffem%femfile%iunit
	inquire(unit=iunit, name=file)
        call fem_file_get_format_description(iformat,fline)
        call fem_file_make_type(ffem%femrec%ntype,2,itype)

        write(6,*) 'file:   ',trim(file)
        write(6,*) 'nvers:  ',ffem%femfile%nvers
        write(6,*) 'format: ',iformat,"  (",trim(fline),")"
        write(6,*) 'np:     ',ffem%femrec%np
        write(6,*) 'lmax:   ',ffem%femrec%lmax
        write(6,*) 'nvar:   ',ffem%femrec%nvar
        write(6,*) 'ntype:  ',ffem%femrec%ntype

        if( ffem%femrec%lmax > 1 ) then
          write(6,*) 'levels: ',ffem%femrec%lmax
          write(6,'(5f12.4)') ffem%femrec%hlv
        end if
        if( itype(1) .gt. 0 ) then
          write(6,*) 'date and time: ',ffem%femrec%datetime
        end if
        if( itype(2) .gt. 0 ) then
          write(6,*) 'regpar: '
          call printreg(ffem%femrec%regpar)
        end if
 
	end subroutine

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine femutil_combine_data_recs(nfile,frec,fout)

	implicit none

	integer nfile
	type(femrec_type) :: frec(nfile)
	type(femrec_type) :: fout

	integer i,nvar,lmax,np,nv,ip

	nvar = 0
	do i=1,nfile
	  nvar = nvar + frec(i)%nvar
	end do

	fout = frec(1)
	fout%nvar = nvar
	lmax = fout%lmax
	np = fout%np

	deallocate(fout%strings)
	deallocate(fout%data)
	allocate(fout%strings(nvar))
	allocate(fout%data(lmax,np,nvar))

	ip = 0
	do i=1,nfile
	  nv = frec(i)%nvar
	  fout%strings(ip+1:ip+nv) = frec(i)%strings(1:nv)
	  fout%data(:,:,ip+1:ip+nv) = frec(i)%data(:,:,1:nv)
	  ip = ip + nv
	end do

	end subroutine

!******************************************************************

	subroutine femutil_add_data_recs(nfile,frec,fout)

	implicit none

	integer nfile
	type(femrec_type) :: frec(nfile)
	type(femrec_type) :: fout

	integer i,nvar,lmax,np
	real, parameter :: flag = -999.

	fout = frec(1)
	nvar = fout%nvar
	lmax = fout%lmax
	np   = fout%np

	deallocate(fout%data)
	allocate(fout%data(lmax,np,nvar))
	fout%data = 0.

	do i=1,nfile
	  where( frec(i)%data /= flag )
	    fout%data = fout%data + frec(i)%data
	  end where
	  !write(6,*) i,frec(i)%data(1,1,1),frec(i)%data(lmax,np,nvar)
	end do
	!write(6,*) i,fout%data(1,1,1),fout%data(lmax,np,nvar)

	where( frec(1)%data == flag )
	  fout%data = flag
	end where

	end subroutine

!==================================================================
	end module fem_util
!==================================================================

!------------------------------------------------------------------
! utility routines
!------------------------------------------------------------------

	subroutine write_regular_2d_1var_record(iunit
     +			,string
     +			,regpar
     +			,np,data)

	implicit none

	integer iunit
	character*(*) string
	real regpar(7)
	integer np
	real data(np)

	integer nvers,lmax,nvar,ntype,iformat
	integer datetime(2)
	integer ilhkv(1)
	real hlv(1)
	real hd(1)
	double precision dtime

	dtime = 0
	datetime = 0
	nvers = 0
	lmax = 1
	nvar = 1
	ntype = 10
	iformat = 1
	hlv = 10000.
	ilhkv = 1
	hd = 10000.

	call fem_file_write_header(iformat,iunit
     +				,dtime
     +                          ,nvers
     +				,np
     +				,lmax
     +                          ,nvar
     +				,ntype
     +                          ,lmax
     +				,hlv
     +				,datetime
     +				,regpar)

	!do iv=1,nvar
	  call fem_file_write_data(iformat,iunit
     +                          ,nvers
     +				,np
     +				,lmax
     +                          ,string
     +                          ,ilhkv
     +				,hd
     +                          ,lmax
     +				,data)
	!end do

	end subroutine

!******************************************************************

	subroutine write_regular_2d_nodes(file
     +			,string
     +			,dreg
     +			,data)

	implicit none

	character*(*) file,string
	real dreg
	real data(*)	!must be nkn

	integer nvers,lmax,nvar,ntype,iformat
	integer, parameter :: iunit = 1
	integer np,nx,ny
	real xmin,ymin,dx,dy,flag
	real regpar(7)
	real, allocatable :: rdata(:)

        call make_reg_box(dreg,regpar)
        call getreg(regpar,nx,ny,xmin,ymin,dx,dy,flag)
        call setgeo(xmin,ymin,dx,dy,flag)

        np = nx*ny
        allocate( rdata(np) )

        call av2am(data,rdata,nx,ny)

        open(iunit,file=file,status='unknown',form='formatted')
        call write_regular_2d_1var_record(iunit,string,regpar,np,rdata)
        close(iunit)

	end

!******************************************************************
!
!	program fem_util_test_main
!	end program
!
!******************************************************************

