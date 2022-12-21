
!--------------------------------------------------------------------------
!
!    Copyright (C) 2012-2020  Georg Umgiesser
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

! FEM file management routines
!
! contents :
!
! revision log :
!
! 02.10.2012	ggu	created from scratch
! 05.11.2012	ggu	changed VERS_6_1_60
! 17.12.2012	ggu	changed VERS_6_1_61a
! 16.05.2013	ggu	better documentation
! 13.06.2013	ggu	changed VERS_6_1_65
! 12.09.2013	ggu	changed VERS_6_1_67
! 24.04.2014	ggu	use nvar>0 as indication of good read
! 05.05.2014	ggu	changed VERS_6_1_74
! 30.05.2014	ggu	restructured
! 16.06.2014	ggu	time is now double precision
! 27.06.2014	ggu	changed VERS_6_1_78
! 07.07.2014	ggu	first version consolidated
! 20.10.2014	ggu	second version (date record is just after first record)
! 29.10.2014	ggu	new routine fem_file_is_fem_file()
! 05.11.2014	ggu	changed VERS_7_0_5
! 26.11.2014	ggu	changed VERS_7_0_7
! 05.12.2014	ggu	changed VERS_7_0_8
! 12.12.2014	ggu	changed VERS_7_0_9
! 09.01.2015	ggu	new routine fem_file_get_format_description()
! 14.01.2015	ggu	new routine fem_file_string2time()
! 26.02.2015	ggu	changed VERS_7_1_5
! 21.05.2015	ggu	changed VERS_7_1_11
! 10.07.2015	ggu	changed VERS_7_1_50
! 05.11.2015	ggu	changed VERS_7_3_12
! 18.12.2015	ggu	changed VERS_7_3_17
! 21.01.2016	ggu	read and write string without leading blanks
! 22.03.2016	ggu	changed VERS_7_5_6
! 15.04.2016	ggu	changed VERS_7_5_8
! 13.05.2016	ggu	nvers = 3 -> add data size to records
! 14.05.2016	ggu	new module to collect global parameters
! 25.05.2016	ggu	changed VERS_7_5_10
! 10.06.2016	ggu	changed VERS_7_5_13
! 14.06.2016	ggu	changed VERS_7_5_14
! 27.06.2016	ggu	changed VERS_7_5_16
! 05.10.2016	ggu	routine to clean data from NaNs
! 12.01.2017	ggu	changed VERS_7_5_21
! 09.05.2017	ggu	changed VERS_7_5_26
! 02.09.2017	ggu	changed VERS_7_5_31
! 04.11.2017	ggu	changed VERS_7_5_34
! 14.11.2017	ggu	changed VERS_7_5_36
! 05.12.2017	ggu	changed VERS_7_5_39
! 03.04.2018	ggu	changed VERS_7_5_44
! 16.02.2019	ggu	changed VERS_7_5_60
! 29.03.2020	ggu	better ntype handling
! 11.11.2020	ggu	bug fix for nvers < 3 (define lmax in data read)
! 23.06.2021	ggu	more documentation
! 27.01.2022	ggu	new routine for skipping record, find last record
!
! notes :
!
! versions
!
! nvers == 1		no regular grid allowed
! nvers == 2		complete specification
! nvers == 3		write np,lmax for each record -> can mix 2d/3d records
!
! format for file (nvers == 3)
!
!	time record 1
!	time record 2
!	time record ...
!
! format for time record
!
!	header record
!	data record for variable 1
!	data record for variable 2
!	data record for variable ...
!	data record for variable nvar
!
! format for header record
!
!	dtime,nvers,idfem,np,lmax,nvar,ntype
!	date,time				for ntype == 1
!	(hlv(l),l=1,lmax)			only if( lmax > 1 )
!	regpar					for ntype == 10
!	other lines depending on ntype
!
! format for data record
!
!	if( lmax == 1 )
!		string
!		np,lmax				only for nvers > 2
!		(data(1,k),k=1,np)
!	if( lmax > 1 )
!		string
!		np,lmax				only for nvers > 2
!		do k=1,np
!		  lm,hd(k),(data(l,k),l=1,lm)
!		end do
!
! legend
!
! dtime		time stamp (double precision, seconds)
! nvers		version of file format
! idfem		id to identify fem file (must be 957839)
! np		number of horizontal points given
! lmax		maximum number of layers given (1 for 2D)
! nvar		number of variables in time record
! ntype		type of data, defines extra data to follow
! date		reference date (integer, YYYYMMDD)
! time		reference time (integer, hhmmss)
! hlv		layer depths (the bottom of each layer)
! string 	string with description of data (character*80)
! ilhkv(k)	total number of levels of node k (1 for 2D)
! hd(k)		total depth in node k (real, -999 if unknown)
! data(l,k)	data for variable at level l and node k (real)
! lm		total number of vertical data for point k
! k,l		index for horizontal/vertical dimension
! regpar	regular grid info: nx,ny,x0,y0,dx,dy,flag
! nx,ny		size of regular grid
! x0,y0		origin of regular grid (real)
! dx,dy		space increment of regular grid (real)
! flag		flag for invalid data of regular grid (real)
!
! file type (ntype)
!
! 0		no other lines in header
! 1		give date/time of reference on extra line
! 10		regular grid, info on extra line (regpar)
! 20		rotated regular grid, information on extra line (not yet ready)
!
! combinations are possible, example
!
! 11		date/time and regular grid
!
! routines to write and read fem files
!
! fem_file_write_header()
! fem_file_write_data()
!
! fem_file_read_params()
! fem_file_read_2header()
! fem_file_read_data()
!
! for writing one time record the calling sequence is
!
!	call fem_file_write_header()	!write header
!	call fem_file_write_data()	!write first data record
!	call fem_file_write_data()	!write second data record
!	...
!	call fem_file_write_data()	!write nvar data record
!
! for reading one time record the calling sequence is
!
!	call fem_file_read_params()	!read parameters of header
!	... allocate arrays
!	call fem_file_read_2header()	!read arrays of header
!	call fem_file_read_data()	!read first data record
!	call fem_file_read_data()	!read second data record
!	...
!	call fem_file_read_data()	!read nvar data record

!==================================================================
        module fem_file
!==================================================================

	implicit none

	integer, save :: nvmax = 3	!newest version
	integer, save :: nvmin = 1	!oldest supported version
	integer, save :: idfem = 957839	!fem file id - do not change

	real, save :: femflag = -999.	!default flag for no data

!==================================================================
        end module fem_file
!==================================================================

c************************************************************
c************************************************************
c************************************************************

	subroutine fem_file_write_header(iformat,iunit,dtime
     +				,nvers,np,lmax
     +				,nvar,ntype
     +				,nlvddi,hlv,datetime,regpar)

c writes header of fem file

        implicit none

	integer iformat		!formatted or unformatted
	integer iunit		!file unit
	double precision dtime	!time stamp
	integer nvers		!version of file format
	integer np		!size of data (horizontal, nodes or elements)
	integer lmax		!maximum vertical values (1 for 2d)
	integer nvar		!number of variables to write
	integer ntype		!type of information contained
	integer nlvddi		!vertical dimension of data
	real hlv(nlvddi)	!depth at bottom of layer
	integer datetime(2)	!date and time parameters
	real regpar(7)		!parameters for regular field

	call fem_file_write_params(iformat,iunit,dtime
     +				,nvers,np,lmax
     +				,nvar,ntype,datetime)

	call fem_file_write_2header(iformat,iunit,ntype,lmax
     +			,hlv,regpar)

	end

c************************************************************

	subroutine fem_file_write_params(iformat,iunit,dtime
     +				,nvers,np,lmax
     +				,nvar,ntype,datetime)

c writes first header of fem file

	use fem_file

        implicit none

	integer iformat		!formatted or unformatted
	integer iunit		!file unit
	double precision dtime	!time stamp
	integer nvers		!version of file format
	integer np		!size of data (horizontal, nodes or elements)
	integer lmax		!maximum vertical values (1 for 2d)
	integer nvar		!number of variables to write
	integer ntype		!type of information contained
	integer datetime(2)	!date and time information

	integer l,nv
	integer itype(2)
	integer(kind=8) :: itlong

	nv = nvers
	if( nv .eq. 0 ) nv = nvmax	!default
	if( nv .ne. nvmax ) goto 99
	if( lmax < 1 ) goto 98

	if( iformat == 1 ) then
	  itlong = dtime
	  write(iunit,1000) itlong,nv,idfem,np,lmax,nvar,ntype
	else
	  write(iunit) dtime,nv,idfem,np,lmax,nvar,ntype
	end if

	call fem_file_make_type(ntype,2,itype)

	if( itype(1) .eq. 1 ) then
	  if( iformat == 1 ) then
	    write(iunit,*) datetime
	  else
	    write(iunit) datetime
	  end if
	end if

	return
 1000	format(i20,i4,i8,i12,i6,i4,i6)
   98	continue
	write(6,*) 'lmax = ',lmax
	stop 'error stop fem_file_write_header: lmax < 1'
   99	continue
	write(6,*) 'nvers = ',nvers
	stop 'error stop fem_file_write_header: nvers'
	end

c************************************************************

	subroutine fem_file_write_2header(iformat,iunit,ntype,lmax
     +			,hlv,regpar)

	implicit none

	integer iformat		!formatted or unformatted
	integer iunit		!file unit
	integer ntype		!type of information contained
	integer lmax		!maximum vertical values (1 for 2d)
	real hlv(lmax)		!vertical structure
	real regpar(7)		!regular array params

	integer l,i
	integer itype(2)

	call fem_file_make_type(ntype,2,itype)

	if( lmax .gt. 1 ) then
	  if( iformat == 1 ) then
	    write(iunit,1000) (hlv(l),l=1,lmax)
	  else
	    write(iunit) (hlv(l),l=1,lmax)
	  end if
	end if

	if( itype(2) .gt. 0 ) then
	  if( iformat == 1 ) then
	    write(iunit,*) regpar
	  else
	    write(iunit) regpar
	  end if
	end if

	return
 1000	format(5g14.6)
	end

c************************************************************

	subroutine fem_file_write_data(iformat,iunit
     +				,nvers,np,lmax
     +				,string
     +				,ilhkv,hd
     +				,nlvddi,data)

c writes data of the file

	use fem_file

        implicit none

	integer iformat		!formatted or unformatted
	integer iunit		!file unit
	integer nvers		!version of file format
	integer np		!size of data (horizontal, nodes or elements)
	integer lmax		!maximum vertical values (1 for 2d)
	character*(*) string	!string explanation
	integer ilhkv(np)	!number of layers in node k
	real hd(np)		!total depth in node k
	integer nlvddi		!vertical dimension of data
	real data(nlvddi,np)	!data

	logical b2d
	integer k,lm,l,nv
	character*60 text
	character*80 textu	!we need 80 chars for unformatted write

	nv = nvers
	if( nv .eq. 0 ) nv = nvmax	!default

	text = adjustl(string)
	textu = text
	b2d = lmax .le. 1

	if( iformat == 1 ) then
	  write(iunit,'(a)') text
	  if( nv >= 3 ) write(iunit,*) np,lmax
	  if( b2d ) then
	    write(iunit,1000) (data(1,k),k=1,np)
	  else
	    do k=1,np
	      lm = min(lmax,ilhkv(k))
	      write(iunit,1010) lm,hd(k),(data(l,k),l=1,lm)
	    end do
	  end if
	else
	  write(iunit) textu
	  if( nv >= 3 ) write(iunit) np,lmax
	  if( b2d ) then
	    write(iunit) (data(1,k),k=1,np)
	  else
	    do k=1,np
	      lm = ilhkv(k)
	      write(iunit) lm,hd(k),(data(l,k),l=1,lm)
	    end do
	  end if
	end if

	return
 1000	format(5g14.6)
 1010	format(i8,(5g14.6))
	end

c************************************************************
c************************************************************
c************************************************************

	function fem_file_is_fem_file(file)

c checks if file is fem file

	implicit none

	logical fem_file_is_fem_file
	character*(*) file	!file name

	integer iformat		!is formatted? -1 for no fem file (return)

	call fem_file_test_fem_file(file,iformat)

	!write(6,*) 'fem_file_is_fem_file: ',trim(file),iformat

	fem_file_is_fem_file = iformat >= 0

	end

c************************************************************

	subroutine fem_file_test_fem_file(file,iformat)

c tries to open fem file for read
c
c returns -1 in iformat if no fem file, else iformat indicates format

	implicit none

	character*(*) file	!file name
	integer iformat		!is formatted? -1 for no fem file (return)

	integer nvar,np,ntype
	logical filex

	iformat = -1

	if( .not. filex(file) ) then
	  write(6,*) 'file does not exist: ',file
	  return
	end if

	call fem_file_test_formatted(file,np,nvar,ntype,iformat)

	end

c************************************************************

	subroutine fem_file_write_open(file,iformat,iunit)

c opens fem file for write

	implicit none

	character*(*) file	!file name
	integer iformat		!formatted write?
	integer iunit		!unit of opened file (return) (0 for error)

	character*80 form
	integer ifileo

        if( iformat == 0 ) then
          form='unformatted'
        else
          form='formatted'
        end if

        iunit = ifileo(60,file,form,'unknown')

	end

c************************************************************

	subroutine fem_file_read_open(file,nexp,iformat,iunit)

c opens fem file for read

	implicit none

	character*(*) file	!file name
	integer nexp		!expected size of data (0 if unknown)
	integer iformat		!is formatted? (return)
	integer iunit		!unit of opened file (return) (0 for error)

	integer nvar,np,ntype
	integer itype(2)
	logical filex,breg

	iunit = 0
	iformat = 0

	if( .not. filex(file) ) then
	  write(6,*) 'file does not exist: ',file
	  return
	end if

	call fem_file_test_formatted(file,np,nvar,ntype,iformat)

	call fem_file_make_type(ntype,2,itype)
	breg = itype(2) > 0

	if( nvar .gt. 0 ) then
	  if( .not. breg .and. nexp .gt. 0 .and. np .ne. nexp ) then
	    write(6,*) 'fem_file_read_open: data not of expected size'
	    write(6,*) 'nvar,nexp,np: ',nvar,nexp,np
	    call fem_file_write_info(file,iformat)
	  else
	    call find_unit(iunit)
	    if( iformat == 1 ) then
	      open(iunit,file=file,form='formatted',status='old')
	    else
	      open(iunit,file=file,form='unformatted',status='old')
	    end if
	  end if
	else
	  write(6,*) 'fem_file_read_open: error opening file '
	  write(6,*) 'nvar = ',nvar
	  call fem_file_write_info(file,1)
	  call fem_file_write_info(file,0)
	end if

	end

c************************************************************

	subroutine fem_file_write_info(file,iformat)

c writes information on file from header

	use fem_file

	implicit none

	character*(*) file	!file name
	integer iformat		!is formatted?

	integer iunit
	double precision dtime
	integer nvers,np,lmax,nvar,ntype
	integer id
	integer ios

	dtime = 0
	np = 0
	lmax = 0
	nvar = 0
	nvers = 0
	ntype = 0

	iunit = 90
	call find_unit(iunit)

	write(6,*) 'debug info for file: ',file

	if( iformat == 1 ) then
	  open(iunit,file=file,form='formatted',status='old')
	  read(iunit,*,iostat=ios) dtime,nvers,id,np,lmax,nvar,ntype
	  write(6,*) 'formatted read: '
	else
	  open(iunit,file=file,form='unformatted',status='old')
	  read(iunit,iostat=ios) dtime,nvers,id,np,lmax,nvar,ntype
	  write(6,*) 'unformatted read: '
	end if

	if( ios .gt. 0 ) then
	  write(6,*) 'fem_file_write_info: error reading file'
	else if( ios .lt. 0 ) then
	  write(6,*) 'fem_file_write_info: EOF found'
	else if( id .ne. idfem ) then
	  write(6,*) 'file is not a FEM file: ',id,idfem
	else
	  write(6,'(g14.2,6i10)') dtime,nvers,id,np,lmax,nvar,ntype
	end if

	close(iunit)

	end

c************************************************************

	subroutine fem_file_test_formatted(file,np,nvar,ntype,iformat)

c checks if file is readable and formatted or unformatted

	implicit none

	character*(*) file	!file name
	integer np		!size of data (return)
	integer nvar		!successful read => nvar>0 (return)
	integer ntype		!type of data
	integer iformat		!is formatted? (return)

	integer iunit
	integer nvers,np0,lmax
	integer datetime(2)
	double precision dtime
	integer ierr
	logical bdebug,bopen,bexist

c------------------------------------------------------
c initialize parameters
c------------------------------------------------------

	bdebug = .true.
	bdebug = .false.

	nvers = 0
	dtime = 0
	np0 = 0
	lmax = 0
	nvar = 0
	ntype = 0

c------------------------------------------------------
c find unit to open file
c------------------------------------------------------

	iunit = 90
	call find_unit(iunit)
	inquire(file=file,exist=bexist)
	if( .not. bexist ) then
	  write(6,*) 'file does not exist: ',trim(file)
	  goto 9
	end if
	!write(6,*) 'try to open file on unit ',iunit
	!write(6,*) 'file name: ',trim(file)

c------------------------------------------------------
c first try unformatted
c------------------------------------------------------

	ierr = 77
	open(iunit,file=file,form='unformatted',status='old',err=2)

	iformat = 0
	call fem_file_read_params(iformat,iunit,dtime
     +				,nvers,np,lmax,nvar,ntype,datetime,ierr)

	close(iunit)
	!write(6,*) 'ierrrr ',ierr

	if( ierr .ne. 0 ) then
	  if( bdebug ) write(6,*) 'unformatted read error ',ierr
	else	!ok, probably unformatted
	  return
	end if

    2	continue
	!inquire(iunit,opened=bopen)
	!write(6,*) 'unformatted open error...',bopen,ierr

c------------------------------------------------------
c now try formatted
c------------------------------------------------------

	ierr = 77
	open(iunit,file=file,form='formatted',status='old',err=8)

	iformat = 1
	call fem_file_read_params(iformat,iunit,dtime
     +				,nvers,np,lmax,nvar,ntype,datetime,ierr)

	close(iunit)
	!write(6,*) 'formatted ierr = ',ierr

	if( ierr .ne. 0 ) then
	  if( bdebug ) write(6,*) 'formatted read error ',ierr
	else	!ok, probably formatted
	  return
	end if

    8	continue
	!inquire(iunit,opened=bopen)
	!write(6,*) 'formatted open error...',bopen,ierr

c------------------------------------------------------
c no successful opening
c------------------------------------------------------

    9	continue

	np = 0
	nvar = 0
	ntype = 0
	iformat = -1
	if( ierr == 77 ) iformat = -77

c------------------------------------------------------
c end of routine
c------------------------------------------------------

	end

c************************************************************

	subroutine fem_file_get_data_description(file
     +			,strings,ierr)

c returns data description for first record

	implicit none

	character*(*) file		!file name
	character*80 strings(*)		!return - must have dimension nvar
	integer ierr

	integer iformat
	integer np0,iunit,i
	integer nvers,np,lmax,nvar,ntype
	integer datetime(2)
	real regpar(7)
	double precision dtime
	character*80 string

	np0 = 0
	ierr = 1

	call fem_file_read_open(file,np0,iformat,iunit)
	if( iunit .le. 0 ) return

	call fem_file_read_params(iformat,iunit,dtime
     +				,nvers,np,lmax
     +				,nvar,ntype,datetime,ierr)
	if( ierr .ne. 0 ) return

	call fem_file_skip_2header(iformat,iunit
     +				,ntype,lmax,regpar,ierr)
	if( ierr .ne. 0 ) return

	do i=1,nvar
	  call fem_file_skip_data(iformat,iunit
     +				,nvers,np,lmax
     +				,string,ierr)
	  if( ierr .ne. 0 ) return
	  strings(i) = adjustl(string)
	end do

	close(iunit)
	ierr = 0

	end

c************************************************************
c************************************************************
c************************************************************

	subroutine fem_file_read_params(iformat,iunit,dtime
     +				,nvers,np,lmax
     +				,nvar,ntype,datetime,ierr)

c reads and checks params of next header

        implicit none

	integer iformat		!formatted or unformatted
	integer iunit		!file unit
	double precision dtime	!time stamp
	integer nvers		!version of file format
	integer np		!size of data (horizontal, nodes or elements)
	integer lmax		!vertical values
	integer nvar		!number of variables to write
	integer ntype		!type of information contained
	integer datetime(2)	!date and time information
	integer ierr		!return error code

	integer id,i
	integer itype(2)
	character*80 string

	ierr = 0
	if( iunit < 1 ) goto 99

	if( iformat == 1 ) then
	  read(iunit,*,end=1,err=2) dtime,nvers,id,np,lmax,nvar,ntype
	else
	  read(iunit,end=1,err=2) dtime,nvers,id,np,lmax,nvar,ntype
	end if

	call fem_file_check_params(nvers,id,np,lmax,nvar,ntype,ierr)
	if( ierr .ne. 0 ) return

	call fem_file_make_type(ntype,2,itype)

	if( nvers == 1 .and. itype(1) > 0 ) goto 98

	if( itype(1) > 0 ) then
	  if( iformat == 1 ) then
	    read(iunit,*,err=3) (datetime(i),i=1,2)
	  else
	    read(iunit,err=3) (datetime(i),i=1,2)
	  end if
	  !if( itype(1) == 1 ) then
	  !  read(string,'(2i10)') datetime
	  !else if( itype(1) == 2 ) then
	  !else
	  !  write(6,*) 'itype(1) = ',itype(1)
	  !  stop 'error stop fem_file_read_params: error itype'
	  !end if
	else
	  datetime = 0
	end if

	return
    1	continue
	backspace(iunit)
	ierr = -1
	return
    2	continue
	!write(6,*) 'error reading fem record header on unit ',iunit
	ierr = 1
	return
    3	continue
	write(6,*) 'error reading fem date-time record on unit ',iunit
	ierr = 3
	return
   98	continue
	write(6,*) 'impossible combination: nvers,ntype: ',nvers,ntype
	stop 'error stop fem_file_read_params: nvers and ntype'
   99	continue
	write(6,*) 'impossible unit number: ',iunit
	stop 'error stop fem_file_read_params: iunit'
	end

c************************************************************

	subroutine fem_file_peek_params(iformat,iunit,dtime
     +				,nvers,np,lmax,nvar,ntype,datetime,ierr)

c reads and checks params of next header (non advancing read)

        implicit none

	integer iformat		!formatted or unformatted
	integer iunit		!file unit
	double precision dtime	!time stamp
	integer nvers		!version of file format
	integer np		!size of data (horizontal, nodes or elements)
	integer lmax		!vertical values
	integer nvar		!number of variables to write
	integer ntype		!type of information contained
	integer datetime(2)	!date and time information
	integer ierr		!return error code

	call fem_file_read_params(iformat,iunit,dtime
     +				,nvers,np,lmax,nvar,ntype,datetime,ierr)

	if( ierr .ne. 0 ) return

	call fem_file_back_params(iunit,nvers,ntype)

	end

c************************************************************

	subroutine fem_file_back_params(iunit,nvers,ntype)

! backspaces first header

	implicit none

	integer iunit,nvers,ntype

	integer itype(2)

	backspace(iunit)
	call fem_file_make_type(ntype,2,itype)
	if( itype(1) > 0 ) backspace(iunit)
	
	end

c************************************************************

	subroutine fem_file_check_params(nvers,id,np,lmax,nvar,ntype,ierr)

c reads and checks params of next header

	use fem_file

        implicit none

	integer nvers		!version of file format
	integer id		!id of fem file
	integer np		!size of data (horizontal, nodes or elements)
	integer lmax		!vertical values
	integer nvar		!number of variables to write
	integer ntype		!type of information contained
	integer ierr		!return error code

	ierr = 11
	if( id .ne. idfem ) goto 9
	ierr = 13
	if( nvers .lt. nvmin .or. nvers .gt. nvmax ) goto 9
	ierr = 15
	if( np .le. 0 ) goto 9
	ierr = 17
	if( lmax .le. 0 .or. lmax .gt. 1000 ) goto 9
	ierr = 19
	if( nvar .le. 0 .or. nvar .gt. 100 ) goto 9
	ierr = 21
	if( ntype .lt. 0 .or. ntype .gt. 21 ) goto 9

	ierr = 0
	return

    9	continue

	return
	end

c************************************************************

	subroutine fem_file_read_2header(iformat,iunit,ntype,lmax
     +			,hlv,regpar,ierr)

c reads hlv of header

	use fem_file

        implicit none

	integer iformat		!formatted or unformatted
	integer iunit		!file unit
	integer ntype		!type of second header
	integer lmax		!vertical dimension
	real hlv(lmax)		!vertical structure
	real regpar(7)		!regular array params
	integer ierr		!return error code

	integer l,i
	integer itype(2)

	call fem_file_make_type(ntype,2,itype)

	ierr = 3
	if( lmax .gt. 1 ) then
	  if( iformat == 1 ) then
	    read(iunit,*,err=1) (hlv(l),l=1,lmax)
	  else
	    read(iunit,err=1) (hlv(l),l=1,lmax)
	  end if
	else
	  hlv(1) = 10000.
	end if

	ierr = 7
	if( itype(2) .gt. 0 ) then
	  if( iformat == 1 ) then
	    read(iunit,*,err=1) (regpar(i),i=1,7)
	  else
	    read(iunit,err=1) (regpar(i),i=1,7)
	  end if
	  femflag = regpar(7)
	else
	  regpar = 0.
	end if

	ierr = 0
	return
    1	continue
	end

c************************************************************

	subroutine fem_file_skip_2header(iformat,iunit
     +				,ntype,lmax,regpar,ierr)

c skips additional headers in fem file

        implicit none

	integer iformat		!formatted or unformatted
	integer iunit		!file unit
	integer ntype		!type of second header
	integer lmax		!total number of elements to read
	real regpar(7)		!regular array params
	integer ierr		!return error code

	integer l,i
	integer itype(2)
	real aux

	call fem_file_make_type(ntype,2,itype)

	ierr = 3
	if( lmax .gt. 1 ) then			!read hlv
	  if( iformat  == 1 ) then
	    read(iunit,*,err=1) (aux,l=1,lmax)
	  else
	    read(iunit,err=1) (aux,l=1,lmax)
	  end if
	end if

	ierr = 7
	if( itype(2) .gt. 0 ) then		!read regpar
	  if( iformat == 1 ) then
	    read(iunit,*,err=1) (regpar(i),i=1,7)
	  else
	    read(iunit,err=1) (regpar(i),i=1,7)
	  end if
	else
	  regpar = 0.
	end if

	ierr = 0
	return
    1	continue
	end

c************************************************************
c************************************************************
c************************************************************

	subroutine fem_file_read_data(iformat,iunit
     +				,nvers,np,lmax
     +				,string
     +				,ilhkv,hd
     +				,nlvddi,data
     +				,ierr)

c reads data of the file

	use fem_file

        implicit none

	integer iformat		!formatted or unformatted
	integer iunit		!file unit
	integer nvers		!version of file format
	integer np		!size of data (horizontal, nodes or elements)
	integer lmax		!vertical values (return)
	character*(*) string	!data description
	integer ilhkv(np)	!number of layers in node k
	real hd(np)		!total depth in node k
	integer nlvddi		!vertical dimension of data
	real data(nlvddi,np)	!data
	integer ierr		!return error code

	logical b2d
	integer k,lm,l,npp
	real hdepth
	character*80 text

	ierr = 0
	npp = np
	lmax = nlvddi

	data = femflag

	if( iformat == 1 ) then
	  read(iunit,'(a)',err=13) text
	  if( nvers >= 3 ) read(iunit,*,err=11) npp,lmax
	  if( lmax > nlvddi ) goto 97
	  if( np /= npp ) goto 96
	  b2d = lmax .le. 1
	  if( b2d ) then
	    read(iunit,*,err=15) (data(1,k),k=1,np)
	  else
	    do k=1,np
	      read(iunit,*,err=15) lm,hd(k),(data(l,k),l=1,min(lm,lmax))
	      if( lm .gt. lmax ) goto 99
	      ilhkv(k) = lm
	    end do
	  end if
	else
	  read(iunit,err=13) text
	  if( nvers >= 3 ) read(iunit,err=11) npp,lmax
	  if( lmax > nlvddi ) goto 97
	  if( np /= npp ) goto 96
	  b2d = lmax .le. 1
	  if( b2d ) then
	    read(iunit,err=15) (data(1,k),k=1,np)
	  else
	    do k=1,np
	      read(iunit,err=15) lm,hd(k),(data(l,k),l=1,min(lm,lmax))
	      if( lm .gt. lmax ) goto 99
	      ilhkv(k) = lm
	    end do
	  end if
	end if

	if( b2d ) then
	  ilhkv = 1
	  hd = 10000.
	end if

	call fem_clean_data(nlvddi,np,ilhkv,data)	!delete nans etc..

	string = adjustl(text)

	return
   11	continue
	write(6,*) 'error reading data size'
	ierr = 11
	return
   13	continue
	write(6,*) 'error reading string description'
	ierr = 13
	return
   15	continue
	write(6,*) 'error reading data record'
	ierr = 15
	return
   96	continue
	write(6,*) 'error reading data record: size mismatch'
	write(6,*) 'np,npp: ',np,npp
	ierr = 96
	return
   97	continue
	write(6,*) 'error reading data record: lmax > nlvddi'
	write(6,*) 'lmax,nlvddi: ',lmax,nlvddi
	ierr = 97
	return
   99	continue
	write(6,*) 'error reading data record: too much vertical data'
	write(6,*) 'k,lm,lmax: ',k,lm,lmax
	ierr = 99
	return
	end

c************************************************************

	subroutine fem_file_skip_data(iformat,iunit
     +				,nvers,np,lmax
     +				,string,ierr)

c skips one record of data of the file

        implicit none

	integer iformat		!formatted or unformatted
	integer iunit		!file unit
	integer nvers		!version of file format
	integer np		!size of data (horizontal, nodes or elements)
	integer lmax		!vertical values
	character*(*) string	!data description
	integer ierr		!return error code

	logical b2d
	integer k,lm,l,npp
	real aux
	character*80 text

	ierr = 0
	npp = np

	if( iformat  == 1 ) then
	  read(iunit,'(a)',err=13) text
	  if( nvers >= 3 ) read(iunit,*,err=11) npp,lmax
	  if( np /= npp ) goto 96
	  b2d = lmax .le. 1
	  if( b2d ) then
	    read(iunit,*,err=15) (aux,k=1,np)
	  else
	    do k=1,np
	      read(iunit,*,err=15) lm,aux,(aux,l=1,lm)
	      if( lm .gt. lmax ) goto 99
	    end do
	  end if
	else
	  read(iunit,err=13) text
	  if( nvers >= 3 ) read(iunit,err=11) npp,lmax
	  if( np /= npp ) goto 96
	  b2d = lmax .le. 1
	  if( b2d ) then
	    read(iunit,err=15) (aux,k=1,np)
	  else
	    do k=1,np
	      read(iunit,err=15) lm,aux,(aux,l=1,lm)
	      if( lm .gt. lmax ) goto 99
	    end do
	  end if
	end if

	string = adjustl(text)

	return
   11	continue
	write(6,*) 'error reading data size'
	ierr = 11
	return
   13	continue
	write(6,*) 'error reading string description'
	ierr = 13
	return
   15	continue
	write(6,*) 'error skipping data record'
	ierr = 15
	return
   96	continue
	write(6,*) 'error reading data record: size mismatch'
	write(6,*) 'np,npp: ',np,npp
	ierr = 96
	return
   99	continue
	write(6,*) 'error reading data record: too much vertical data'
	write(6,*) 'k,lm,lmax: ',k,lm,lmax
	ierr = 99
	return
	end

c************************************************************

	subroutine fem_file_skip_record(iformat,iunit,atime,ierr)

	implicit none

	integer iformat		!formatted or unformatted
	integer iunit		!file unit
	double precision atime	!absolute time stamp
	integer ierr		!return error code

	integer nvers		!version of file format
	integer np		!size of data (horizontal, nodes or elements)
	integer lmax		!vertical values
	integer nvar		!number of variables to write
	integer ntype		!type of information contained
	real regpar(7)		!regular array params
	integer datetime(2)	!date and time information
	double precision dtime	!time stamp
	character*80 string	!data description

	integer iv

	call fem_file_read_params(iformat,iunit,dtime
     +				,nvers,np,lmax
     +				,nvar,ntype,datetime,ierr)
	if( ierr /= 0 ) return

	call fem_file_convert_time0(datetime,dtime,atime)

	call fem_file_skip_2header(iformat,iunit
     +				,ntype,lmax,regpar,ierr)
	if( ierr /= 0 ) return

	do iv=1,nvar
	  call fem_file_skip_data(iformat,iunit
     +				,nvers,np,lmax
     +				,string,ierr)
	  if( ierr /= 0 ) return
	end do

	end

c************************************************************

	subroutine fem_get_first_and_last(iformat,iunit
     +						,nrec,afirst,alast)

	implicit none

	integer iformat,iunit
	integer nrec
	double precision afirst,alast

	integer ierr
	double precision atime

	rewind(iunit)

	nrec = 0
	do
	  call fem_file_skip_record(iformat,iunit,atime,ierr)
	  if( ierr /= 0 ) exit
	  nrec = nrec + 1
	  if( nrec == 1 ) afirst = atime
	  alast = atime
	end do

	if( ierr > 0 ) then
	  write(6,*) 'ierr = ',ierr,'  nrec = ',nrec
	  stop 'error stop fem_get_first_and_last: read error'
	end if

	rewind(iunit)

	end

c************************************************************

	subroutine fem_clean_data(nlvddi,np,ilhkv,data)	

c delete nans etc..

	implicit none

	integer nlvddi,np
	integer ilhkv(np)
	real data(nlvddi,np)

	integer k,l,lmax
	logical fem_is_nan

	do k=1,np
	  lmax = ilhkv(k)
	  do l=1,lmax
	    if( fem_is_nan(data(l,k)) ) exit
	  end do
	  ilhkv(k) = l-1
	end do

	end

c************************************************************

        function fem_is_nan(val)

c tests val for NaN

        implicit none

        real val
        logical fem_is_nan
        integer itot

        itot = 0

        if( val .gt. 0. ) itot = itot + 1
        if( val .lt. 0. ) itot = itot + 1
        if( val .eq. 0. ) itot = itot + 1

        fem_is_nan = itot .ne. 1

	end

c************************************************************
c************************************************************
c************************************************************
c next three routines should be deleted 	FIXME
c************************************************************
c************************************************************
c************************************************************

	subroutine fem_file_string2time0(string,atime)

c converts string to time stamp
c
c string can be an absolute date (YYYY-MM-DD[::hh:mm:ss])
c or a relative time (integer)

	implicit none

	character*(*) string		!string indicating date
	double precision atime		!absolute time (return)

	integer year,month,day,hour,min,sec
	integer date,time,it
	integer ierr

	call dtsunform(year,month,day,hour,min,sec,string,ierr)

	if( ierr > 0 ) then
	  read(string,'(i10)',err=9) it
	  atime = it
	else
          call packdate(date,year,month,day)
          call packtime(time,hour,min,sec)
	  call dts_to_abs_time(date,time,atime)
	end if

	return
    9	continue
        write(6,*) '*** cannot parse date: ',ierr,string(1:20)
        write(6,*) '    format should be YYYY-MM-DD::hh:mm:ss'
        write(6,*) '    possible also YYYY-MM-DD[::hh[:mm[:ss]]]'
        write(6,*) '    or it should be an integer (relative time)'
	stop 'error stop fem_file_string2time: conversion error'
	end

c************************************************************

	subroutine fem_file_convert_time0(datetime,dtime,atime)

	implicit none

	integer datetime(2)		!reference date
	double precision dtime		!relative time
	double precision atime		!absolute time (return)

	double precision atime0

	if( datetime(1) > 0 ) then
	  call dts_to_abs_time(datetime(1),datetime(2),atime0)
	  atime = atime0 + dtime
	else
	  atime = dtime
	end if

	end

c************************************************************

	subroutine fem_file_convert_atime0(datetime,dtime,atime)

c converts from atime to datetime and dtime (only if in atime is real date)

	implicit none

	integer datetime(2)		!reference date (return)
	double precision dtime		!relative time (return)
	double precision atime		!absolute time (in)

	double precision atime1000	!1000 year limit
	parameter (atime1000 = 1000*365*86400.)

	if( atime > atime1000 ) then	!real date
	  call dts_from_abs_time(datetime(1),datetime(2),atime)
	  dtime = 0.
	else				!no date - keep relative time
	  datetime = 0
	  dtime = atime
	end if

	end

c************************************************************
c************************************************************
c************************************************************

	subroutine fem_file_make_type(ntype,imax,itype)
        integer ntype,imax,itype(imax)
        call fem_file_to_itype(ntype,imax,itype)
        end

c************************************************************

	subroutine fem_file_to_itype(ntype,imax,itype)

c ntype to itype

	implicit none

	integer ntype
	integer imax
	integer itype(imax)

	integer i,j

	j = ntype
	do i=1,imax
	  itype(i) = j - 10*(j/10)
	  j=j/10
	end do

	end

c************************************************************

	subroutine fem_file_from_itype(ntype,imax,itype)

c itype to ntype

	implicit none

	integer ntype
	integer imax
	integer itype(imax)

	integer i,j,ifact

	ntype = 0
        ifact = 1
	do i=1,imax
	  ntype = ntype + ifact*itype(i)
	  ifact=ifact*10
	end do

	end

c************************************************************

        subroutine fem_file_set_ntype(ntype,ipos,ival)

	integer ntype
	integer ipos
	integer ival

        integer imax
	integer, allocatable :: itype(:)

        imax = nint(log10(float(ntype)))
        imax = max(imax+2,ipos)
        allocate(itype(imax))
	call fem_file_to_itype(ntype,imax,itype)
        itype(ipos) = ival
	call fem_file_from_itype(ntype,imax,itype)

        end

c************************************************************

	function fem_file_regular(ntype)

	implicit none

	integer fem_file_regular
	integer ntype

	integer itype(2)

	call fem_file_make_type(ntype,2,itype)

	fem_file_regular = itype(2)

	end

c************************************************************

	subroutine fem_file_get_format_description(iformat,line)

	implicit none

	integer iformat
	character*(*) line

	if( iformat == 0 ) then
	  line = "unformatted"
	else if( iformat == 1 ) then
	  line = "formatted"
	else if( iformat == 2 ) then
	  line = "binary direct"
	else if( iformat == 3 ) then
	  line = "time series"
	else
	  line = "unknown"
	end if

	end

c************************************************************
c************************************************************
c************************************************************

	subroutine test_type

	integer imax,ntype
	integer itype(4)

	imax = 4
	ntype = 101
	call fem_file_make_type(ntype,imax,itype)
	write(6,*) ntype,itype
	ntype = 320
	call fem_file_make_type(ntype,imax,itype)
	write(6,*) ntype,itype
	end

c************************************************************
c	program subfile_main
c	call test_type
c	end
c************************************************************

