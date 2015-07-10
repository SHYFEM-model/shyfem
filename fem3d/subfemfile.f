c
c $Id: newpr.f,v 1.24 2010-02-22 15:38:36 georg Exp $
c
c administration of external files
c
c contents :
c
c revision log :
c
c 02.10.2012	ggu	created from scratch
c 16.05.2013	ggu	better documentation
c 24.04.2014	ggu	use nvar>0 as indication of good read
c 30.05.2014	ggu	restructured
c 16.06.2014	ggu	time is now double precision
c 07.07.2014	ggu	first version consolidated
c 20.10.2014	ggu	second version (date record is just after first record)
c 29.10.2014	ggu	new routine fem_file_is_fem_file()
c 09.01.2015	ggu	new routine fem_file_get_format_description()
c 14.01.2015	ggu	new routine fem_file_string2time()
c
c notes :
c
c format for file (nvers == 2)
c
c	time record 1
c	time record 2
c	time record ...
c
c format for time record
c
c	header record
c	data record for variable 1
c	data record for variable 2
c	data record for variable ...
c	data record for variable nvar
c
c format for header record
c
c	dtime,nvers,id,np,lmax,nvar,ntype
c	date,time				only if ntype odd
c	(hlv(l),l=1,lmax)			only if( lmax > 1 )
c	other lines depending on ntype
c
c format for data record
c
c	if( lmax == 1 )
c		string
c		(data(1,k),k=1,np)
c	if( lmax > 1 )
c		string
c		do k=1,np
c		  lm,hd(k),(data(l,k),l=1,lm)
c		end do
c
c format for regpar
c
c	nx,ny,x0,y0,dx,dy,flag
c
c legend
c
c dtime		time stamp (double precision, seconds)
c nvers		version of file format
c id		id to identify fem file (must be 957839)
c np		number of horizontal points given
c lmax		maximum number of layers given
c nvar		number of variables in time record
c ntype		type of data, defines extra data to follow
c hlv		layer depths
c string 	string with description of data
c hd(k)		total depth in node k
c data(l,k)	data for variable
c lm		total number of vertical data provided for point k
c k,l		index for horizontal/vertical dimension
c nx,ny		size of regular grid
c x0,y0		origin of regular grid
c dx,dy		space increment of regular grid
c flag		flag for invalid data of regular grid
c
c file type (ntype)
c
c 0		no other lines in header
c 1		give date/time of reference in extra line
c 10		regular grid, information on extra line
c 20		rotated regular grid, information on extra line (not yet ready)
c
c combinations are possible, example:
c
c 21		date/time and regular rotated grid
c
c************************************************************
c************************************************************
c************************************************************

	subroutine fem_file_write_header(iformat,iunit,it
     +				,nvers,np,lmax
     +				,nvar,ntype
     +				,nlvddi,hlv,datetime,regpar)

c writes header of fem file

        implicit none

	integer iformat		!formatted or unformatted
	integer iunit		!file unit
	double precision it	!time stamp
	integer nvers		!version of file format
	integer np		!size of data (horizontal, nodes or elements)
	integer lmax		!maximum vertical values (1 for 2d)
	integer nvar		!number of variables to write
	integer ntype		!type of information contained
	integer nlvddi		!vertical dimension of data
	real hlv(nlvddi)	!depth at bottom of layer
	integer datetime(2)	!date and time parameters
	real regpar(7)		!parameters for regular field

	call fem_file_write_params(iformat,iunit,it
     +				,nvers,np,lmax
     +				,nvar,ntype,datetime)

	call fem_file_write_2header(iformat,iunit,ntype,lmax
     +			,hlv,regpar)

	end

c************************************************************

	subroutine fem_file_write_params(iformat,iunit,it
     +				,nvers,np,lmax
     +				,nvar,ntype,datetime)

c writes first header of fem file

        implicit none

	integer idfem
	parameter ( idfem = 957839 )

	integer iformat		!formatted or unformatted
	integer iunit		!file unit
	double precision it	!time stamp
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
	if( nv .eq. 0 ) nv = 2	!default
	if( nv .lt. 2 .or. nv .gt. 2 ) goto 99
	if( lmax < 1 ) goto 98

	if( iformat == 1 ) then
	  itlong = it
	  write(iunit,1000) itlong,nv,idfem,np,lmax,nvar,ntype
	else
	  write(iunit) it,nv,idfem,np,lmax,nvar,ntype
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
	    write(iunit,*) (hlv(l),l=1,lmax)
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

	end

c************************************************************

	subroutine fem_file_write_data(iformat,iunit
     +				,nvers,np,lmax
     +				,string
     +				,ilhkv,hd
     +				,nlvddi,data)

c writes data of the file

        implicit none

	integer iformat		!formatted or unformatted
	integer iunit		!file unit
	integer nvers		!version of file format
	integer np		!size of data (horizontal, nodes or elements)
	integer lmax		!maximum vertical values (1 for 2d)
	character*(*) string	!string explanation
	integer ilhkv(np)	!number of layers in point k (node)
	real hd(np)		!total depth
	integer nlvddi		!vertical dimension of data
	real data(nlvddi,np)	!data

	logical b2d
	integer k,lm,l,nv
	character*60 text
	character*80 textu	!we need 80 chars for unformatted write

	nv = nvers
	if( nv .eq. 0 ) nv = 1	!default

	text = string
	textu = string
	b2d = lmax .le. 1

	if( iformat == 1 ) then
	  write(iunit,*) text
	  if( b2d ) then
	    write(iunit,*) (data(1,k),k=1,np)
	  else
	    do k=1,np
	      lm = ilhkv(k)
	      write(iunit,*) lm,hd(k),(data(l,k),l=1,lm)
	    end do
	  end if
	else
	  write(iunit) textu
	  if( b2d ) then
	    write(iunit) (data(1,k),k=1,np)
	  else
	    do k=1,np
	      lm = ilhkv(k)
	      write(iunit) lm,hd(k),(data(l,k),l=1,lm)
	    end do
	  end if
	end if

	end

c************************************************************
c************************************************************
c************************************************************

	subroutine fem_file_is_fem_file(file,iformat)

c tries to open fem file for read
c
c returns -1 in iformat if no fem file, else iformat indicates format

	implicit none

	character*(*) file	!file name
	integer iformat		!is formatted? -1 for no fem file (return)

	integer nvar,np,ntype
	logical filex

	iformat = 0

	if( .not. filex(file) ) then
	  write(6,*) 'file does not exist: ',file
	  return
	end if

	call fem_file_test_formatted(file,np,nvar,ntype,iformat)

	end

c************************************************************

	subroutine fem_file_read_open(file,nexp,iunit,iformat)

c opens fem file for read

	implicit none

	character*(*) file	!file name
	integer nexp		!expected size of data (0 if unknown)
	integer iunit		!unit of opened file (return) (0 for error)
	integer iformat		!is formatted? (return)

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

	implicit none

	character*(*) file	!file name
	integer iformat		!is formatted?

	integer idfem
	parameter ( idfem = 957839 )

	integer iunit
	integer it,nvers,np,lmax,nvar,ntype
	integer id
	integer ios

	it = 0
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
	  read(iunit,*,iostat=ios) it,nvers,id,np,lmax,nvar,ntype
	  write(6,*) 'formatted read: '
	else
	  open(iunit,file=file,form='unformatted',status='old')
	  read(iunit,iostat=ios) it,nvers,id,np,lmax,nvar,ntype
	  write(6,*) 'unformatted read: '
	end if

	if( ios .gt. 0 ) then
	  write(6,*) 'fem_file_write_info: error reading file'
	else if( ios .lt. 0 ) then
	  write(6,*) 'fem_file_write_info: EOF found'
	else if( id .ne. idfem ) then
	  write(6,*) 'file is not a FEM file: ',id,idfem
	else
	  write(6,*) it,nvers,id,np,lmax,nvar,ntype
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
	double precision it
	integer ierr
	logical bdebug

c------------------------------------------------------
c initialize parameters
c------------------------------------------------------

	bdebug = .true.
	bdebug = .false.

	nvers = 0
	it = 0
	np0 = 0
	lmax = 0
	nvar = 0
	ntype = 0

c------------------------------------------------------
c find unit to open file
c------------------------------------------------------

	iunit = 90
	call find_unit(iunit)

c------------------------------------------------------
c first try unformatted
c------------------------------------------------------

	ierr = 77
	open(iunit,file=file,form='unformatted',status='old',err=2)

	iformat = 0
	call fem_file_read_params(iformat,iunit,it
     +				,nvers,np,lmax,nvar,ntype,datetime,ierr)

	close(iunit)

	if( ierr .ne. 0 ) then
	  if( bdebug ) write(6,*) 'unformatted read error ',ierr
	else	!ok, probably unformatted
	  return
	end if

    2	continue

c------------------------------------------------------
c now try formatted
c------------------------------------------------------

	ierr = 77
	open(iunit,file=file,form='formatted',status='old',err=8)

	iformat = 1
	call fem_file_read_params(iformat,iunit,it
     +				,nvers,np,lmax,nvar,ntype,datetime,ierr)

	close(iunit)

	if( ierr .ne. 0 ) then
	  if( bdebug ) write(6,*) 'formatted read error ',ierr
	else	!ok, probably formatted
	  return
	end if

    8	continue

c------------------------------------------------------
c no successful opening
c------------------------------------------------------

	np = 0
	nvar = 0
	ntype = 0
	iformat = -1
	iformat = -ierr

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
	character*80 strings(1)		!return - must have dimension nvar
	integer ierr

	integer iformat
	integer np0,iunit,i
	integer nvers,np,lmax,nvar,ntype
	integer datetime(2)
	double precision it
	character*80 string

	np0 = 0
	ierr = 1

	call fem_file_read_open(file,np0,iunit,iformat)
	if( iunit .le. 0 ) return

	call fem_file_read_params(iformat,iunit,it
     +				,nvers,np,lmax
     +				,nvar,ntype,datetime,ierr)
	if( ierr .ne. 0 ) return

	call fem_file_skip_2header(iformat,iunit
     +				,ntype,lmax,ierr)
	if( ierr .ne. 0 ) return

	do i=1,nvar
	  call fem_file_skip_data(iformat,iunit
     +				,nvers,np,lmax
     +				,string,ierr)
	  if( ierr .ne. 0 ) return
	  strings(i) = string
	end do

	close(iunit)
	ierr = 0

	end

c************************************************************
c************************************************************
c************************************************************

	subroutine fem_file_read_params(iformat,iunit,it
     +				,nvers,np,lmax
     +				,nvar,ntype,datetime,ierr)

c reads and checks params of next header

        implicit none

	integer iformat		!formatted or unformatted
	integer iunit		!file unit
	double precision it	!time stamp
	integer nvers		!version of file format
	integer np		!size of data (horizontal, nodes or elements)
	integer lmax		!vertical values
	integer nvar		!number of variables to write
	integer ntype		!type of information contained
	integer datetime(2)	!date and time information
	integer ierr		!return error code

	integer id,i
	integer itype(2)

	ierr = 0
	if( iunit < 1 ) goto 99

	if( iformat == 1 ) then
	  read(iunit,*,end=1,err=2) it,nvers,id,np,lmax,nvar,ntype
	else
	  read(iunit,end=1,err=2) it,nvers,id,np,lmax,nvar,ntype
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
	else
	  datetime = 0
	end if

	return
    1	continue
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

	subroutine fem_file_peek_params(iformat,iunit,it
     +				,nvers,np,lmax,nvar,ntype,datetime,ierr)

c reads and checks params of next header (non advancing read)

        implicit none

	integer iformat		!formatted or unformatted
	integer iunit		!file unit
	double precision it	!time stamp
	integer nvers		!version of file format
	integer np		!size of data (horizontal, nodes or elements)
	integer lmax		!vertical values
	integer nvar		!number of variables to write
	integer ntype		!type of information contained
	integer datetime(2)	!date and time information
	integer ierr		!return error code

	integer itype(2)

	call fem_file_read_params(iformat,iunit,it
     +				,nvers,np,lmax,nvar,ntype,datetime,ierr)

	if( ierr .ne. 0 ) return

	backspace(iunit)

	call fem_file_make_type(ntype,2,itype)
	if( itype(1) > 0 ) backspace(iunit)

	end

c************************************************************

	subroutine fem_file_check_params(nvers,id,np,lmax,nvar,ntype,ierr)

c reads and checks params of next header

        implicit none

	integer nvers		!version of file format
	integer id		!id of fem file
	integer np		!size of data (horizontal, nodes or elements)
	integer lmax		!vertical values
	integer nvar		!number of variables to write
	integer ntype		!type of information contained
	integer ierr		!return error code

	integer idfem
	parameter ( idfem = 957839 )

	ierr = 11
	if( id .ne. idfem ) goto 9
	ierr = 13
	if( nvers .lt. 1 .or. nvers .gt. 2 ) goto 9
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

        implicit none

	integer iformat		!formatted or unformatted
	integer iunit		!file unit
	integer ntype		!type of second header
	integer lmax		!total number of elements to read
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
	else
	  regpar = 0.
	end if

	ierr = 0
	return
    1	continue
	end

c************************************************************

	subroutine fem_file_skip_2header(iformat,iunit
     +				,ntype,lmax,ierr)

c skips additional headers in fem file

        implicit none

	integer iformat		!formatted or unformatted
	integer iunit		!file unit
	integer lmax		!total number of elements to read
	integer ntype		!type of second header
	integer ierr		!return error code

	integer l,i
	integer itype(2)
	real aux

	call fem_file_make_type(ntype,2,itype)

	ierr = 3
	if( lmax .gt. 1 ) then
	  if( iformat  == 1 ) then
	    read(iunit,*,err=1) (aux,l=1,lmax)
	  else
	    read(iunit,err=1) (aux,l=1,lmax)
	  end if
	end if

	ierr = 7
	if( itype(2) .gt. 0 ) then
	  if( iformat == 1 ) then
	    read(iunit,*,err=1) (aux,i=1,7)
	  else
	    read(iunit,err=1) (aux,i=1,7)
	  end if
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

        implicit none

	integer iformat		!formatted or unformatted
	integer iunit		!file unit
	integer nvers		!version of file format
	integer np		!size of data (horizontal, nodes or elements)
	integer lmax		!vertical values
	character*(*) string	!string explanation
	integer ilhkv(np)	!number of layers in point k (node)
	real hd(np)		!total depth
	integer nlvddi		!vertical dimension of data
	real data(nlvddi,np)	!data
	integer ierr		!return error code

	logical b2d
	integer k,lm,l
	real hdepth
	character*80 text

	ierr = 0
	b2d = lmax .le. 1

	if( iformat == 1 ) then
	  read(iunit,'(a)',err=13) text
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
	  do k=1,np
	    ilhkv(k) = 1
	    hd(k) = 10000.
	  end do
	end if

	string = text

	return
   13	continue
	write(6,*) 'error reading string description'
	ierr = 13
	return
   15	continue
	write(6,*) 'error reading data record'
	ierr = 15
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
	character*(*) string	!string explanation
	integer ierr		!return error code

	logical b2d
	integer k,lm,l
	real aux
	character*80 text

	ierr = 0
	b2d = lmax .le. 1

	if( iformat  == 1 ) then
	  read(iunit,'(a)',err=13) text
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
	  if( b2d ) then
	    read(iunit,err=15) (aux,k=1,np)
	  else
	    do k=1,np
	      read(iunit,err=15) lm,aux,(aux,l=1,lm)
	      if( lm .gt. lmax ) goto 99
	    end do
	  end if
	end if

	string = text

	return
   13	continue
	write(6,*) 'error reading string description'
	ierr = 13
	return
   15	continue
	write(6,*) 'error skipping data record'
	ierr = 15
	return
   99	continue
	write(6,*) 'error reading data record: too much vertical data'
	write(6,*) 'k,lm,lmax: ',k,lm,lmax
	ierr = 99
	return
	end

c************************************************************
c************************************************************
c************************************************************

	subroutine fem_file_string2time(string,atime)

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

	subroutine fem_file_convert_time(datetime,dtime,atime)

	implicit none

	integer datetime(2)		!reference date
	double precision dtime		!relative time
	double precision atime		!absolute time (return)

	double precision dtime0

	if( datetime(1) > 0 ) then
	  call dts_to_abs_time(datetime(1),datetime(2),dtime0)
	  atime = dtime0 + dtime
	else
	  atime = dtime
	end if

	end

c************************************************************

	subroutine fem_file_convert_atime(datetime,dtime,atime)

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

	subroutine fem_file_make_type(ntype,imax,itype)

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

	function fem_file_regular(ntype)

	implicit none

	integer ntype

	integer fem_file_regular
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






