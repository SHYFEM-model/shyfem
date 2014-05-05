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
c
c notes :
c
c format for file (nvers == 1)
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
c	it,nvers,np,lmax,nvar,ntype
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
c legend
c
c it		time stamp (integer, seconds)
c nvers		version of file format
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
c
c file type (ntype)
c
c 0		no other lines in header
c 1		give date/time of reference in extra line (not yet ready)
c 10		regular grid, information on extra line (not yet ready)
c 20		rotated regular grid, information on extra line (not yet ready)
c
c combinations are possible, example:
c
c 21		date/time and regular rotated grid
c
c************************************************************
c************************************************************
c************************************************************

	subroutine fem_file_write_header(bformat,iunit,it
     +				,nvers,np,lmax,nvar,ntype,nlvdim,hlv)

c writes header of the file

        implicit none

	logical bformat		!formatted or unformatted
	integer iunit		!file unit
	integer it		!time stamp
	integer nvers		!version of file format
	integer np		!size of data (horizontal, nodes or elements)
	integer lmax		!maximum vertical values
	integer nvar		!number of variables to write
	integer ntype		!type of information contained
	integer nlvdim		!vertical dimension of data
	real hlv(nlvdim)	!depth at bottom of layer

	integer l,nv

	nv = nvers
	if( nv .eq. 0 ) nv = 1	!default
	if( nv .lt. 1 .or. nv .gt. 1 ) goto 99

	if( bformat ) then
	  write(iunit,1000) it,nv,np,lmax,nvar,ntype
	  if( lmax .gt. 1 ) write(iunit,*) (hlv(l),l=1,lmax)
	else
	  write(iunit) it,nv,np,lmax,nvar,ntype
	  if( lmax .gt. 1 ) write(iunit) (hlv(l),l=1,lmax)
	end if

	return
 1000	format(i14,5i12)
   99	continue
	write(6,*) 'nvers = ',nvers
	stop 'error stop fem_file_write_header: nvers'
	end

c************************************************************

	subroutine fem_file_write_data(bformat,iunit
     +				,nvers,np,lmax
     +				,nlvdim,ilhkv
     +				,string,hd,data)

c writes data of the file

        implicit none

	logical bformat		!formatted or unformatted
	integer iunit		!file unit
	integer nvers		!version of file format
	integer np		!size of data (horizontal, nodes or elements)
	integer lmax		!maximum vertical values (1 for 2d)
	integer nlvdim		!vertical dimension of data
	integer ilhkv(np)	!number of layers in point k (node)
	character*(*) string	!string explanation
	real hd(np)		!total depth
	real data(nlvdim,np)	!data

	logical b2d
	integer k,lm,l,nv
	character*60 text
	character*80 textu	!we need 80 chars for unformatted write

	nv = nvers
	if( nv .eq. 0 ) nv = 1	!default

	text = string
	textu = string
	b2d = lmax .le. 1

	if( bformat ) then
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

	subroutine fem_file_write_3d(bformat,iunit,it
     +				,np,lmax,nlvdim,hlv
     +				,ilhkv,string,hd,data)

c writes 1 variable of a 3d field

	implicit none

	logical bformat		!formatted or unformatted
	integer iunit		!file unit
	integer it		!time stamp
	integer np		!size of data (horizontal, nodes or elements)
	integer lmax		!vertical values
	integer nlvdim		!vertical dimension
	real hlv(nlvdim)	!depth at bottom of layer
	integer ilhkv(np)	!number of layers in point k (node)
	character*(*) string	!string explanation
	real hd(np)		!total depth
	real data(nlvdim,np)	!data

	integer nvers,ntype
	integer nvar

	nvers = 1
	ntype = 0
	nvar = 1

	call fem_file_write_header(bformat,iunit,it,nvers
     +				,np,lmax,nvar,ntype,nlvdim,hlv)
	call fem_file_write_data(bformat,iunit,nvers
     +				,np,lmax,nlvdim,ilhkv,string,hd,data)

	end

c************************************************************

	subroutine fem_file_write_2d(bformat,iunit,it
     +				,np,string,data)

c writes 1 variable of a 2d field

	implicit none

	logical bformat		!formatted or unformatted
	integer iunit		!file unit
	integer it		!time stamp
	integer np		!size of data (horizontal, nodes or elements)
	character*(*) string	!string explanation
	real data(np)		!data

	integer nvers,ntype
	integer nvar,lmax,nlvdim
	real hlv(1),hd(1)
	integer ilhkv(1)

	nvers = 1
	ntype = 0
	nvar = 1
	lmax = 1
	nlvdim = 1
	hlv(1) = 10000.
	hd(1) = 10000.
	ilhkv(1) = 1

	call fem_file_write_header(bformat,iunit,it,nvers
     +				,np,lmax,nvar,ntype,nlvdim,hlv)
	call fem_file_write_data(bformat,iunit,nvers
     +				,np,lmax,nlvdim,ilhkv,string,hd,data)

	end

c************************************************************
c************************************************************
c************************************************************

	subroutine fem_file_read_open(name,np,iunit,bformat)

c opens file for read

	implicit none

	character*(*) name	!file name
	integer np		!expected size of data
	integer iunit		!unit of opened file (in/out) (0 for error)
	logical bformat		!is formatted? (return)

	integer nvar
	logical filex

	iunit = 0
	bformat = .true.

	if( .not. filex(name) ) then
	  write(6,*) 'file does not exist: ',name
	  return
	end if

	call fem_file_test_formatted(name,np,nvar,bformat)

	if( nvar .gt. 0 ) then
	  call find_unit(iunit)
	  if( bformat ) then
	    open(iunit,file=name,form='formatted',status='old')
	  else
	    open(iunit,file=name,form='unformatted',status='old')
	  end if
	else
	  write(6,*) 'fem_file_read_open: error opening file ',name
	  call fem_file_write_info(name,.true.)
	  call fem_file_write_info(name,.false.)
	end if

	end

c************************************************************

	subroutine fem_file_write_info(name,bformat)

c writes information on file from header

	character*(*) name	!file name
	logical bformat		!is formatted?

	integer iunit
	integer it,nvers,np,lmax,nvar,ntype
	integer ios

	it = 0
	np = 0
	lmax = 0
	nvar = 0
	nvers = 0
	ntype = 0

	iunit = 90
	call find_unit(iunit)

	if( bformat ) then
	  open(iunit,file=name,form='formatted',status='old')
	  read(iunit,*,iostat=ios) it,nvers,np,lmax,nvar,ntype
	  write(6,*) 'formatted read: '
	else
	  open(iunit,file=name,form='unformatted',status='old')
	  read(iunit,iostat=ios) it,nvers,np,lmax,nvar,ntype
	  write(6,*) 'unformatted read: '
	end if

	if( ios .gt. 0 ) then
	  write(6,*) 'fem_file_write_info: error reading file'
	else if( ios .lt. 0 ) then
	  write(6,*) 'fem_file_write_info: EOF found'
	else
	  write(6,*) it,nvers,np,lmax,nvar,ntype
	end if

	close(iunit)

	end

c************************************************************

	subroutine fem_file_test_formatted(name,np,nvar,bformat)

c checks if file is readable and formatted or unformatted

	implicit none

	character*(*) name	!file name
	integer np		!expected size of data, 0 if no idea
	integer nvar		!successful read => nvar>0 (return)
	logical bformat		!is formatted? (return)

	integer iunit
	integer it,nvers,np0,lmax,ntype
	integer ierr
	logical bdebug

c------------------------------------------------------
c initialize parameters
c------------------------------------------------------

	bdebug = .false.
	bdebug = .true.

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

	open(iunit,file=name,form='unformatted',status='old',err=2)

	bformat = .false.
	call fem_file_get_params(bformat,iunit,it
     +				,nvers,np0,lmax,nvar,ntype,ierr)

	close(iunit)

	if( ierr .ne. 0 ) then
	  if( bdebug ) write(6,*) 'unformatted read error'
	else if( np .gt. 0 .and. np0 .gt. 1 .and. np .ne. np0 ) then
	  if( bdebug ) then
	    write(6,*) 'unformatted error in np,np0: ',np,np0
	  end if
	else if( nvar .le. 0 .or. lmax .lt. 0 ) then
	  if( bdebug ) then
	    write(6,*) 'unformatted error in nvar,lmax: ',nvar,lmax
	  end if
	else	!ok, probably unformatted
	  return
	end if

    2	continue

c------------------------------------------------------
c now try formatted
c------------------------------------------------------

	open(iunit,file=name,form='formatted',status='old',err=8)

	bformat = .true.
	call fem_file_get_params(bformat,iunit,it
     +				,nvers,np0,lmax,nvar,ntype,ierr)

	close(iunit)

	if( ierr .ne. 0 ) then
	  if( bdebug ) write(6,*) 'formatted read error'
	else if( np .gt. 0 .and. np0 .gt. 1 .and. np .ne. np0 ) then
	  if( bdebug ) then
	    write(6,*) 'formatted error in np,np0: ',np,np0
	  end if
	else if( nvar .le. 0 .or. lmax .lt. 0 ) then
	  if( bdebug ) then
	    write(6,*) 'formatted error in nvar,lmax: ',nvar,lmax
	  end if
	else	!ok, probably formatted
	  return
	end if

    8	continue

c------------------------------------------------------
c no successful opening
c------------------------------------------------------

	nvar = 0

c------------------------------------------------------
c end of routine
c------------------------------------------------------

	end

c************************************************************

	subroutine fem_file_get_params(bformat,iunit,it
     +				,nvers,np,lmax,nvar,ntype,ierr)

c reads and checks params of next header

        implicit none

	logical bformat		!formatted or unformatted
	integer iunit		!file unit
	integer it		!time stamp
	integer nvers		!version of file format
	integer np		!size of data (horizontal, nodes or elements)
	integer lmax		!vertical values
	integer nvar		!number of variables to write
	integer ntype		!type of information contained
	integer ierr		!return error code

	ierr = 0

	if( bformat ) then
	  read(iunit,*,end=1,err=2) it,nvers,np,lmax,nvar,ntype
	else
	  read(iunit,end=1,err=2) it,nvers,np,lmax,nvar,ntype
	end if

	call fem_file_check_params(nvers,np,lmax,nvar,ntype,ierr)
	if( ierr .ne. 0 ) goto 9

	ierr = 0
	backspace(iunit)
	return

    9	continue
	write(6,*) 'fem_file_get_params problem gguuu ',bformat
	ierr = 9
	!backspace(iunit)
	return

    1	continue
	ierr = -1
	backspace(iunit)
	return

    2	continue
	ierr = 1
	backspace(iunit)
	return

	end

c************************************************************

	subroutine fem_file_get_hlv(bformat,iunit,nlvdim,nlv,hlv,ierr)

c reads hlv of next header

        implicit none

	logical bformat		!formatted or unformatted
	integer iunit		!file unit
	integer nlvdim		!dimensions of vertical array
	integer nlv		!number of levels
	real hlv(nlvdim)	!vertical structure
	integer ierr		!return error code

	integer it,nvers,np,lmax,nvar,ntype
	integer l

	ierr = 0

	if( bformat ) then
	  read(iunit,*,end=1,err=2) it,nvers,np,lmax,nvar,ntype
	else
	  read(iunit,end=1,err=2) it,nvers,np,lmax,nvar,ntype
	end if

	if( lmax .gt. nlvdim ) goto 9

	if( lmax .gt. 1 ) then
	  if( bformat ) then
	    read(iunit,*,err=3) (hlv(l),l=1,lmax)
	  else
	    read(iunit,err=3) (hlv(l),l=1,lmax)
	  end if
	end if

	nlv = lmax
	if( lmax .eq. 1 ) hlv(1) = 10000.

	ierr = 0
	backspace(iunit)
	if( lmax .gt. 1 ) backspace(iunit)
	return

    1	continue
	ierr = -1
	backspace(iunit)
	return

    2	continue
	ierr = 1
	backspace(iunit)
	return

    3	continue
	ierr = 3
	backspace(iunit)
	backspace(iunit)
	return

    9	continue
	ierr = 9
	backspace(iunit)
	return

	end

c************************************************************

	subroutine fem_file_check_params(nvers,np,lmax,nvar,ntype,ierr)

c reads and checks params of next header

        implicit none

	integer nvers		!version of file format
	integer np		!size of data (horizontal, nodes or elements)
	integer lmax		!vertical values
	integer nvar		!number of variables to write
	integer ntype		!type of information contained
	integer ierr		!return error code

	ierr = 0

	if( nvers .lt. 1 .or. nvers .gt. 1 ) goto 9
	if( np .le. 0 ) goto 9
	if( lmax .le. 0 .or. lmax .gt. 1000 ) goto 9
	if( nvar .le. 0 .or. nvar .gt. 100 ) goto 9
	if( ntype .lt. 0 .or. ntype .gt. 2 ) goto 9

	return

    9	continue
	ierr = 9

	return
	end

c************************************************************

	function fem_std_func(lmax,hlv,lm_data,hd_data,data_in
     +			,nlvdim,lm_fem,hd_fem,data_out)

c standard function for copying data to new structure

	implicit none

	logical fem_std_func
	integer lmax		!max layers in file
	real hlv(lmax)		!layer depth in file
	integer lm_data		!vertical levels read for point
	real hd_data		!depth value read for point
	real data_in(lm_data)	!data read for point
	integer nlvdim		!vertical dimension of data_out
	integer lm_fem		!number of levels of converted data
	real hd_fem		!depth of converted data
	real data_out(nlvdim)	!converted data

	integer l,lin

	lin = min(lm_data,nlvdim)

	do l=1,lin
	  data_out(l) = data_in(l)
	end do

	do l=lin+1,nlvdim
	  data_out(l) = data_in(lin)
	end do

	lm_fem = lm_data
	hd_fem = hd_data

	fem_std_func = .true.

	end

c************************************************************
c************************************************************
c************************************************************

	subroutine fem_file_read_header(bformat,iunit,it
     +				,nvers,np,lmax,nvar,ntype,nlvdim,hlv,ierr)

c reads header of the file

        implicit none

	logical bformat		!formatted or unformatted
	integer iunit		!file unit
	integer it		!time stamp
	integer nvers		!version of file format
	integer np		!size of data (horizontal, nodes or elements)
	integer lmax		!vertical values
	integer nvar		!number of variables to write
	integer ntype		!type of information contained
	integer nlvdim		!vertical dimension of data
	real hlv(nlvdim)	!depth at bottom of layer
	integer ierr		!return error code

	integer l,np0

	ierr = -1
	hlv(1) = 10000.

	if( bformat ) then
	  read(iunit,*,end=1) it,nvers,np0,lmax,nvar,ntype
	  if( lmax .gt. 1 ) read(iunit,*) (hlv(l),l=1,min(lmax,nlvdim))
	else
	  read(iunit,end=1) it,nvers,np0,lmax,nvar,ntype
	  if( lmax .gt. 1 ) read(iunit) (hlv(l),l=1,min(lmax,nlvdim))
	end if

	if( lmax .gt. nlvdim) goto 98

	call fem_file_check_params(nvers,np0,lmax,nvar,ntype,ierr)
	if( ierr .ne. 0 ) goto 99
	if( np .gt. 0 .and. np0 .gt. 1 .and. np .ne. np0 ) goto 99

	np = np0
	ierr = 0

    1	continue

	return
   98	continue
	write(6,*) 'nlvdim,lmax: ',nlvdim,lmax
	backspace(iunit)
	ierr = 98
	return
	!stop 'error stop fem_file_read_header: dimension nlvdim'
   99	continue
	write(6,*) it,nvers,np0,lmax,nvar,ntype
	write(6,*) np,np0
	stop 'error stop fem_file_read_header: parameters'
	end

c************************************************************

	subroutine fem_file_read_data(bformat,iunit
     +				,nvers,np,lmax
     +				,ilhkv,hd
     +				,string,nlvdim,data)

c reads data of the file

        implicit none

	logical bformat		!formatted or unformatted
	integer iunit		!file unit
	integer nvers		!version of file format
	integer np		!size of data (horizontal, nodes or elements)
	integer lmax		!vertical values
	integer ilhkv(np)	!number of layers in point k (node)
	real hd(np)		!total depth
	character*(*) string	!string explanation
	integer nlvdim		!vertical dimension of data
	real data(nlvdim,np)	!data

	logical b2d
	integer k,lm,l
	real hdepth
	character*80 text

	b2d = lmax .le. 1

	if( bformat ) then
	  read(iunit,'(a)') text
	  if( b2d ) then
	    read(iunit,*) (data(1,k),k=1,np)
	  else
	    do k=1,np
	      read(iunit,*) lm,hd(k),(data(l,k),l=1,min(lm,lmax))
	      if( lm .gt. lmax ) goto 99
	      ilhkv(k) = lm
	    end do
	  end if
	else
	  read(iunit) text
	  if( b2d ) then
	    read(iunit) (data(1,k),k=1,np)
	  else
	    do k=1,np
	      read(iunit) lm,hd(k),(data(l,k),l=1,min(lm,lmax))
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
   99	continue
	write(6,*) 'k,lm,lmax: ',k,lm,lmax
	stop 'error stop fem_file_read_data: dimension lmax'
	end

c************************************************************

	subroutine fem_file_skip_data(bformat,iunit
     +				,nvers,np,lmax
     +				,string)

c skips one record of data of the file

        implicit none

	logical bformat		!formatted or unformatted
	integer iunit		!file unit
	integer nvers		!version of file format
	integer np		!size of data (horizontal, nodes or elements)
	integer lmax		!vertical values
	character*(*) string	!string explanation

	logical b2d
	integer k,lm,l
	real aux
	character*80 text

	b2d = lmax .le. 1

	if( bformat ) then
	  read(iunit,'(a)') text
	  if( b2d ) then
	    read(iunit,*) (aux,k=1,np)
	  else
	    do k=1,np
	      read(iunit,*) lm,aux,(aux,l=1,lm)
	      if( lm .gt. lmax ) goto 99
	    end do
	  end if
	else
	  read(iunit) text
	  if( b2d ) then
	    read(iunit) (aux,k=1,np)
	  else
	    do k=1,np
	      read(iunit) lm,aux,(aux,l=1,lm)
	      if( lm .gt. lmax ) goto 99
	    end do
	  end if
	end if

	string = text

	return
   99	continue
	write(6,*) 'k,lm,lmax: ',k,lm,lmax
	stop 'error stop fem_file_read_data: dimension lmax'
	end

c************************************************************

	subroutine fem_file_read_3d(bformat,iunit,it
     +				,np,lmax,nlvdim,hlv
     +				,ilhkv,string,hd,data,ierr)

c reads 1 variable of a 3d field

	implicit none

	logical bformat		!formatted or unformatted
	integer iunit		!file unit
	integer it		!time stamp
	integer np		!size of data (horizontal, nodes or elements)
	integer lmax		!vertical values
	integer nlvdim		!vertical dimension
	real hlv(nlvdim)	!depth at bottom of layer
	integer ilhkv(np)	!number of layers in point k (node)
	character*(*) string	!string explanation
	real hd(np)		!total depth
	real data(nlvdim,np)	!data
	integer ierr		!return error code

	integer nvers,nvar,ntype

	call fem_file_read_header(bformat,iunit,it
     +				,nvers,np,lmax,nvar,ntype,nlvdim,hlv,ierr)

	if( ierr .ne. 0 ) return
	if( nvar .ne. 1 ) goto 99

	call fem_file_read_data(bformat,iunit
     +				,nvers,np,lmax
     +				,ilhkv,hd
     +				,string,nlvdim,data)

	return
   99	continue
	write(6,*) 'nvar: ',nvar
	stop 'error stop fem_file_read_3d: can read only one variable'
	end

c************************************************************

	subroutine fem_file_read_2d(bformat,iunit,it
     +				,np,string,data,ierr)

c reads 1 variable of a 2d field

	implicit none

	logical bformat		!formatted or unformatted
	integer iunit		!file unit
	integer it		!time stamp
	integer np		!size of data (horizontal, nodes or elements)
	character*(*) string	!string explanation
	real data(np)		!data
	integer ierr		!return error code

	integer nvers,lmax,nvar,ntype,nlvdim
	real hlv(1),hd(1)
	integer ilhkv(1)

	logical fem_std_func
	external fem_std_func

	nlvdim = 1
	hlv(1) = 10000.
	hd(1) = 10000.
	ilhkv(1) = 1

	call fem_file_read_header(bformat,iunit,it
     +				,nvers,np,lmax,nvar,ntype,nlvdim,hlv,ierr)

	if( ierr .ne. 0 ) return
	if( nvar .ne. 1 ) goto 99
	if( lmax .ne. 1 ) goto 98

	call fem_file_read_data(bformat,iunit
     +				,nvers,np,lmax
     +				,ilhkv,hd
     +				,string,nlvdim,data)

	return
   98	continue
	write(6,*) 'lmax: ',lmax
	stop 'error stop fem_file_read_2d: data is not 2d'
   99	continue
	write(6,*) 'nvar: ',nvar
	stop 'error stop fem_file_read_2d: can read only one variable'
	end

c************************************************************
c************************************************************
c************************************************************

c	subroutine find_unit(iunit)
c	iunit = 77
c	end

c	program subfile_main
c	end

c************************************************************






