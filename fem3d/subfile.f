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
c
c************************************************************

	subroutine write_fem_file_header(bformat,iunit,it
     +				,nkn,lmax,nvar,nvers,ntype,hlv)

c writes header of the file

        implicit none

	logical bformat		!formatted or unformatted
	integer iunit		!file unit
	integer it		!time stamp
	integer nkn		!size of data (horizontal, nodes or elements)
	integer lmax		!vertical values
	integer nvar		!number of variables to write
	integer nvers		!version of file format
	integer ntype		!type of information contained
	real hlv(lmax)		!depth at bottom of layer

	integer l,nv

	nv = nvers
	if( nv .eq. 0 ) nv = 2	!default
	if( nv .lt. 1 .and. nv .gt. 2 ) goto 99

	if( bformat ) then
	  write(iunit,1000) it,nkn,lmax,nvar,nv,ntype
	  write(iunit,*) (hlv(l),l=1,lmax)
	else
	  write(iunit) it,nkn,lmax,nvar,nv,ntype
	  write(iunit) (hlv(l),l=1,lmax)
	end if

	return
 1000	format(2i14,4i8)
   99	continue
	write(6,*) 'nvers = ',nvers
	stop 'error stop write_fem_file_header: nvers'
	end

c************************************************************

	subroutine write_fem_file_data(bformat,iunit
     +				,nkn,lmax,nvers
     +				,nlvdim,ilhkv
     +				,string,data)

c writes data of the file

        implicit none

	logical bformat		!formatted or unformatted
	integer iunit		!file unit
	integer nkn		!size of data (horizontal, nodes or elements)
	integer lmax		!vertical values (1 for 2d)
	integer nvers		!version of file format
	integer nlvdim		!vertical dimension of data
	integer ilhkv(nkn)	!number of layers in point k (node)
	character*(*) string	!string explaination
	real data(nlvdim,nkn)	!data

	logical bcompact,b2d
	integer k,lm,l,nv
	character*80 text

	nv = nvers
	if( nv .eq. 0 ) nv = 2	!default

	text = string
	bcompact = nv .eq. 2
	b2d = lmax .eq. 1 .or. nlvdim .eq. 1

	if( bformat ) then
	  write(iunit,*) text
	  if( bcompact .and. b2d ) then
	    write(iunit,*) (data(1,k),k=1,nkn)
	  else
	    do k=1,nkn
	      lm = ilhkv(k)
	      if( b2d ) lm = 1
	      write(iunit,*) k,lm,(data(l,k),l=1,lm)
	    end do
	  end if
	else
	  write(iunit) text
	  if( bcompact .and. b2d ) then
	    write(iunit) (data(1,k),k=1,nkn)
	  else
	    do k=1,nkn
	      lm = ilhkv(k)
	      if( b2d ) lm = 1
	      write(iunit) k,lm,(data(l,k),l=1,lm)
	    end do
	  end if
	end if

	end

c************************************************************

	subroutine write_fem_file_3d(bformat,iunit,it
     +				,nkn,nlv,nlvdim,hlv
     +				,ilhkv,string,data)

c writes 1 variable of a 3d field

	implicit none

	logical bformat		!formatted or unformatted
	integer iunit		!file unit
	integer it		!time stamp
	integer nkn		!size of data (horizontal, nodes or elements)
	integer nlv		!vertical values
	integer nlvdim		!vertical dimension
	real hlv(nlv)		!depth at bottom of layer
	integer ilhkv(nkn)	!number of layers in point k (node)
	character*(*) string	!string explaination
	real data(nlvdim,nkn)	!data

	integer nvers,ntype
	integer nvar

	nvers = 2
	ntype = 0
	nvar = 1

	call write_fem_file_header(bformat,iunit,it
     +				,nkn,nlv,nvar,nvers,ntype,hlv)
	call write_fem_file_data(bformat,iunit
     +				,nkn,nlv,nvers
     +				,nlvdim,ilhkv,string,data)

	end

c************************************************************

	subroutine write_fem_file_2d(bformat,iunit,it
     +				,nkn,string,data)

c writes 1 variable of a 2d field

	implicit none

	logical bformat		!formatted or unformatted
	integer iunit		!file unit
	integer it		!time stamp
	integer nkn		!size of data (horizontal, nodes or elements)
	character*(*) string	!string explaination
	real data(nkn)		!data

	integer nvers,ntype
	integer nvar,nlv,nlvdim
	real hlv(1)
	integer ilhkv(1)

	nvers = 2
	ntype = 0
	nvar = 1
	nlv = 1
	nlvdim = 1
	hlv(1) = 10000.
	ilhkv(1) = 1

	call write_fem_file_header(bformat,iunit,it
     +				,nkn,nlv,nvar,nvers,ntype,hlv)
	call write_fem_file_data(bformat,iunit
     +				,nkn,nlv,nvers
     +				,nlvdim,ilhkv,string,data)

	end

c************************************************************
c************************************************************
c************************************************************

	subroutine fem_file_read_open(name,nkn,iunit,bformat)

c checks if file is formatted or unformatted and opens it for read

	implicit none

	character*(*) name	!string explaination
	integer nkn		!expected size of data
	integer iunit		!unit of opened file (in/out)
	logical bformat		!is formatted? (return)

	logical bok

	call fem_file_test_open(name,nkn,bok,bformat)

	if( bok ) then
	  if( bformat ) then
	    open(iunit,file=name,form='formatted',status='old')
	  else
	    open(iunit,file=name,form='unformatted',status='old')
	  end if
	else
	  call fem_file_write_info(name,.true.)
	  call fem_file_write_info(name,.false.)
	  stop 'error stop fem_file_read_open: file open or read'
	end if

	end

c************************************************************

	subroutine fem_file_write_info(name,bformat)

c writes information on file from header

	character*(*) name	!string explaination
	logical bformat		!is formatted?

	integer iunit
	integer it,nkn,lmax,nvar,nvers,ntype

	it = 0
	nkn = 0
	lmax = 0
	nvar = 0
	nvers = 0
	ntype = 0

	iunit = 90
	call find_unit(iunit)

	if( bformat ) then
	  open(iunit,file=name,form='formatted',status='old')
	  read(iunit,*,err=1) it,nkn,lmax,nvar,nvers,ntype
	  write(6,*) 'formatted read: '
	else
	  open(iunit,file=name,form='unformatted',status='old')
	  read(iunit,err=1) it,nkn,lmax,nvar,nvers,ntype
	  write(6,*) 'unformatted read: '
	end if
    1	continue

	write(6,*) it,nkn,lmax,nvar,nvers,ntype

	close(iunit)

	end

c************************************************************

	subroutine fem_file_test_open(name,nkn,bok,bformat)

c checks if file is formatted or unformatted and opens it for read

	implicit none

	character*(*) name	!string explaination
	integer nkn		!expected size of data, 0 if no idea
	logical bok		!successful read? (return)
	logical bformat		!is formatted? (return)

	integer iunit
	integer it,nkn0,nkn1,lmax,nvar,nvers,ntype

c------------------------------------------------------
c initialize parameters
c------------------------------------------------------

	nkn0 = 0
	nkn1 = 0
	it = 0
	lmax = 0
	nvar = 0
	nvers = 0

	bok = .true.

c------------------------------------------------------
c find unit to open file
c------------------------------------------------------

	iunit = 90
	call find_unit(iunit)

c------------------------------------------------------
c first try unformatted
c------------------------------------------------------

	open(iunit,file=name,form='unformatted',status='old',err=2)
	read(iunit,err=1) it,nkn0,lmax,nvar,nvers,ntype

	if( nkn .gt. 0 .and. nkn .ne. nkn0 ) goto 1
	if( lmax .le. 0 .or. lmax .gt. 1000 ) goto 1
	if( nvar .le. 0 .or. nvar .gt. 100 ) goto 1
	if( nvers .lt. 1 .or. nvers .gt. 2 ) goto 1
	if( ntype .lt. 0 .or. ntype .gt. 2 ) goto 1

c       -----------------------------------------------
c	we arrived here... this means the file is (probably) unformatted
c       -----------------------------------------------

	close(iunit)
	bformat = .false.
	return

    1	continue
	close(iunit)
    2	continue

c------------------------------------------------------
c now try formatted
c------------------------------------------------------

	open(iunit,file=name,form='formatted',status='old',err=8)
	read(iunit,*,err=1) it,nkn1,lmax,nvar,nvers

	if( nkn .gt. 0 .and. nkn .ne. nkn0 ) goto 9
	if( lmax .le. 0 .or. lmax .gt. 1000 ) goto 9
	if( nvar .le. 0 .or. nvar .gt. 100 ) goto 9
	if( nvers .le. 0 .or. nvers .gt. 1 ) goto 9
	if( ntype .lt. 0 .or. ntype .gt. 2 ) goto 9

c       -----------------------------------------------
c	we arrived here... this means the file is (probably) formatted
c       -----------------------------------------------

	close(iunit)
	bformat = .true.
	return

    9	continue
	close(iunit)
    8	continue

c------------------------------------------------------
c no successful opening
c------------------------------------------------------

	bok = .false.

c------------------------------------------------------
c end of routine
c------------------------------------------------------

	end

c************************************************************
c************************************************************
c************************************************************

	subroutine read_fem_file_header(bformat,iunit,it
     +				,nkn,lmax,nlvdim,nvar,hlv)

c reads header of the file

        implicit none

	logical bformat		!formatted or unformatted
	integer iunit		!file unit
	integer it		!time stamp
	integer nkn		!size of data (horizontal, nodes or elements)
	integer lmax		!vertical values
	integer nlvdim		!vertical dimension
	integer nvar		!number of variables to write
	real hlv(nlvdim)	!depth at bottom of layer

	integer nvers		!version of file format
	integer ntype
	integer l,nkn0

	if( bformat ) then
	  read(iunit,*) it,nkn0,lmax,nvar,nvers,ntype
	  read(iunit,*) (hlv(l),l=1,min(lmax,nlvdim))
	else
	  read(iunit) it,nkn0,lmax,nvar,nvers,ntype
	  read(iunit) (hlv(l),l=1,min(lmax,nlvdim))
	end if

	if( lmax .gt. nlvdim) goto 98

	if( nkn .gt. 0 .and. nkn .ne. nkn0 ) goto 99
	if( lmax .le. 0 .or. lmax .gt. 1000 ) goto 99
	if( nvar .le. 0 .or. nvar .gt. 100 ) goto 99
	if( nvers .lt. 1 .or. nvers .gt. 2 ) goto 99

	return
   98	continue
	write(6,*) 'nlvdim,lmax: ',nlvdim,lmax
	stop 'error stop read_fem_file_header: dimension nlvdim'
   99	continue
	write(6,*) it,nkn0,lmax,nvar,nvers
	write(6,*) nkn,nkn0
	stop 'error stop read_fem_file_header: parameters'
	end

c************************************************************

	subroutine read_fem_file_data(bformat,iunit,string
     +				,nkn,nlv,nlvdim,hlv,ilhkv
     +				,lmax,hlvdata,data)

c reads data of the file

        implicit none

	logical bformat		!formatted or unformatted
	integer iunit		!file unit
	character*(*) string	!string explaination
	integer nkn		!size of data (horizontal, nodes or elements)
	integer nlv		!vertical FE levels
	integer nlvdim		!vertical dimension
	real hlv(nlvdim)	!depth of layers in FEM
	integer ilhkv(nkn)	!number of layers in point k (node)
	integer lmax		!maximum layers in data
	real hlvdata(lmax)	!depth of layers in data
	real data(nlvdim,nkn)	!data

	integer ndim
	parameter (ndim=1000)

	logical bcons
	integer k,i,lm,l,lfem
	character*80 text
	real value(ndim+1)
	real hfem(0:nlvdim)
	real hdata(0:ndim+1)

	if( lmax .gt. ndim ) goto 97

	bcons = .false.		!conserve total quantity?

	hfem(0) = 0.
	do l=1,nlv
	  hfem(l) = hlv(l)
	end do

	hdata(0) = 0.
	do l=1,lmax
	  hdata(l) = hlvdata(l)
	end do

	if( bformat ) then
	  read(iunit,*) text
	  do k=1,nkn
	    read(iunit,*) i,lm,(value(l),l=1,min(lm,lmax))
	    if( lm .gt. lmax ) goto 99
	    if( i .ne. k ) goto 98
	    lfem = ilhkv(k)
	    if( nlvdim .eq. 1 ) lfem = 1
	    call intp_vert(bcons,lm,hdata,value,lfem,hfem,data(1,k))
	  end do
	else
	  read(iunit) text
	  do k=1,nkn
	    read(iunit) i,lm,(value(l),l=1,min(lm,lmax))
	    if( lm .gt. lmax ) goto 99
	    if( i .ne. k ) goto 98
	    lfem = ilhkv(k)
	    if( nlvdim .eq. 1 ) lfem = 1
	    call intp_vert(bcons,lm,hdata,value,lfem,hfem,data(1,k))
	  end do
	end if

	string = text

	return
   97	continue
	write(6,*) 'ndim,lmax: ',ndim,lmax
	stop 'error stop read_fem_file_data: dimension ndim'
   98	continue
	write(6,*) 'i,k: ',i,k
	stop 'error stop read_fem_file_data: index mismatch'
   99	continue
	write(6,*) 'lm,ndim: ',lm,ndim
	stop 'error stop read_fem_file_data: dimension lmax'
	end

c************************************************************

	subroutine read_fem_file_3d(bformat,iunit,it
     +				,nkn,nlv,nlvdim,hlv
     +				,ilhkv,string,data)

c reads 1 variable of a 3d field

	implicit none

	logical bformat		!formatted or unformatted
	integer iunit		!file unit
	integer it		!time stamp
	integer nkn		!size of data (horizontal, nodes or elements)
	integer nlv		!vertical values
	integer nlvdim		!vertical dimension
	real hlv(nlv)		!depth at bottom of layer
	integer ilhkv(nkn)	!number of layers in point k (node)
	character*(*) string	!string explaination
	real data(nlvdim,nkn)	!data

	integer ndim
	parameter(ndim=1000)
	real hlvdata(ndim)

	integer nvar,lmax

	call read_fem_file_header(bformat,iunit,it
     +				,nkn,lmax,ndim,nvar,hlvdata)

	if( nvar .ne. 1 ) goto 99

	call read_fem_file_data(bformat,iunit,string
     +				,nkn,nlv,nlvdim,hlv,ilhkv
     +				,lmax,hlvdata,data)

	return
   99	continue
	write(6,*) 'nvar: ',nvar
	stop 'error stop read_fem_file_data: can read only one variable'
	end

c************************************************************

	subroutine read_fem_file_2d(bformat,iunit,it
     +				,nkn,string,data)

c reads 1 variable of a 2d field

	implicit none

	logical bformat		!formatted or unformatted
	integer iunit		!file unit
	integer it		!time stamp
	integer nkn		!size of data (horizontal, nodes or elements)
	character*(*) string	!string explaination
	real data(nkn)		!data

	integer nvar,lmax,nlvdim,nlv
	real hlv(1)
	real hlvdata(1)
	integer ilhkv(1)

	nlv = 1
	nlvdim = 1
	hlv(1) = 10000.
	ilhkv(1) = 1

	call read_fem_file_header(bformat,iunit,it
     +				,nkn,lmax,nlvdim,nvar,hlvdata)

	if( nvar .ne. 1 ) goto 99
	if( lmax .ne. 1 ) goto 98

	call read_fem_file_data(bformat,iunit,string
     +				,nkn,nlv,nlvdim,hlv,ilhkv
     +				,lmax,hlvdata,data)

	return
   98	continue
	write(6,*) 'lmax: ',lmax
	stop 'error stop read_fem_file_2d: lmax'
   99	continue
	write(6,*) 'nvar: ',nvar
	stop 'error stop read_fem_file_2d: can read only one variable'
	end

c************************************************************






