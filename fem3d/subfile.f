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
     +				,nkn,lmax,nvar,hlv)

c writes header of the file

        implicit none

	logical bformat		!formatted or unformatted
	integer iunit		!file unit
	integer it		!time stamp
	integer nkn		!size of data (horizontal, nodes or elements)
	integer lmax		!vertical values
	integer nvar		!number of variables to write
	real hlv(lmax)		!depth at bottom of layer

	integer nvers		!version of file format
	integer l

	nvers = 1

	if( bformat ) then
	  write(iunit,*) it,nkn,lmax,nvar,nvers
	  write(iunit,*) (hlv(l),l=1,lmax)
	else
	  write(iunit) it,nkn,lmax,nvar,nvers
	  write(iunit) (hlv(l),l=1,lmax)
	end if

	end

c************************************************************

	subroutine write_fem_file_data(bformat,iunit,string
     +				,nkn,nlvdim,ilhkv,data)

c writes data of the file

        implicit none

	logical bformat		!formatted or unformatted
	integer iunit		!file unit
	character*(*) string	!string explaination
	integer nkn		!size of data (horizontal, nodes or elements)
	integer nlvdim		!vertical dimension
	integer ilhkv(nkn)	!number of layers in point k (node)
	real data(nlvdim,nkn)	!data

	integer k,lm,l
	character*80 text

	text = string

	if( bformat ) then
	  write(iunit,*) text
	  do k=1,nkn
	    lm = ilhkv(k)
	    if( nlvdim .eq. 1 ) lm = 1
	    write(iunit,*) k,lm,(data(l,k),l=1,lm)
	  end do
	else
	  write(iunit) text
	  do k=1,nkn
	    lm = ilhkv(k)
	    if( nlvdim .eq. 1 ) lm = 1
	    write(iunit) k,lm,(data(l,k),l=1,lm)
	  end do
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

	integer nvar

	nvar = 1

	call write_fem_file_header(bformat,iunit,it
     +				,nkn,nlv,nvar,hlv)
	call write_fem_file_data(bformat,iunit,string
     +				,nkn,nlvdim,ilhkv,data)

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

	integer nvar,lmax,nlvdim
	real hlv(1)
	integer ilhkv(1)

	nvar = 1
	lmax = 1
	nlvdim = 1
	hlv(1) = 10000.
	ilhkv(1) = 1

	call write_fem_file_header(bformat,iunit,it
     +				,nkn,lmax,nvar,hlv)
	call write_fem_file_data(bformat,iunit,string
     +				,nkn,nlvdim,ilhkv,data)

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

	integer it,nkn0,nkn1,lmax,nvar,nvers

	nkn0 = 0
	nkn1 = 0
	it = 0
	lmax = 0
	nvar = 0
	nvers = 0

c------------------------------------------------------
c first try unformatted
c------------------------------------------------------

	open(iunit,file=name,form='unformatted',status='old',err=2)
	read(iunit,err=1) it,nkn0,lmax,nvar,nvers

	if( nkn .ne. nkn0 ) goto 1
	if( lmax .le. 0 .or. lmax .gt. 1000 ) goto 1
	if( nvar .le. 0 .or. nvar .gt. 100 ) goto 1
	if( nvers .le. 0 .or. nvers .gt. 1 ) goto 1

c       -----------------------------------------------
c	we arrived here... this means the file is unformatted
c       -----------------------------------------------

	backspace(iunit)
	bformat = .false.
	return

    1	continue
	close(iunit)
    2	continue

c------------------------------------------------------
c now try formatted
c------------------------------------------------------

	open(iunit,file=name,form='formatted',status='old',err=9)
	read(iunit,*,err=1) it,nkn1,lmax,nvar,nvers

	if( nkn .ne. nkn1 ) goto 9
	if( lmax .le. 0 .or. lmax .gt. 1000 ) goto 9
	if( nvar .le. 0 .or. nvar .gt. 100 ) goto 9
	if( nvers .le. 0 .or. nvers .gt. 1 ) goto 9

c       -----------------------------------------------
c	we arrived here... this means the file is formatted
c       -----------------------------------------------

	backspace(iunit)
	bformat = .true.
	return

    9	continue

c------------------------------------------------------
c error handling
c------------------------------------------------------

	write(6,*) 'nkn: ',nkn,nkn0,nkn1
	write(6,*) 'params: ',it,lmax,nvar,nvers
	write(6,*) 'cannot open or read file: ',name
	stop 'error stop fem_file_open: file open or read'

c------------------------------------------------------
c end of routine
c------------------------------------------------------

	end

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
	integer l,nkn0

	if( bformat ) then
	  read(iunit,*) it,nkn0,lmax,nvar,nvers
	  read(iunit,*) (hlv(l),l=1,min(lmax,nlvdim))
	else
	  read(iunit) it,nkn0,lmax,nvar,nvers
	  read(iunit) (hlv(l),l=1,min(lmax,nlvdim))
	end if

	if( lmax .gt. nlvdim) goto 98

	if( nkn .ne. nkn0 ) goto 99
	if( lmax .le. 0 .or. lmax .gt. 1000 ) goto 99
	if( nvar .le. 0 .or. nvar .gt. 100 ) goto 99
	if( nvers .le. 0 .or. nvers .gt. 1 ) goto 99

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






