c
c $Id: subsse.f,v 1.25 2009-04-03 16:38:23 georg Exp $
c
c interpolation routines
c
c contents :
c
c subroutine intp_lsqr(nn,ndiff,xx,xee)	least square interpolation
c function intp_lagr(n,x,y,xe)		lagrangian interpolation
c
c function tcomp(nintp,t)		t where new value has to be read
c subroutine intp_ts(iunit,nintp,nvar,t,vars,rint,b3d,bformat)
c		interpolation directly from file (more files, multiple columns)
c
c subroutine read_time_series(unit,nvar,b3d,time,values,ierr)
c		time series read
c subroutine read_0_time_series(unit,nvar,time,values,ierr)
c		classic time series read:  time,val1,val2,val3,...valn
c subroutine read_3_time_series(unit,nvar,time,values,ierr)
c		3D time series read
c
c function exxpp(nintp,nmax,x,y,xe,istart)	interpolation from array
c function rintq(...)			interpolation in a square
c
c exf routines (see later)
c
c revision log :
c
c 29.06.1998	ggu	error check in exxqq9
c 14.07.1998	ggu	new routine exxqq
c 21.08.1998	ggu	xv renamed to xx
c 31.05.2001	ggu	revised exxp (style) and exxpp (error fix and grade)
c 12.02.2003	ggu	return constant if t is out of range in exxqq
c 11.03.2005	ggu	new interpolation also for 3D time series
c 16.05.2005	ggu	more debug and new error messages
c 01.02.2006	ggu	new routine exfpres() for pointer to results
c 17.02.2006	ggu	set iunit to -1 if no file given (exffil)
c 29.02.2008	ggu	name change: exxp -> intp_lagr, extrp -> intp_lsqr
c 29.02.2008	ggu	deleted exxqq, exxqq9
c 17.03.2008	ggu	completely restructured
c 18.04.2008	ggu	new exffild, adjusted exffil, bugfix in exfini, exfintp
c 03.09.2008	ggu	bug fix in intp_ts, bug fix in exfini (aux array)
c 08.10.2008	ggu	introduced matinv0 for independence
c 08.11.2008	ggu	better error handling
c 02.04.2009	ggu	if less data given lower interpolation (REDINT)
c 03.04.2009	ggu	new routine intp_neville() (stable lagrange interpol.)
c
c*************************************************************

c*************************************************************
c*************************************************************
c*************************************************************

	subroutine intp_ts(iunit,nintp,ndata,t,vars,rint,b3d,bformat)

c interpolation directly from file (more files and multiple columns)
c
c interpolation of values read directly from formatted file
c t values must be in increasing order
c rt values for which interpolation has to be performed
c ...must be called for in increasing order
c works also for one point interpolation
c
c iunit		file number from which data is read (if negative reset file)
c nintp		number of points used for interpolation
c ndata		number of data in time step
c t		t value for which y values have to be found
c vars		array with working variables already read
c rint		array of interpolated return values at time t
c b3d		if file format is 3D or 0D
c bformat	if file is formatted (only for 3d)
c
c calls: read_time_series, intp_neville, tcomp

	implicit none

c arguments
	integer iunit
	integer nintp
	integer ndata
	real t
	!real vars(0:ndata,nintp+1)	!error -> only data, not results
	real vars(0:ndata,nintp)
	real rint(ndata)
        logical b3d
        logical bformat
c local
	integer ndim
	parameter (ndim=5)
	integer unit
	integer i,j,ierr
	real tc,time
	character*70 name
	real x(ndim), y(ndim)
c functions
	real intp_neville
	real tcomp
c save
	logical bdebug
	real eps
	save eps,bdebug
	data eps / 1.e-5 /
	data bdebug / .false. /

c----------------------------------------------------------
c some checks
c----------------------------------------------------------

	if( nintp .le. 0 ) goto 90
	if( nintp .gt. ndim ) goto 90

        !b3d = .false.

c----------------------------------------------------------
c rewind if necessary
c----------------------------------------------------------

	unit = iunit

	if( unit .lt. 0 ) then	!reset file
		unit = -unit
		rewind(unit)

		if( bdebug ) then
		  write(6,*) 'intp_ts: (Initializing unit) ',unit
		  write(6,*) unit,ndata,b3d,t
		end if

		do i=1,nintp
		  if(bdebug) write(6,*) 'intp_ts: (reading initial data) ',i
                  call read_time_series(unit,ndata,b3d,time,rint,ierr)
                  if( ierr .gt. 0 ) goto 95
                  !if( ierr .lt. 0 ) goto 94
                  if( ierr .lt. 0 ) goto 2	!REDINT
                  vars(0,i) = time
		  do j=1,ndata
		    vars(j,i) = rint(j)
                  end do
		  if( bdebug ) write(6,*) (vars(j,i),j=0,ndata)
		end do
    2		continue
		if( i .ne. nintp+1 ) then	!reduce intpol if less data
		  write(6,*) '------------------------------'
		  write(6,*) 'file contains not enough data'
		  write(6,*) 'unit,nintp: ',unit,nintp
		  write(6,*) 'at least ',nintp,' data needed'
		  nintp = i - 1
		  write(6,*) 'using lower interpolation: ',nintp
		  write(6,*) '------------------------------'
		end if

		if( nintp .gt. 1 ) then		!check time values
		  do i=2,nintp
		    if( vars(0,i) .le. vars(0,i-1) ) goto 87
		  end do
		end if

		return
	else if( unit .eq. 0 ) then
		goto 93
	end if

c----------------------------------------------------------
c get new t value
c----------------------------------------------------------

	tc = tcomp(ndata,nintp,vars)	!critical t when to read new values

	do while( t .gt. tc )
		!write(6,*) 'reading data for boundary: ',i,t,tc,b3d
                call read_time_series(unit,ndata,b3d,time,rint,ierr)
                if( ierr .gt. 0 ) goto 96
                if( ierr .lt. 0 ) goto 1

		call obc_check_time(unit,nintp,ndata,vars,time)
		call obc_insert_record(nintp,ndata,vars,time,rint)

		tc = tcomp(ndata,nintp,vars)     !pass in time column
	end do
    1	continue

c----------------------------------------------------------
c do the interpolation for every column
c----------------------------------------------------------

	call obc_interpolate(nintp,ndata,vars,t,rint)

c----------------------------------------------------------
c end of routine
c----------------------------------------------------------

	return
   87	continue
	write(6,*) 'time values not ascending: '
	write(6,*) (vars(0,i),i=1,nintp)
	stop 'error stop : intp_ts'
   90	continue
	write(6,*) 'Value for nintp out of range: ',nintp
	write(6,*) 'Possible min/max values: ',1,ndim
	stop 'error stop : intp_ts'
   91	continue
	write(6,*) 'No extrapolation possible'
	write(6,*) 'Actual t value :',t
	write(6,*) 'Available time levels :'
	write(6,*) (vars(0,i),i=1,nintp)
	stop 'error stop : intp_ts'
   93	continue
	write(6,*) 'Cannot read from unit 0'
	stop 'error stop : intp_ts'
   94	continue
	write(6,*) 'End of file while initializing on unit :',unit
	write(6,*) 'i,nintp,ndata: ',i,nintp,ndata
	write(6,*) 'time: ',time
	write(6,*) 'File must consist of at least ',nintp,' data'
	stop 'error stop : intp_ts'
   95	continue
	write(6,*) 'Read error while initializing from file :',unit
	write(6,*) 'i,nintp,ndata: ',i,nintp,ndata
	write(6,*) 'time: ',time
	write(6,*) (rint(i),i=1,ndata)
	write(6,*) 'Attention : Cannot read unformatted file'
	stop 'error stop : intp_ts'
   96	continue
	write(6,*) 'Read error from file :',unit
	write(6,*) 'nintp,ndata: ',nintp,ndata
	write(6,*) 't,tc,time: ',t,tc,time
	write(6,*) (rint(i),i=1,ndata)
	stop 'error stop : intp_ts'
	end

c*************************************************************

	subroutine obc_interpolate(nintp,ndata,vars,time,rint)

	implicit none

	integer nintp		!order of interpolation (4=cubic)
	integer ndata		!total number of data
	real vars(0:ndata,nintp)!values of variables, 0 column is time
	real time		!time for desired interpolated values
	real rint(0:ndata)	!values for time interpolated

	real eps
	parameter (eps=1.e-5)
	integer ndim
	parameter (ndim=10)
	real x(ndim), y(ndim)
	integer i,j

	real intp_neville

	if( nintp .gt. ndim ) stop 'error stop obc_interpolate: ndim'

c	----------------------------------------------------------
c	check if we are really doing an interpolation
c	----------------------------------------------------------

	i = 0
	if( time .lt. vars(0,1)-eps ) i = 1
	if( time .gt. vars(0,nintp)+eps ) i = nintp
	if( nintp .le. 0 ) i = 1

	if( i .gt. 0 ) then	!keep constant -> extrapolation
	  do j=1,ndata
	    rint(j) = vars(j,i)
	  end do
	  return
	end if

c	----------------------------------------------------------
c	check if we are really doing an interpolation
c	----------------------------------------------------------

	do i=1,nintp
	  x(i) = vars(0,i)
	end do

	do j=1,ndata
	  do i=1,nintp
	    y(i) = vars(j,i)
	  end do
	  rint(j) = intp_neville(nintp,x,y,time)
	end do

c	----------------------------------------------------------
c	end of routine
c	----------------------------------------------------------

	end

c*************************************************************

	subroutine obc_check_time(iunit,nintp,ndata,vars,time)

	implicit none

	integer iunit
	integer nintp,ndata
	real vars(0:ndata,nintp)
	real time

	character*60 name

	if( time .gt. vars(0,nintp) ) return

	write(6,*) 'time values not in ascending order'
	write(6,*) 'unit = ',iunit
	write(6,*) 'nintp = ',nintp
	write(6,*) 'time: ',vars(0,nintp),time
	call filna(iunit,name)
	write(6,'(a,a)') 'file = ',name
	stop 'error stop : obc_check_time'
	end

c*************************************************************

	subroutine obc_init_record(nintp,ndata,vars)

	implicit none

	integer nintp,ndata
	real vars(0:ndata,nintp)

	integer i,j

	do i=1,nintp
	  do j=0,ndata
	    vars(j,i) = 0.
          end do
	end do

	end

c*************************************************************

	subroutine obc_insert_record(nintp,ndata,vars,time,rint)

	implicit none

	integer nintp,ndata
	real vars(0:ndata,nintp)
	real time
	real rint(ndata)

	integer i,j

	do i=1,nintp-1
	  do j=0,ndata
	    vars(j,i) = vars(j,i+1)
          end do
	end do

        vars(0,nintp) = time
	do j=1,ndata
	  vars(j,nintp) = rint(j)
        end do

	end

c*************************************************************

	function tcomp(ndata,nintp,t)

c returns value of t where new value has to be read

	implicit none

	real tcomp
	integer ndata			!total size of data
	integer nintp			!grade of interpolation
	real t(0:ndata,nintp+1)		!data, time is at position 0

	integer nold,n1,n2
	save nold,n1,n2
	data nold / 0 /		!impossible value

c----------------------------------------------------------
c if value of nintp has changed -> compute new n1,n2
c
c we could compute this every time, but this slightly more efficient
c----------------------------------------------------------

	if( nintp .ne. nold ) then

	    if( mod(nintp,2) .eq. 0 ) then	!even
		n1=1+nintp/2
		n2=n1
	    else
		n1=1+nintp/2
		if(nintp.gt.1) then
			n2=n1+1
		else
			n2=n1
		end if
	    end if

	    nold = nintp

	end if

c----------------------------------------------------------
c return compare value
c----------------------------------------------------------

	tcomp = 0.5 * ( t(0,n1) + t(0,n2) )

c----------------------------------------------------------
c end of routine
c----------------------------------------------------------

	end

c*************************************************************
c*************************************************************
c*************************************************************

        subroutine read_time_series(unit,ndata,b3d,time,values,ierr)

c time series read

        implicit none

        integer unit
        integer ndata
	logical b3d			!true if 3d read
        real time
        real values(ndata)
        integer ierr

	if( b3d ) then
          call read_3_time_series(unit,ndata,time,values,ierr)
	else
          call read_0_time_series(unit,ndata,time,values,ierr)
	end if

	end

c*************************************************************

        subroutine read_0_time_series(unit,nvar,time,values,ierr)

c classic time series read:  time,val1,val2,val3,...valn

        implicit none

        integer unit
        integer nvar
        real time
        real values(nvar)
        integer ierr

        integer j

	read(unit,*,iostat=ierr) time,(values(j),j=1,nvar)

        end

c*************************************************************

        subroutine read_3_time_series(unit,ndata,time,values,ierr)

c 3D time series read:  
c
c       time,lmax,nk,nvar
c       values

        implicit none

        integer unit
        integer ndata
        real time
        real values(ndata)
        integer ierr

        integer j,n,lmax,nk,nvar
        integer ip,ivar,k,kn
	logical bdebug
	character*80 name

	bdebug = .true.
	bdebug = .false.

	if( bdebug ) write(6,*) '3d TS read : ',unit,ndata

	read(unit,*,iostat=ierr) time,lmax,nk,nvar
	if( bdebug ) write(6,*) time,lmax,nk,nvar
	n = 0
        if( ierr .gt. 0 ) goto 97
        if( ierr .lt. 0 ) return
        n = lmax * nk * nvar
        if( n .ne. ndata ) goto 99

        ip = 0
        do ivar=1,nvar
          do k=1,nk
	    read(unit,*,iostat=ierr) kn,(values(ip+j),j=1,lmax)
	    !write(6,*) k,kn,(values(ip+j),j=1,lmax)
            if( kn .ne. k ) goto 98
            ip = ip + lmax
          end do
        end do

	if( bdebug ) write(6,*) '3d TS read (last ip) : ',unit,ip

        return
   97   continue
	call filna(unit,name)
        write(6,*) 'read error on unit: ',unit,' with file name: '
	write(6,*) name
        write(6,*) time,lmax,nk,nvar
	write(6,*) 'error reading header of data set'
        stop 'error stop read_3_time_series: header'
   98   continue
	call filna(unit,name)
        write(6,*) 'read error on unit: ',unit,' with file name: '
	write(6,*) name
        write(6,*) time,lmax,nk,nvar
        write(6,*) k,kn
        write(6,*) '(node number not compatible)'
        stop 'error stop read_3_time_series: nk'
   99   continue
	call filna(unit,name)
        write(6,*) 'read error on unit: ',unit,' with file name: '
	write(6,*) name
        write(6,*) time,lmax,nk,nvar
        write(6,*) n,ndata
        write(6,*) '(data set may not be compatible)'
        write(6,*) 'calling routine wants ',ndata,' data'
        write(6,*) 'data file provides    ',n,' data'
        write(6,*) 'There could be a mismatch between nbdim and lmax'
        write(6,*) 'In this case please set nbdim in STR to'
        write(6,*) 'lmax given in data file'
        stop 'error stop read_3_time_series: ndata'
        end

c***************************************************************
c***************************************************************
c***************************************************************
c
c exffil	opens file and sets up array
c exffils	short version of exffil
c exffild	short version of exffil with default setting
c exfini	initializes array with file already open
c exfintp	interpolates value
c exfget	get last interpolated values
c exfgetvar	get last interpolated values (only for variable ivar)
c exfset	set values
c exfsetdef	set default values
c exfinfo	info on array
c
c exfunit
c exfnvar
c exfsize
c exfpres
c exfcheck
c
c all information needed is stored in vector array(*)
c
c called routines:
c
c exffil		exfini
c exffils		exffil
c exffild		exffils, 
c exfini		intp_ts
c exfintp		intp_ts
c exfget		-
c exfset		-
c exfinfo		-
c
c intp_ts		read_time_series, intp_neville, tcomp
c read_time_series	read_0_time_series, read_3_time_series
c
c used variable:
c
c       iunit           unit number of file
c       nintp           degree of interpolation (2 linear, 4 cubic)
c       nvar            number of variables
c       nsize           number of data per variable (may be 0 -> 1)
c	ndata		total data per time step (normally nvar)
c       ndim            dimension of array
c       nextra          extra information at beginning of array
c       ires            pointer into array to keep interpolated values
c       nspace          space needed to hold all information in array
c
c       rguard          guard value to check for out-of-bound access
c       
c if nsize = 0 or 1, then normal read, else 3D read
c
c iunit is 0 before initialization
c if file has been opened, iunit > 0, else iunit = -1
c
c formula for computing needed space:
c
c	nspace = 1 + nextra + (nintp+1) * (ndata+1)
c
c       1               guard value at end of array
c       nextra          extra header information
c       ndata+1         to keep variables (ndata) and time
c       nintp+1         to keep nintp time steps and one for result
c
c formula for computing pointer ires (pointer to results):
c
c       ires   = 1 + nextra + nintp * (ndata+1)
c
c filling of array:
c
c       header data				(nextra)
c	time_1,data_1				(1+ndata)
c	time_2,data_2				(1+ndata)
c	...
c	time_nintp,data_nintp			(1+ndata)
c	time_intp,data_intp			(1+ndata)
c       rguard					(1)
c
c header data is:
c
c	iunit,nintp,nvar,nsize,ndata,ndim,nextra,ires,nspace,rguard
c
c Usage: (easy)
c
c call exffils('ps.dat',ndim,array)	!opens and initializes
c ...
c call exfintp(array,t,value)		!interpolates for time t -> value(s)
c
c***************************************************************
c***************************************************************
c***************************************************************

        subroutine exffile(file,nintp,nvar,np,ndim,array)

c opens file and inititializes array
c
c everything needed is in array (unit, vars etc...)

        implicit none

        character*(*) file      !file name
	integer nintp		!grade of interpolation (2=linear,4=cubic)
	integer nvar		!how many vars (columns) to read/interpolate
        integer np              !number of points (horizontal) expected
        integer ndim            !dimension of array
        real array(ndim)        !array with all information

        integer iunit
        integer ifileo

        iunit = 0
	if( file .ne. ' ' ) then
          iunit = ifileo(iunit,file,'form','old')
          if( iunit .le. 0 ) goto 99
	end if

        call exfinit(iunit,nintp,nvar,np,ndim,array)

	return
   99	continue
	write(6,*) 'file = ',file
	stop 'error stop exffil: cannot open file'
        end

c***************************************************************
c***************************************************************
c***************************************************************

        subroutine exffil(file,nintp,nvar,nsize,ndim,array)

c opens file and inititializes array
c
c everything needed is in array (unit, vars etc...)

        implicit none

        character*(*) file      !file name
	integer nintp		!grade of interpolation (2=linear,4=cubic)
	integer nvar		!how many vars (columns) to read/interpolate
        integer nsize           !number of data per variable (may be 0 -> 1)
        integer ndim            !dimension of array
        real array(ndim)        !array with all information

        integer iunit
        integer ifileo

        iunit = 0
	if( file .ne. ' ' ) then
          iunit = ifileo(iunit,file,'form','old')
          if( iunit .le. 0 ) goto 99
	end if

        call exfini(iunit,nintp,nvar,nsize,ndim,array)

	return
   99	continue
	write(6,*) 'file = ',file
	stop 'error stop exffil: cannot open file'
        end

c***************************************************************

        subroutine exffils(file,ndim,array)

c opens file and inititializes array - simplified version

	implicit none

        character*(*) file      !file name
        integer ndim            !dimension of array
        real array(ndim)        !array with all information

        integer nintp           !grade of interpolation (2=linear,4=cubic)
        integer nvar            !how many columns to read/interpolate
        integer np              !number of horizontal points expected

	nintp=2
	nvar=1
        np=0			!unsure about it

	call exffile(file,nintp,nvar,np,ndim,array)

	end

c***************************************************************

	subroutine exffild(file,ndim,array,default)

c opens file and inititializes array with default - simplified version

	implicit none

        character*(*) file      !file name
        integer ndim            !dimension of array
        real array(ndim)        !array with all information
        real default(1)         !default to use in case no file exists

        call exffils(file,ndim,array)
	call exfsetdef(array,default)

	end

c***************************************************************
c***************************************************************
c***************************************************************

	!subroutine obc_init(iunit,nintp,nvar,np,ndim,array)
	subroutine exfinit(iunit,nintp,nvar,np,ndim,array)

c sets up interpolation from file -> all information is in array
c
c       space is computed as follows:
c
c       nvar variables + time -> nvar+1
c       nintp values per variable + interpolated values -> nintp+1
c       nextra is extra information heading the other data
c       one guard value at end of array

	implicit none

	include 'subobc.h'

	integer iunit		!unit of file
	integer nintp		!grade of interpolation (2=linear, 4=cubic)
	integer nvar		!how many columns to read/interpolate
	integer np		!number of horizontal points
	integer ndim		!dimension of array, on return space used
	real array(ndim)	!array with information

	logical b3d,debug,bformat
	integer iformat,lmax,np0,nvar0
	integer ires,nspace,ndata,nnintp,nsize
        integer i
	real time
	character*80 file

	debug = .false.

c	-------------------------------------------------------------
c	set-up parameters
c	-------------------------------------------------------------

        nnintp = nintp
	if( iunit .le. 0 ) nnintp = 0   !reserve some space only for results

	call exfcompsize(iunit,np0,nvar0,lmax,bformat,b3d)

c	if 0d
c		np0 = 1
c		lmax = 1
c		b3d = false

	if( nvar .gt. 0 .and. nvar .ne. nvar0 ) goto 91
	if( np .gt. 0 .and. np0 .gt. 1 .and. np .ne. np0 ) goto 92

	nsize = np0 * lmax
	if( .not. b3d ) nsize = 0
	ndata = nvar0 * max(1,nsize)
	iformat = 0
	if( bformat ) iformat = 1

        ires   = 1 + nextra + nnintp * (ndata+1)	!pointer to results
	nspace = 1 + nextra + (nnintp+1) * (ndata+1)	!total space needed

c	-------------------------------------------------------------
c	check space and exit with error in case
c	-------------------------------------------------------------

	if( nspace .gt. ndim ) then
                write(6,*) '*** error in exfini'
		write(6,*) 'Space in array is not sufficient'
		write(6,*) 'dimension of array: ',ndim
		write(6,*) 'space needed:       ',nspace
		write(6,*) 'formula: 1 + nextra + (nnintp+1) * (ndata+1)'
		write(6,*) 'with nextra = ',nextra
		write(6,*) 'iunit,nnintp,nvar:  ',iunit,nnintp,nvar
		write(6,*) 'nsize,ndata:  ',nsize,ndata
		call filna(iunit,file)
		write(6,*) 'filename: ',file
		stop 'error stop exfini: ndim'
	end if

c	-------------------------------------------------------------
c	initialization of array, using result section as aux array
c	-------------------------------------------------------------

	if(debug) write(6,*) 'exfini: (initializing) ',iunit,nspace

        do i=1,nspace
          array(i) = 0.
        end do

	if( iunit .gt. 0 ) then
	  time = 0.     !is not used
	  call intp_ts(-iunit,nnintp,ndata,time,array(nextra+1)
     +			,array(ires+1),b3d,bformat)
        end if

	if(debug) then
	  write(6,*) 'exfini (end of intp_ts) : ',iunit,nintp,nnintp,nvar
	end if

c	-------------------------------------------------------------
c	write parameters to header of array
c	-------------------------------------------------------------

	array(ip_iunit) = iunit
	array(ip_nintp) = nnintp
	array(ip_nvar) = nvar0
	array(ip_nsize) = nsize
	array(ip_ndata) = ndata
	array(ip_ndim) = ndim
	array(ip_nextra) = nextra
	array(ip_ires) = ires
	array(ip_nspace) = nspace
	array(ip_np) = np0
	array(ip_lmax) = lmax
	array(ip_iformat) = iformat
	array(nextra) = rguard
	array(nspace) = rguard

c	-------------------------------------------------------------
c	in case flag unit as not used
c	-------------------------------------------------------------

	if( iunit .eq. 0 ) array(ip_iunit) = -1	!flag unit as not used

	if(debug) then
	  write(6,*) 'exfini: (finished initializing) ',iunit,nvar
	end if

c	-------------------------------------------------------------
c	end of routine
c	-------------------------------------------------------------

	return
   91	continue
	stop 'error stop exfinit: 91'
   92	continue
	stop 'error stop exfinit: 92'
	end

c***************************************************************

	subroutine exfintp(array,t,rint)

c interpolation from file -> all information is in array
c
c interpolated values are in last part of array

	implicit none

	include 'subobc.h'

	real array(*)		!array with information from set-up
	real t			!t value for which to interpolate
	real rint(1)		!interpolated values

	logical b3d,bformat
	integer iformat
	integer iunit,nintp,nsize,ndata
        integer ires,nspace,i

	iunit   = nint(array(ip_iunit))
	nintp   = nint(array(ip_nintp))
	nsize   = nint(array(ip_nsize))
	ndata   = nint(array(ip_ndata))
	ires    = nint(array(ip_ires))
	iformat = nint(array(ip_iformat))

	!b3d = nsize .gt. 1
	b3d = nsize .gt. 0
	bformat = iformat .gt. 0

	if( iunit .gt. 0 ) then
	  call intp_ts(iunit,nintp,ndata,t,array(nextra+1)
     +				,rint,b3d,bformat)
          array(ires) = t
          do i=1,ndata
            array(ires+i) = rint(i)
          end do
	else			!default data is already in place
          do i=1,ndata
            rint(i) = array(ires+i)
          end do
	end if

	call exfcheck(array)

	end

c***************************************************************

	subroutine exfget(array,t,rint)

c get last interpolated values (ndata values)

	implicit none

	include 'subobc.h'

	real array(*)		!array with information from set-up
	real t			!t value for which to interpolate
	real rint(1)		!interpolated values

        integer ndata,ires,i

	ndata  = nint(array(ip_ndata))
	ires   = nint(array(ip_ires))

        t = array(ires)
        do i=1,ndata
          rint(i) = array(ires+i)
        end do

	call exfcheck(array)

	end

c***************************************************************

	subroutine exfgetvar(array,ivar,t,rint)

c get last interpolated values (only for variable ivar, nsize values)

	implicit none

	include 'subobc.h'

	real array(*)		!array with information from set-up
	integer ivar		!number of variable needed
	real t			!t value for which to interpolate
	real rint(1)		!interpolated values

        integer nvar,nsize
	integer ires,ip,i

	nvar   = nint(array(ip_nvar))
	nsize  = max(nint(array(ip_nsize)),1)
	ires   = nint(array(ip_ires))

	if( ivar .gt. nvar ) then
	  write(6,*) 'ivar = ',ivar,'   nvar = ',nvar
	  stop 'error stop exfgetvar: ivar too high'
	end if

        t = array(ires)
	ip = ires + (ivar-1)*nsize	!pointer to results of variable ivar
        do i=1,nsize
          rint(i) = array(ip+i)
        end do

	call exfcheck(array)

	end

c***************************************************************

	subroutine exfset(array,t,rint)

c sets new actual values

	implicit none

	include 'subobc.h'

	real array(*)		!array with information from set-up
	real t			!t value for which to interpolate
	real rint(1)		!interpolated values

        integer ndata,ires,i

	ndata  = nint(array(ip_ndata))
	ires   = nint(array(ip_ires))

        array(ires) = t
        do i=1,ndata
          array(ires+i) = rint(i)
        end do

	call exfcheck(array)

	end

c***************************************************************

	subroutine exfsetdef(array,rdef)

c sets default value (one for each variable)

	implicit none

	include 'subobc.h'

	real array(*)		!array with information from set-up
	real rdef(*)		!default values for every variable

        integer nsize,nvar
	integer ires,i
	integer ivar,ip

	nvar   = nint(array(ip_nvar))
	nsize  = max(nint(array(ip_nsize)),1)
	ires   = nint(array(ip_ires))

        array(ires) = 0.	!time - not important
	do ivar=1,nvar
	  ip = ires + (ivar-1)*nsize	!pointer to results of variable ivar
          do i=1,nsize
            array(ip+i) = rdef(ivar)
	  end do
        end do

	call exfcheck(array)

	end

c***************************************************************

	subroutine exfinfo(ipunit,array)

c interpolation from file -> info

	implicit none

	include 'subobc.h'

	integer ipunit		!unit where to print on (<0 -> 6)
	real array(*)		!array with information from set-up

	logical bdebug
	integer iunit,nintp,nvar,nsize,ndata,ndim
        integer ires,nspace
	integer iformat,lmax,np
	integer ipu
	integer ip,in,i

	bdebug = .false.
	bdebug = .true.

	ipu = ipunit
	if( ipunit .le. 0 ) ipu = 6

	iunit   = nint(array(ip_iunit))
	nintp   = nint(array(ip_nintp))
	nvar    = nint(array(ip_nvar))
	nsize   = nint(array(ip_nsize))
	ndata   = nint(array(ip_ndata))
	ndim    = nint(array(ip_ndim))
	ires    = nint(array(ip_ires))
	nspace  = nint(array(ip_nspace))
	np      = nint(array(ip_np))
	lmax    = nint(array(ip_lmax))
	iformat = nint(array(ip_iformat))

        write(ipu,*) 'info on array interpolation:'
        write(ipu,*) 'unit     : ',iunit
        write(ipu,*) 'nintp    : ',nintp
        write(ipu,*) 'nvar     : ',nvar
        write(ipu,*) 'nsize    : ',nsize
        write(ipu,*) 'ndata    : ',ndata
        write(ipu,*) 'ndim     : ',ndim
        write(ipu,*) 'nextra   : ',nextra
        write(ipu,*) 'ires     : ',ires
        write(ipu,*) 'nspace   : ',nspace
        write(ipu,*) 'np       : ',np
        write(ipu,*) 'lmax     : ',lmax
        write(ipu,*) 'iformat  : ',iformat
        write(ipu,*) 'rguard   : ',array(nspace)

	if( bdebug ) then
	  ip = 1 + nextra
	  do in=1,nintp
	    write(ipu,*) 'intp: level,pointer,time: ',in,ip,array(ip)
	    write(ipu,*) (array(ip+i),i=1,ndata)
	    ip = ip + ndata + 1
	  end do
	  write(ipu,*) 'result: level,pointer,time: ',0,ip,array(ip)
	  write(ipu,*) (array(ires+i),i=1,ndata)
	end if

        end

c***************************************************************

	subroutine exfunit(array,iunit)

c returns information on unit number

	implicit none

	include 'subobc.h'

	real array(*)		!array with information from set-up
        integer iunit           !unit number of file, 0 if not initialized

	iunit = nint(array(ip_iunit))

        end

c***************************************************************

	subroutine exfnvar(array,nvar)

c returns information on number of variables

	implicit none

	include 'subobc.h'

	real array(*)		!array with information from set-up
        integer nvar            !total number of variables

	nvar = nint(array(ip_nvar))

	end

c***************************************************************

	subroutine exfsize(array,nvar,nsize,ndata)

c returns information on size of data

	implicit none

	include 'subobc.h'

	real array(*)		!array with information from set-up
        integer nvar            !total number of variables
        integer nsize           !data per variable
        integer ndata           !total data in array

	nvar  = nint(array(ip_nvar))
	nsize = nint(array(ip_nsize))
	ndata = nint(array(ip_ndata))

        end

c***************************************************************

	subroutine exfpres(array,ires)

c returns pointer to interpolated result
c
c array(ires)	actual time of interpolated results
c array(ires+i)	actual interpolated result of variable i

	implicit none

	include 'subobc.h'

	real array(*)		!array with information from set-up
        integer ires            !pointer to results in array

	ires   = nint(array(ip_ires))

        end

c***************************************************************

	subroutine exfcheck(array)

c checks array for guard values

	real array(*)

	include 'subobc.h'

	integer nspace

	nspace = nint(array(ip_nspace))

        if( array(nextra) .ne. rguard ) then
            stop 'error stop exfintp: first guard value altered'
     	else if( array(nspace) .ne. rguard ) then
            stop 'error stop exfintp: last guard value altered'
        end if

	end

c***************************************************************

	subroutine exfcompsize(iunit,np,nvar,lmax,bformat)

c gets info on file

	implicit none

	integer iunit		!unit number (input)
	integer np		!number of points in file (0 for time series)
	integer nvar		!number of variables (normally 1)
	integer lmax		!number of levels
	logical bformat		!is file formatted?

	integer it,nvers,ntype,ierr
	real f(10)

	integer iret
	character*1000 line
	integer iscanf

	np = 0
	nvar = 1
	lmax = 1
	bformat = .true.

	if( iunit .le. 0 ) return

c try 3D (fem-file) format

        call fem_file_get_params(bformat,iunit,it
     +                          ,nvers,np,lmax,nvar,ntype,ierr)
	if( ierr .eq. 0 ) return

c try 0D (time series) format

	np = 0
	read(iunit,'(a)',end=1,err=1) line
	iret = iscanf(line,f,0)
	if( iret .lt. 0 ) iret = -iret - 1
	nvar = iret
	backspace(iunit)
	return

    1	continue
	write(6,*) 'Cannot read file on unit: ',iunit
	stop 'error stop exfcompsize: read error'
	end

c***************************************************************
c***************************************************************
c***************************************************************
c***************************************************************
c***************************************************************

	subroutine exfsetdebug(debug)
	logical debug
	logical bdebug
	common /exfdebug/bdebug
	save /exfdebug/
	bdebug=debug
	end

	subroutine exfwritedebug(n,a)
	integer n
	real a(n)
	logical bdebug
	common /exfdebug/bdebug
	save /exfdebug/
	if( bdebug ) then
	  write(91,*) '------exfwritedebug------'
	  write(91,*) n,(a(i),i=1,n)
	  write(91,*) '-------------------------'
	end if
	end
	
c***************************************************************

