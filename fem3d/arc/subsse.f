c
c $Id: subsse.f,v 1.16 2007-09-27 14:22:53 georg Exp $
c
c interpolation routines
c
c contents :
c
c subroutine extrp(nn,ndiff,xx,xee)	least square interpolation
c function exxp(n,x,y,xe)		regular interpolation
c function tcomp(nintp,t)		t where new value has to be read
c subroutine exxqq(iunit,nintp,nvar,t,vars,rint)
c	interpolation directly from file (more files and multiple columns)
c function exxqq9(iunit,nintp,rx,x,y)	interpolation directly from file (more)
c function exxpp(nintp,nmax,x,y,xe,istart)	interpolation from array
c function rintq(...)			interpolation in a square
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
c
c**************************************************
c
	subroutine extrp(nn,ndiff,xx,xee)
c
c interpolation by least square method
c
c nn		number of points used (4=cubic)
c ndiff		number of grades the interpolating polynome
c		...is lowered (normally 1)
c xx		on entry  : x-values used for interpolation
c		on return : c-values, that have to be used
c		...in computing interpolated y-value -->
c		...ye=c(i)*y(i) ; i=1,n
c xee		x-value for which y-value has to be computed
c
	parameter(iexpmx=5,idim=iexpmx+1)
	real x(idim,idim),a(idim,idim)
	real b(idim,idim)
	real xx(idim),xe(idim)
c
	if(nn.le.ndiff) goto 99
	if(ndiff.le.0) goto 98
	if(nn.gt.iexpmx+1) goto 97
c
	n=nn-ndiff
c
	do ii=1,nn
	x(ii,1)=1.
	do i=2,n
	x(ii,i)=x(ii,i-1)*xx(ii)
	end do
	end do
c
	do i=1,n
	do j=1,n
	h=0.
	do ii=1,nn
	h=h+x(ii,j)*x(ii,i)
	end do
	a(i,j)=h
	end do
	end do
c
	call matinv(a,xx,xe,n,idim)	!xx,xe aux vectors
c
	do i=1,n
	do ii=1,nn
	h=0.
	do j=1,n
	h=h+a(i,j)*x(ii,j)
	end do
	b(i,ii)=h
	end do
	end do
c
	xe(1)=1.
	do i=2,n
	xe(i)=xe(i-1)*xee
	end do
c
	do ii=1,nn
	h=0.
	do i=1,n
	h=h+b(i,ii)*xe(i)
	end do
	xx(ii)=h
	end do
c
	return
c
   97	continue
	write(6,*) 'dimension error'
	write(6,*) 'iexpmx,nn :',iexpmx,nn
	stop 'error stop : extrp'
c
   98	continue
	write(6,*) 'ndiff must be greater than 0'
	write(6,*) 'ndiff =',ndiff
	stop 'error stop : extrp'
c
   99	continue
	write(6,*) 'nn must be greater than ndiff'
	write(6,*) 'nn,ndiff =',nn,ndiff
	stop 'error stop : extrp'
c
	end
c
c**********************************************

	function exxp(nintp,x,y,xe)

c interpolation routine (Lagrangian interpolation)
c
c exxp		interpolated y-value for xe
c nintp		number of points used for interpolation (4=cubic)
c x,y		x/y-values (vector) ( x(i) != x(j) for i != j )
c xe		x-value for which y-value has to computed

	implicit none

	real exxp
	integer nintp
	real x(nintp),y(nintp)
	real xe

	integer i,ii
	real f,g

	f=0.
	do i=1,nintp
	  g=1.
	  do ii=1,nintp
	    if(i.ne.ii) then
		g=g*(xe-x(ii))/(x(i)-x(ii))
	    end if
	  end do
	  f=f+y(i)*g
	end do

	exxp=f

	end

c*************************************************************

	function tcomp(nintp,t)

c returns value of t where new value has to be read

	implicit none

	real tcomp
	integer nintp
	real t(1)

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

	tcomp = 0.5 * ( t(n1) + t(n2) )

c----------------------------------------------------------
c end of routine
c----------------------------------------------------------

	end

c*************************************************************

        subroutine exxqq(iunit,nintp,nvar,t,vars,rint)
        call intp_ts(iunit,nintp,nvar,t,vars,rint,.false.)
        end

        subroutine intp_0_ts(iunit,nintp,nvar,t,vars,rint)
        call intp_ts(iunit,nintp,nvar,t,vars,rint,.false.)
        end

        subroutine intp_3_ts(iunit,nintp,nvar,t,vars,rint)
        call intp_ts(iunit,nintp,nvar,t,vars,rint,.true.)
        end

c*************************************************************

	subroutine intp_ts(iunit,nintp,nvar,t,vars,rint,b3d)

c interpolation directly from file (more files and multiple columns)
c
c interpolation of values read directly from formatted file
c records of file must be written as : write(iunit,*) x y_1 y_2 ... y_nvar
c x values must be in increasing order
c rx values for which interpolation has to be performed
c ...must be called for in increasing order
c works also for one point interpolation
c
c iunit		file number from which data is read (if negative reset file)
c nintp		number of points used for interpolation
c nvar		number of columns to interpolate (normal time series: 1)
c t		t value for which y values have to be found
c vars		array with working variables already read
c rint		array of interpolated return values at time t
c b3d		if file format is 3D or 0D

	implicit none

c arguments
	integer iunit
	integer nintp,nvar
	real t
	real vars(nintp,0:nvar)
	real rint(nvar)
        logical b3d
c local
	integer unit
	integer i,j,ierr
	real tc,time
	character*70 name
c functions
	real exxp
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
		  write(6,*) unit,nvar,b3d,t
		end if

		do i=1,nintp
		  if(bdebug) write(6,*) 'intp_ts: (reading initial data) ',i
                  call read_time_series(unit,nvar,b3d,time,rint,ierr)
                  if( ierr .gt. 0 ) goto 95
                  if( ierr .lt. 0 ) goto 94
		  !read(unit,*,end=94,err=95) time,(rint(j),j=1,nvar)
                  vars(i,0) = time
		  do j=1,nvar
		    vars(i,j) = rint(j)
                  end do
		  if( bdebug ) write(6,*) (vars(i,j),j=0,nvar)
		end do

		return
	else if( unit .eq. 0 ) then
		goto 93
	end if

c----------------------------------------------------------
c get new t value
c----------------------------------------------------------

	tc = tcomp(nintp,vars(1,0))	!critical t when to read new values

	do while( t .gt. tc )
		!write(6,*) 'reading data for boundary: ',i,t,tc
                call read_time_series(unit,nvar,b3d,time,rint,ierr)
                if( ierr .gt. 0 ) goto 96
                if( ierr .lt. 0 ) goto 1
		!read(iunit,*,end=1,err=96) time,(rint(j),j=1,nvar)

		if( bdebug ) then
                        write(6,*) 'new... ',time,(rint(j),j=1,nvar)
                end if

		do j=0,nvar
		  do i=1,nintp-1
		    vars(i,j) = vars(i+1,j)
		  end do
		end do

                vars(nintp,0) = time
		do j=1,nvar
		  vars(nintp,j) = rint(j)
                end do

		tc = tcomp(nintp,vars(1,0))     !pass in time column
	end do
    1	continue

c----------------------------------------------------------
c debug output
c----------------------------------------------------------

	if( bdebug ) then
	  write(6,*) 'debug for intp_ts: '
	  do i=1,nintp
	    write(6,*) (vars(i,j),j=0,nvar)
	  end do
	end if

c----------------------------------------------------------
c time series must have t value monotonically increasing
c----------------------------------------------------------

	do i=2,nintp
	    if( vars(i,0) .le. vars(i-1,0) ) goto 88
	end do

c----------------------------------------------------------
c check if we are really doing an interpolation
c----------------------------------------------------------

	if( nintp .gt. 1 ) then			!no check for nintp = 1
		!if( t .lt. vars(1,0)-eps ) goto 91
		!if( t .gt. vars(nintp,0)+eps ) goto 91
		i = 0
		if( t .lt. vars(1,0)-eps ) i = 1
		if( t .gt. vars(nintp,0)+eps ) i = nintp
		if( i .gt. 0 ) then	!keep constant
		  do j=1,nvar
		    rint(j) = vars(i,j)
		  end do
		  return
		end if
	end if

c----------------------------------------------------------
c do the interpolation for every column
c----------------------------------------------------------

c	real vars(nintp,0:nvar)

	do j=1,nvar
	  rint(j) = exxp(nintp,vars(1,0),vars(1,j),t)
	  !write(87,*) j,t,rint(j)
	  !write(87,*) (vars(i,0),i=1,nintp)
	  !write(87,*) (vars(i,j),i=1,nintp)
	end do

c----------------------------------------------------------
c end of routine
c----------------------------------------------------------

	return
   88	continue
	write(6,*) 't values are not in increasing order:'
	do i=1,nintp
	    write(6,*) (vars(i,j),j=0,nvar)
	end do
	write(6,*) 'interpolation grade = ',nintp
	write(6,*) 't = ',t
	write(6,*) 'unit = ',unit
	call filna(iunit,name)
	write(6,'(a,a)') 'file = ',name
	stop 'error stop : exxqq'
   90	continue
	write(6,*) 'Value for nintp not allowed : ',nintp
	stop 'error stop : exxqq'
   91	continue
	write(6,*) 'No extrapolation possible'
	write(6,*) 'Actual t value :',t
	write(6,*) 'Available time levels :'
	write(6,*) (vars(i,0),i=1,nintp)
	stop 'error stop : exxqq'
   93	continue
	write(6,*) 'Cannot read from unit 0'
	stop 'error stop : exxqq'
   94	continue
	write(6,*) 'End of file while initializing on unit :',unit
	write(6,*) 'i,nintp,nvar: ',i,nintp,nvar
	write(6,*) 'time: ',time
	write(6,*) 'File must consist of at least ',nintp,' data'
	stop 'error stop : exxqq'
   95	continue
	write(6,*) 'Read error while initializing from file :',unit
	write(6,*) 'i,nintp,nvar: ',i,nintp,nvar
	write(6,*) 'time: ',time
	write(6,*) (rint(i),i=1,nvar)
	write(6,*) 'Attention : Cannot read unformatted file'
	stop 'error stop : exxqq'
   96	continue
	write(6,*) 'Read error from file :',unit
	write(6,*) 'nintp,nvar: ',nintp,nvar
	write(6,*) 't,tc,time: ',t,tc,time
	write(6,*) (rint(i),i=1,nvar)
	stop 'error stop : exxqq'
	end

c*************************************************************

        subroutine read_time_series(unit,nvar,b3d,time,values,ierr)

c time series read

        implicit none

        integer unit
        integer nvar
	logical b3d	!true if 3d read
        real time
        real values(nvar)
        integer ierr

	if( b3d ) then
          call read_3_time_series(unit,nvar,time,values,ierr)
	else
          call read_0_time_series(unit,nvar,time,values,ierr)
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

        subroutine read_3_time_series(unit,nvar,time,values,ierr)

c 3D time series read:  
c       time,lmax,nk,nvars
c       values

        implicit none

        integer unit
        integer nvar
        real time
        real values(nvar)
        integer ierr

        integer j,n,lmax,nk,nvars
        integer ip,ivar,k,kn
	logical bdebug

	bdebug = .true.
	bdebug = .false.

	if( bdebug ) write(6,*) '3d TS read : ',unit,nvar

	read(unit,*,iostat=ierr) time,lmax,nk,nvars
	if( bdebug ) write(6,*) time,lmax,nk,nvars
	n = 0
        if( ierr .gt. 0 ) goto 99
        if( ierr .lt. 0 ) return
        n = lmax * nk * nvars
        if( n .ne. nvar ) goto 99

        ip = 0
        do ivar=1,nvars
          do k=1,nk
	    read(unit,*,iostat=ierr) kn,(values(ip+j),j=1,lmax)
	    !write(6,*) k,kn,(values(ip+j),j=1,lmax)
            if( kn .ne. k ) goto 98
            ip = ip + lmax
          end do
        end do

	if( bdebug ) write(6,*) '3d TS read (last ip) : ',unit,ip

        return
   98   continue
        write(6,*) 'read error on unit: ',unit
        write(6,*) time,lmax,nk,nvars
        write(6,*) k,kn
        write(6,*) '(node number not compatible)'
        stop 'error stop read_3_time_series: nk'
   99   continue
        write(6,*) 'read error on unit: ',unit
        write(6,*) time,lmax,nk,nvars
        write(6,*) n,nvar
        write(6,*) '(data set may not be compatible)'
        write(6,*) 'calling routine wants ',nvar,' data'
        write(6,*) 'I could read          ',n,' data'
        write(6,*) 'There could be a mismatch between nbdim and lmax'
        stop 'error stop read_3_time_series: nvar'
        end

c*************************************************************

	function exxqq9(iunit,nintp,rx,x,y)

c this routine should be substituted by a more general one
c
c interpolation directly from file (more files contemporaneous)
c
c interpolation of values read directly from formatted file
c records of file must be written as : write(iunit,*) x,y
c x values must be in increasing order
c rx values for which interpolation has to be performed
c ...must be called for in increasing order
c works also for one point interpolation
c
c iunit		file number from which data is read
c nintp		number of points used for interpolation
c		...if negative : reset file
c rx		x value for which interpolated y value has to be found
c x,y		arrays x,y read by subroutine but to be conserved
c		...by calling procedure
c exxqq9	interpolated y value (return value)
c
	implicit none
c
c arguments
	real exxqq9
	integer iunit,nintp
	real rx
        real x(1),y(1)
c local
	integer ninter,i
	integer ncomp1,ncomp2
	real xcomp,xh,yh
	character*70 name
c functions
	real exxp
c save
	logical bdebug
	real eps
	save eps,bdebug
	data eps / 1.e-5 /
	data bdebug / .false. /
c
	ninter=nintp
c
	if(ninter.lt.0) then	!reset file
		ninter=-ninter
		rewind(iunit)
		do i=1,ninter
		  read(iunit,*,err=95,end=94) x(i),y(i)
		  if( bdebug ) write(6,*) x(i),y(i)
		end do
		do i=2,ninter
		  if( x(i) .le. x(i-1) ) goto 88
		end do
	else if(ninter.eq.0) then
		goto 90
	end if
c
	if(mod(ninter,2).eq.0) then	!even
		ncomp1=1+ninter/2
		ncomp2=ncomp1
	else				!odd
		ncomp1=1+ninter/2
		if(ninter.gt.1) then
			ncomp2=ncomp1+1
		else
			ncomp2=ncomp1
		end if
	end if
c
	xcomp=(x(ncomp1)+x(ncomp2))/2.
c
	do while(rx.gt.xcomp)
		read(iunit,*,end=1,err=95) xh,yh
		if( xh .le. x(ninter) ) goto 88
		if( bdebug ) write(6,*) 'new___   ',xh,yh
		do i=1,ninter-1
		   x(i)=x(i+1)
		   y(i)=y(i+1)
		   if( bdebug ) write(6,*) x(i),y(i)
		end do
		x(ninter)=xh
		y(ninter)=yh
		if( bdebug ) write(6,*) x(ninter),y(ninter)
		xcomp=(x(ncomp1)+x(ncomp2))/2.
	end do
c
    1	continue
c
	if(ninter.gt.1) then	!no check for ninter = 1
		if(rx.lt.x(1)-eps) goto 92
		if(rx.gt.x(ninter)+eps) goto 91
	end if
c
	if( bdebug ) then
	  do i=2,ninter
	    if( x(i) .le. x(i-1) ) goto 88
	  end do
	end if

	exxqq9=exxp(ninter,x,y,rx)
c
	return
   88	continue
	write(6,*) 'x values are not in increasing order:'
	write(6,*) (x(i),i=1,ninter)
	write(6,*) (y(i),i=1,ninter)
	write(6,*) 'interpolation grade = ',ninter
	write(6,*) 'rx = ',rx
	write(6,*) 'unit = ',iunit
	call filna(iunit,name)
	write(6,'(a,a)') 'file = ',name
	stop 'error stop : exxqq9'
   90	continue
	write(6,*) 'Value for ninter not allowed : ',ninter
	stop 'error stop : exxqq9'
   91	continue
	write(6,*) 'No extrapolation possible'
	write(6,*) 'Actual x level :',rx
	write(6,*) 'Last available time levels :'
	write(6,*) (x(i),i=1,ninter)
	stop 'error stop : exxqq9'
   92	continue
	write(6,*) 'No extrapolation possible'
	write(6,*) 'Actual x level :',rx
	write(6,*) 'First available time levels :'
	write(6,*) (x(i),i=1,ninter)
	stop 'error stop : exxqq9'
   94	continue
	write(6,*) 'End of file detected on unit :',iunit
	write(6,*) 'File must consist of at least ',ninter,' data'
	stop 'error stop : exxqq9'
   95	continue
	write(6,*) 'Read error from file :',iunit
	write(6,*) 'close to :'
	write(6,*) (x(i),i=1,ninter)
	write(6,*) 'Attention : Cannot read unformatted file'
	stop 'error stop : exxqq9'
	end
c
c*********************************************************
 
	function exxpp(nintp,nmax,x,y,xe,iact)

c interpolation from array
c
c from given values x,y a value ye corresponding
c to xe is interpolated. a cubic interpolation is used
c
c the program is looking for the closest x-value
c only in foreward direction. for this reason
c xe-values must be passed in an increasing sequence
c
c the program does not check, if the value of xe
c is in the bounds of x(1) - x(nmax)
c eventually an extrapolated value is returned in exxpp

	implicit none

	real exxpp	!extrapolated values
	integer nintp	!number of values to use (4 for cubic, 2 for linear)
	integer nmax	!length of data arrays
        real x(1),y(1)	!data arrays
	real xe		!x-value for which y-value has to be interpolated
	integer iact	!element closest to xe (of last call on entry)
			!must be 0 for initialization

	logical bdebug
	integer nanf,nend,i
	integer min,max
	real xlow,xhigh
	real ye

	real exxp

	bdebug = .true.
	bdebug = .false.

c----------------------------------------------------------
c start searching from first element in x
c----------------------------------------------------------

        if( iact .le. 0 ) iact=1

c----------------------------------------------------------
c find to xe closest x-value in vector x
c----------------------------------------------------------

	do while( iact .lt. nmax )
          xlow  = abs(x(iact)-xe)
          xhigh = abs(x(iact+1)-xe)
          if( xhigh .ge. xlow ) goto 1
          iact = iact + 1
	end do
    1   continue

c----------------------------------------------------------
c x(iact) is closest value to xe ...now get closest points around xe
c----------------------------------------------------------

	if( mod(nintp,2) .eq. 0 ) then	!even
		max = nintp / 2
		min = max - 1
	else
		max = nintp / 2
		min = max
	end if

        if( x(iact) .gt. xe ) then
                nanf=iact-max
                nend=iact+min
        else
                nanf=iact-min
                nend=iact+max
        end if

c----------------------------------------------------------
c handling for the beginning or the end of array x
c----------------------------------------------------------

        if( nanf .lt. 1 ) then
                nanf=1
                nend=nintp
        else if(nend.gt.nmax) then
                nanf=nmax-nintp+1
                nend=nmax
        end if

c----------------------------------------------------------
c interpolation
c----------------------------------------------------------

        ye=exxp(nintp,x(nanf),y(nanf),xe)

c----------------------------------------------------------
c debug
c----------------------------------------------------------

	if( bdebug ) then
	  write(6,*) '-------------------------------- debug exxpp'
	  write(6,*) iact,nanf,nend,nintp,nmax
	  write(6,*) (x(i),i=nanf,nend)
	  write(6,*) (y(i),i=nanf,nend)
	  write(6,*) xe,ye
	  write(6,*) '--------------------------------'
	end if

c----------------------------------------------------------
c in ye is interpolated value
c----------------------------------------------------------

        exxpp=ye
	
        end

c***************************************************************
c***************************************************************
c***************************************************************

	function rintq(z,idim,jdim,ipos,jpos,xdelta,ydelta,switch,ier)
c
c interpolation in a square
c
c z 		matrix containing the values at the nodes
c (idim,jdim)	dimension of z
c ipos,jpos	position of local node (0,0) in z
c xdelta	relative x-coordinate (0...1) of interpolation point in square
c ydelta	   "     y-    "         "            "          "        "
c		...(e.g. : xdelta=0.5=ydelta is centre of square)
c switch	value at node that has not to be used for interpolation
c		...these values are extrapolated from the other nodes
c		...in case all 4 nodes equal to switch ==> rintq = switch
c ier		error status (return value)
c		... 0 : ok
c		... 1 : coordinates have been adjusted because out of grid
c rintq		interpolated value at return
c
c numeration of square
c			  (0,1)            (1,1)
c				+--------+
c				| 3    4 |
c				|        |
c				| 1    2 |
c				+--------+
c			  (0,0)            (1,0)
c
c diagonal sum = 5
c
c formula for interpolation :
c
c	z = a + b*x + c*y + d*x*y	resulting for a square in
c
c	a = z(0,0)
c	b = z(1,0) - z(0,0)
c	c = z(0,1) - z(0,0)
c	d = z(1,1) + z(0,0) - z(1,0) - z(0,1)
c
	real z(idim,jdim)
c
	integer iv(4),jv(4),iout(4),iin(4)
	real zh(4)
	data iv,jv /0,1,0,1,0,0,1,1/	!translates matrix into vector coord.
c
	i0=ipos
	j0=jpos
	xd=xdelta
	yd=ydelta
c
c ajust coordinates out of grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	ier=0
	if(i0.lt.1) then
		ier=1
		i0=1
		xd=0.
	else if(i0.ge.idim) then
		ier=1
		i0=idim-1
		xd=1.
	end if
	if(j0.lt.1) then
		ier=1
		j0=1
		yd=0.
	else if(j0.ge.jdim) then
		ier=1
		j0=jdim-1
		yd=1.
	end if
c
c find points with good/bad values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	nout=0
	nin=0
	ztot=0.
	do i=1,4
	zhh=z(i0+iv(i),j0+jv(i))
	zh(i)=zhh
	if(zhh.eq.switch) then
		nout=nout+1
		iout(nout)=i
	else
		nin=nin+1
		iin(nin)=i
		ztot=ztot+zhh
	end if
	end do
c
c extrapolate good (inner) to bad (outer) points %%%%%%%%%%%%%%%%%
c
	if(nout.eq.0) then	!no outer points
c		nothing
	else if(nin.eq.0) then	!no inner point
		rintq=switch
		return
	else if(nin.eq.1) then	!only 1 inner point
		zhh=zh(iin(1))
		do ii=1,4
		zh(ii)=zhh
		end do
	else if(nout.eq.1) then	!only 1 outer point
		iih=iout(1)
		idiag=5-iih
		zh(iih)=ztot-2.*zh(idiag)	!extrapolation from inner
						!...triangel to ext. point
	else			!2 inner/outer points
		if(iin(1)+iin(2).eq.5) then	!diagonal out of area
			zzh=(zh(iin(1))+zh(iin(2)))*0.5
			zh(iout(1))=zzh
			zh(iout(2))=zzh
		else				!side out of area
			iih=5-iin(2)		!to find point to be extrapol.
			zh(iih)=zh(iin(1))	!...get second inner point
			iih=5-iin(1)		!...and go to the diagonal
			zh(iih)=zh(iin(2))	!...
		end if
	end if
c
c interpolation in square %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	a=zh(1)
	b=zh(2)-a
	c=zh(3)-a
	d=zh(4)-a-b-c
c
	zzz = a + b*xd + c*yd + d*xd*yd
c
	rintq=zzz
c
	return
	end

c***************************************************************
c***************************************************************
c***************************************************************
c
c exffil	opens file and sets up array
c exffils	short version of exffil
c exfini	initializes array with file already open
c exfintp	interpolates value
c exfinfo	info on array
c
c all information needed is stored in vector array(*)
c
c used variable:
c
c       iunit           unit number of file
c       nintp           degree of interpolation (2 linear, 4 cubic)
c       nvar            number of variables
c       nread           read mode (dimension: 0-3)
c       ndim            dimension of array
c
c       nextra          extra information at beginning of array
c       nspace          space needed to hold all information in array
c       ires            pointer into array to keep interpolated values
c       rguard          guard value to check for out-of-bound access
c       
c formula for computing needed space:
c	nspace = 1 + nextra + (nintp+1) * (nvar+1)
c
c       1               guard value
c       nextra          extra information
c       nvar+1          to keep variables (nvar) and time
c       nintp+1         to keep nintp time steps and one for result
c
c formula for computing pointer ires:
c       ires   = 1 + nextra + nintp * (nvar+1)
c
c filling of array:
c
c       iunit,nintp,nvar,nread			(nextra)
c	var_0(i)		i=1,nintp	(nintp)	 - time
c	var_1(i)		i=1,nintp	(nintp)	 - first variable
c       ...
c	var_nvar(i)		i=1,nintp	(nintp)	 - last variable
c       t_intp,var_intp(j)      j=1,nvar 	(nvar+1) - at interpolated time
c       rguard					(1)
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

        subroutine exffil(file,nintp,nvar,nread,ndim,array)

c opens file and inititializes array
c
c everything needed is in array (unit, vars etc...)

        implicit none

        character*(*) file      !file name
	integer nintp		!grade of interpolation (2=linear,4=cubic)
	integer nvar		!how many columns to read/interpolate
        integer nread           !read mode (0-3)
        integer ndim            !dimension of array
        real array(ndim)        !array with all information

        integer iunit
        integer ifileo

	if( file .ne. ' ' ) then
          iunit = 55
	  !write(6,*) 'opening file : ',file
          iunit = ifileo(iunit,file,'form','old')
	  if( iunit .le. 0 ) stop 'error stop exffil: no such file'
	else
	  iunit = -1
	end if

        call exfini(iunit,nintp,nvar,nread,ndim,array)

        end

c***************************************************************

        subroutine exffils(file,ndim,array)
	implicit none
        character*(*) file      !file name
        integer ndim            !dimension of array
        real array(ndim)        !array with all information
        integer nintp           !grade of interpolation (2=linear,4=cubic)
        integer nvar            !how many columns to read/interpolate
        integer nread           !read mode (0-3)
	nintp=2
	nvar=1
        nread=0
	call exffil(file,nintp,nvar,nread,ndim,array)
	end

c***************************************************************

	subroutine exfini(iunit,nintp,nvar,nread,ndim,array)

c sets up interpolation from file -> all information is in array
c
c       space is computed as follows:
c
c       nvar variables + time -> nvar+1
c       nintp values per variable + interpolated values -> nintp+1
c       nextra is extra information heading the other data
c       one guard value at end of array

	implicit none

	integer iunit		!unit of file
	integer nintp		!grade of interpolation (4 for cubic)
	integer nvar		!how many columns to read/interpolate
	integer nread
	integer ndim		!dimension of array, on return space used
	real array(ndim)	!array with information

        real rguard
        parameter(rguard=1.234543e+20)
	integer nextra
        parameter(nextra=4)
	integer nxdim
        parameter(nxdim=1000)

	logical b3d,debug
	integer nspace,nnintp
        integer i
	real auxv(nxdim)    !aux - not used
	real time
	character*80 file

	debug = .false.

	if( nvar .gt. nxdim ) then
	  write(6,*) 'nvar too big: ',nvar,nxdim
	  write(6,*) 'please increment nxdim in exfini'
	  stop 'error stop exfini: nxdim'
	end if

        nnintp = nintp
	if( iunit .le. 0 ) nnintp = 0   !reserve some space anyway

	b3d = nread .gt. 0
	if(debug) write(6,*) 'exfini: ',nextra,nnintp,nvar

	nspace = 1 + nextra + (nnintp+1) * (nvar+1)

	if( nspace .gt. ndim ) then
                write(6,*) '*** error in exfini'
		write(6,*) 'Space in array is not sufficient'
		write(6,*) 'dimension of array: ',ndim
		write(6,*) 'space needed:       ',nspace
		write(6,*) 'formula: 1 + nextra + (nnintp+1) * (nvar+1)'
		write(6,*) 'with nextra = ',nextra
		write(6,*) 'iunit,nnintp,nvar:  ',iunit,nnintp,nvar
		write(6,*) 'nread:  ',nread
		call filna(iunit,file)
		write(6,*) 'filename: ',file
		stop 'error stop exfini: ndim'
	end if

        do i=1,nspace
          array(i) = 0.
        end do

	if(debug) write(6,*) 'exfini: (initializing) ',iunit,nspace

	if( iunit .gt. 0 ) then
	  time = 0.     !is not used
	  call intp_ts(-iunit,nnintp,nvar,time,array(nextra+1),auxv,b3d)
        end if

	if(debug) then
	  write(6,*) 'exfini (end of intp_ts) : ',iunit,nintp,nnintp,nvar
	end if

	array(1) = iunit
	array(2) = nnintp
	array(3) = nvar
	array(4) = nread
	array(nspace) = rguard

	if(debug) write(6,*) 'exfini: (finished initializing) ',iunit,nvar

	end

c***************************************************************

	subroutine exfintp(array,t,rint)

c interpolation from file -> all information is in array
c
c interpolated values are in last part of array

	implicit none

	real array(9)		!array with information from set-up
	real t			!t value for which to interpolate
	real rint(1)		!interpolated values

        real rguard
        parameter(rguard=1.234543e+20)
	integer nextra
        parameter(nextra=4)

	logical b3d
	integer iunit,nintp,nvar,nread
        integer ires,nspace,i

	iunit = nint(array(1))
	nintp = nint(array(2))
	nvar  = nint(array(3))
	nread = nint(array(4))

	b3d = nread .gt. 0

        ires   = 1 + nextra + nintp * (nvar+1)
	nspace = 1 + nextra + (nintp+1) * (nvar+1)

	if( iunit .gt. 0 ) then
	  !call exxqq(iunit,nintp,nvar,t,array(nextra+1),rint)
	  call intp_ts(iunit,nintp,nvar,t,array(nextra+1),rint,b3d)
          array(ires) = t
          do i=1,nvar
            array(ires+i) = rint(i)
          end do
	end if

        if( array(nspace) .ne. rguard ) then
          stop 'error stop exfintp: guard value altered'
        end if

	end

c***************************************************************

	subroutine exfget(array,t,rint)

c get last interpolated values

	implicit none

	real array(3)		!array with information from set-up
	real t			!t value for which to interpolate
	real rint(1)		!interpolated values

        real rguard
        parameter(rguard=1.234543e+20)
	integer nextra
        parameter(nextra=4)

	integer iunit,nintp,nvar
        integer ires,nspace,i

	iunit = nint(array(1))
	nintp = nint(array(2))
	nvar  = nint(array(3))

        ires   = 1 + nextra + nintp * (nvar+1)
	nspace = 1 + nextra + (nintp+1) * (nvar+1)

        t = array(ires)
        do i=1,nvar
          rint(i) = array(ires+i)
        end do

        if( array(nspace) .ne. rguard ) then
          stop 'error stop exfintp: guard value altered'
        end if

	end

c***************************************************************

	subroutine exfset(array,t,rint)

c sets new actual values

	implicit none

	real array(5)		!array with information from set-up
	real t			!t value for which to interpolate
	real rint(1)		!interpolated values

        real rguard
        parameter(rguard=1.234543e+20)
	integer nextra
        parameter(nextra=4)

	integer iunit,nintp,nvar
        integer ires,nspace,i

	iunit = nint(array(1))
	nintp = nint(array(2))
	nvar  = nint(array(3))

        ires   = 1 + nextra + nintp * (nvar+1)
	nspace = 1 + nextra + (nintp+1) * (nvar+1)

        array(ires) = t
        do i=1,nvar
          array(ires+i) = rint(i)
        end do

        if( array(nspace) .ne. rguard ) then
          stop 'error stop exfintp: guard value altered'
        end if

	end

c***************************************************************

	subroutine exfinfo(ipunit,array)

c interpolation from file -> info

	implicit none

	integer ipunit		!unit where to print on (<0 -> 6)
	real array(9)		!array with information from set-up

        real rguard
        parameter(rguard=1.234543e+20)
	integer nextra
        parameter(nextra=4)

	logical bdebug
	integer iunit,nintp,nvar,nread
        integer ires,nspace,i
	integer ipu
	integer ip,np,ivar

	bdebug = .false.
	bdebug = .true.

	ipu = ipunit
	if( ipunit .le. 0 ) ipu = 6

	iunit = nint(array(1))
	nintp = nint(array(2))
	nvar  = nint(array(3))
	nread = nint(array(4))

        ires   = 1 + nextra + nintp * (nvar+1)
	nspace = 1 + nextra + (nintp+1) * (nvar+1)

        write(ipu,*) 'info on array interpolation:'
        write(ipu,*) 'unit     : ',iunit
        write(ipu,*) 'nintp    : ',nintp
        write(ipu,*) 'nvar     : ',nvar
        write(ipu,*) 'nread    : ',nread
        write(ipu,*) 'ires     : ',ires
        write(ipu,*) 'nspace   : ',nspace
        write(ipu,*) 'rguard   : ',array(nspace)

	if( bdebug ) then
	  np = 1 + nextra
	  do ivar=0,nvar
	    ip = np + ivar*nintp
	    write(ipu,*) ivar,ip
	    write(ipu,*) (array(ip-1+i),i=1,nintp)
	  end do
	  write(ipu,*) (array(ires+i),i=0,nvar)
	end if

        end

c***************************************************************

	subroutine exf_set_nread(array,nread)

c sets dimension of read variable (only 0 or 3)

	implicit none

	real array(9)		!array with information from set-up
        integer nread           !read dimension (0-3)

        if( nread .lt. 0 .and. nread .gt. 3 ) then
          write(6,*) nread,'  (olny 0-3 allowed)'
          stop 'error stop exf_set_nread: nread'
        end if

	array(4) = nread

        end

c***************************************************************

	subroutine exfunit(array,iunit)

c returns information on unit number

	implicit none

	real array(5)		!array with information from set-up
        integer iunit           !unit number of file, 0 if not initialized

	iunit = nint(array(1))

        end

c***************************************************************

	subroutine exfnvar(array,nvar)

c returns information on number of variables

	implicit none

	real array(5)		!array with information from set-up
        integer nvar            !total number of variables

	nvar = nint(array(3))

        end

c***************************************************************

	subroutine exfpres(array,ires)

c returns pointer to interpolated result
c
c array(ires)	actual time of interpolated results
c array(ires+i)	actual interpolated result of variable i

	implicit none

	real array(5)		!array with information from set-up
        integer ires            !pointer to results in array

	integer nextra
        parameter(nextra=4)

	integer nintp,nvar

	nintp = nint(array(2))
	nvar  = nint(array(3))
        ires   = 1 + nextra + nintp * (nvar+1)

        end

c***************************************************************

