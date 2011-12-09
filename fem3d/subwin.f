c
c $Id: subwin.f,v 1.29 2010-02-26 17:35:06 georg Exp $
c
c wind routines
c
c contents :
c
c subroutine windad				administers read of wind file
c
c routines for unformatted wind file
c
c subroutine windrd(iunit,it,n,u,v,p,ierr)	shell for unformatted read
c subroutine wrwin(iunit,it,n,u,v,p)		writes unformatted wind record
c subroutine rdwin(iunit,it,n,u,v,p,ierr)	reads unformatted wind record
c
c routines for read of STR section $wind
c
c subroutine inwnds
c subroutine rdwnds
c subroutine ckwnds
c subroutine prwnds
c subroutine tswnds
c
c function nwndget()				number of wind records avail.
c subroutine nwndset( n )			sets number of wind records
c
c subroutine gtwnds(it,u,v,p,nkn,ierr)		gets wind data from section
c
c function binwin( file , nkn )			test if wind file is unform.
c subroutine rdwnda( nwin )			reads ASCII wind data from file
c subroutine rdwndl( line )			reads one line of wind data
c
c subroutine wstress(n,tx,ty)			converts wind to norm. stress
c
c notes :
c
c use new routines ($wind section) with parameter iwtype > 0
c
c if wind file is given read from file
c
c iwtype = 0 no wind data processed
c iwtype = 1 wind data in format (wx,wy)
c iwtype = 2 wind data in format (tx,ty) (stress)
c iwtype = 3 wind data in format (speed,dir) (speed in m/s)
c		dir = 0  -> wind from north
c		dir = 90 -> wind from east
c iwtype = 4 wind data in format (speed,dir) (speed in nodes)
c
c formulas:
c
c here some info on formulas used
c
c formats:
c
c first unformatted file is read
c then formatted file is tried
c finally section in STR file is used
c
c revision log :
c
c 25.06.1997	ggu	written and assembled from subn11.f
c 22.01.1998	ggu	better warnings
c 25.05.1998	ggu	more documentation, stress implemented
c 26.05.1998	ggu	ASCII format implemeted
c 01.06.1998	ggu	ignore after 3rd value for ASCII read
c 17.06.1998	ggu	no debug write
c 13.07.1998	ggu	compiler error: define nwndget as nwndget()
c 14.07.1998	ggu	function binwin rewritten (gave no error on linux)
c 22.07.1998	ggu	documentation
c 23.07.1998	ggu	documentation
c 11.11.1998	ggu	nwdim increased to 200
c 20.04.1999    ggu     converted to stress instead of wind (tauxnv...)
c 31.05.1999    ggu     bug fix -> n was not set to nkn (twice)
c 08.07.1999    ggu     bug fix in documentation
c 05.04.2000	ggu	nwdim increased to 1000
c 03.12.2001	ggu	write only first 100 records of wind data to stdout
c 09.01.2003	ggu	INTEL unassigned itw
c 03.11.2003	ggu	bug in documentation (stress without rho_0)
c 16.09.2004	ggu	new routine convert_wind()
c 09.02.2010	ggu	less diagnostics, only write wind data if read from STR
c 22.02.2010	ggu&ccf	restructured, new formulas Smith&Banke, Large&Pond
c 16.02.2011	ggu	set std pressure if pressure not given
c 25.02.2011	ggu	new param wsmax to catch errors in wind type
c 23.03.2011	ggu	better error message in binwin()
c 18.08.2011	ggu	bug fix in binwin() -> close file
c
c*************************************************************************

c DOCS	START	S_wind
c
c In this section the wind data can be given directly without
c the creation of an external file. Note, however, that
c a wind file specified in the |name| section takes precedence
c over this section. E.g., if both a section |wind| and a
c wind file in |name| is given, the wind data from the file is used.
c
c The format of the wind data in this section is the same as the
c format in the ASCII wind file, i.e., three columns, with
c the first specifying the time in seconds and the other two columns
c giving the wind data. The interpretation of the wind data
c depends on the value of |iwtype|. For more information please
c see the description of |iwtype| in section |para|.
c
c DOCS	END

c DOCS	START	P_wind
c
c DOCS	WIND		Wind parameters
c
c The next two parameters deal with the wind stress to be
c prescribed at the surface of the basin.
c
c The wind data can either be specified in an external file (ASCII
c or binary) or directly in the parameter file in section |wind|.
c The ASCII file or the wind section contain three columns, the first
c giving the time in seconds, and the others the components of
c the wind speed. Please see below how the last two columns are
c interpreted depending on the value of |iwtype|. For the format
c of the binary file please see the relative section.
c If both a wind file and section |wind| are given, data from the
c file is used.
c
c The wind stress is normally computed with the following formula
c \beq
c \tau^x = \rho_a c_D \vert u \vert u^x \quad
c \tau^y = \rho_a c_D \vert u \vert u^y 
c \eeq
c where $\rho_a,\rho_0$ is the density of air and water respectively,
c $u$ the modulus of wind speed and $u^x,u^y$ the components of
c wind speed in $x,y$ direction. In this formulation $c_D$ is a
c dimensionless drag coefficient that varies between 1.5 \ten{-3} and
c 3.2 \ten{-3}. The wind speed is normally the wind speed measured
c at a height of 10 m.
c
c |iwtype|	The type of wind data given (default 1):
c		\begin{description}
c		\item[0] No wind data is processed
c		\item[1] The components of the wind is given in [m/s]
c		\item[2] The stress ($\tau^x,\tau^y$) is directly specified
c		\item[3] The wind is given in speed [m/s] and direction
c			 [degrees]. A direction of 0\degrees{} specifies
c			 a wind from the north, 90\degrees{} a wind
c			 from the east etc.
c		\item[4] As in 3 but the speed is given in knots
c		\end{description}
c |itdrag|	Formula to compute the drag coefficient. A value of 0
c		uses the constant value given in |dragco|. With 1
c		the Smith and Banke formula is used.
c |dragco|	Drag coefficient used in the above formula. The default value
c		is 0 so it must be specified. Please note also that in case
c		of |iwtype| = 2 this value is of no interest, since the
c		stress is specified directly.
c |wsmax|	Maximum wind speed allowed in [m/s]. This is in order to avoid
c		errors if the wind data is given in a different format
c		from the one spwecified by |iwtype|. (Default 50)
c
c DOCS	END

c*************************************************************************
c*************************************************************************
c administration of wind data
c*************************************************************************
c*************************************************************************

	subroutine windad

c administers read of wind file

	implicit none

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

        real tauxnv(1),tauynv(1)
        common /tauxnv/tauxnv,/tauynv/tauynv
        real wxv(1),wyv(1)
        common /wxv/wxv,/wyv/wyv
        real wxov(1),wyov(1)
        common /wxov/wxov,/wyov/wyov
        real wxnv(1),wynv(1)
        common /wxnv/wxnv,/wynv/wynv
        real ppv(1)
        common /ppv/ppv
        real pov(1), pnv(1)
        common /pov/pov, /pnv/pnv

        real metws(1)
        common /metws/metws

	character*80 name
	real rit
	integer n,i,itw,ierr

	logical binwin
	integer ifileo,iround
	integer nwndget
	real getpar

	logical bfile		!true if reading from file
	logical	bstress		!true if wind data in stress
	save bfile,bstress
	integer icall,nbwin,itwold,itwnew,iwtype
	save icall,nbwin,itwold,itwnew,iwtype
	data icall / 0 /

c-----------------------------------------------------------------
c if no wind data return
c-----------------------------------------------------------------

	if( icall .eq. -1 ) return

c-----------------------------------------------------------------
c initialization
c-----------------------------------------------------------------

	if( icall .eq. 0 ) then

	    do i=1,nkn
	        tauxnv(i) = 0.
	        tauynv(i) = 0.
	        ppv(i) = 0.
	    end do

c           ------------------------------------------------------
c	    find out what data to read
c           ------------------------------------------------------

	    iwtype = iround(getpar('iwtype'))
	    bstress = iwtype .eq. 2
            call getfnm('wind',name)
	    bfile = name .ne. ' '

	    if( iwtype .le. 0 ) icall = -1
	    if( nwndget() .le. 0 .and. .not. bfile ) icall = -1
	    if( icall .eq. -1 ) return

c           ------------------------------------------------------
c	    if file is present open and read it
c           ------------------------------------------------------

	    if( bfile ) then
	      if( nwndget() .gt. 0 ) goto 89
	      if( binwin(name,nkn) ) then	!unformatted original format
	        if( iwtype .gt. 2 ) goto 91
	        nbwin=ifileo(0,name,'unform','old')
	        if(nbwin.le.0) goto 93
	        write(6,*) 'reading binary wind-file :'
	        write(6,*) name
	      else				!formatted new format
	        nbwin=ifileo(0,name,'form','old')
	        if(nbwin.le.0) goto 93
	        write(6,*) 'reading ASCII wind-file :'
	        write(6,*) name
		call rdwnda(nbwin)		!read in all data
		close(nbwin)
		bfile = .false.			!fake no read from file anymore
	      end if
	    end if

c           ------------------------------------------------------
c	    read first wind record
c           ------------------------------------------------------

	    n=nkn
	    if( bfile ) then					!unformatted
	      call rdwin(nbwin,itw,n,wxov,wyov,pov,ierr)
	    else						!formatted
	      call gtwnds(itw,wxov,wyov,pov,nkn,ierr)
	    end if
	    if(ierr.ne.0) goto 90

	    itwold=itw
	    itwnew=itw

c           ------------------------------------------------------
c	    copy wind field read to new wind records
c           ------------------------------------------------------

	    do i=1,nkn
	       wxnv(i)=wxov(i)
	       wynv(i)=wyov(i)
	       pnv(i)=pov(i)
	    end do

c           ------------------------------------------------------
c	    if no wind field for time it present use first one
c           ------------------------------------------------------

	    if(itanf.lt.itwold) then
		write(6,*) 'No wind field for time ',itanf
		write(6,*) 'Wind field used from time ',itwold
		itwold=itanf
	    end if
	end if

c-----------------------------------------------------------------
c normal call
c-----------------------------------------------------------------

	icall = icall + 1

c-----------------------------------------------------------------
c loop for next wind field
c-----------------------------------------------------------------

	do while(it.gt.itwnew)

c           ------------------------------------------------------
c	    copy old wind field 
c           ------------------------------------------------------

	    itwold=itwnew
	    do i=1,nkn
	       wxov(i)=wxnv(i)
	       wyov(i)=wynv(i)
	       pov(i)=pnv(i)
	    end do

c           ------------------------------------------------------
c	    read from file (unformatted) or from array (formatted)
c           ------------------------------------------------------

	    n=nkn
	    if( bfile ) then					!unformatted
	      call rdwin(nbwin,itw,n,wxnv,wynv,pnv,ierr)
	    else						!formatted
	      call gtwnds(itw,wxnv,wynv,pnv,nkn,ierr)
	    end if

c           ------------------------------------------------------
c	    no more data -> use last wind field read
c           ------------------------------------------------------

	    if(ierr.lt.0) then					!no more data
	    	do i=1,nkn
	    	   wxnv(i)=wxov(i)
	    	   wynv(i)=wyov(i)
	    	   pnv(i)=pov(i)
	    	end do
	    	write(6,*) 'No wind field for time ',it
	    	write(6,*) 'Wind field used from time ',itwold
	    	write(6,*) 'until end of simulation'
	    	itwnew=itend
	    else if(ierr.gt.0) then
		goto 92
	    else
	        itwnew=itw					!INTEL
	    end if

	end do

c-----------------------------------------------------------------
c interpolation
c-----------------------------------------------------------------

c	write(6,*) it,itwold,itwnew

	if( itwnew .gt. itwold ) then			!HACK
	  rit=float(it-itwold)/float(itwnew-itwold)
	else
	  rit = 0
	end if

	do i=1,nkn
	  wxv(i)=rit*(wxnv(i)-wxov(i))+wxov(i)
	  wyv(i)=rit*(wynv(i)-wyov(i))+wyov(i)
	  ppv(i)=rit*(pnv(i)-pov(i))+pov(i)
	  metws(i) = sqrt( wxv(i)**2 + wyv(i)**2 )
	end do

c-----------------------------------------------------------------
c convert wind data to normalized stress
c-----------------------------------------------------------------

	call wstress(nkn,wxv,wyv,tauxnv,tauynv)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	return
   89   continue
	write(6,*) 'Both wind file and STR wind data given.'
	write(6,*) 'Please give wind data only in one section!'
	stop 'error stop windad: too many sources of wind data'
   90   continue
        write(6,*) 'No data in wind file :'
        write(6,*) name
        stop 'error stop : windad'
   91   continue
        write(6,*) 'wind file not possible with type : ',iwtype
        stop 'error stop : windad'
   92   continue
	write(6,*) 'Error reading wind file : ',ierr
	stop 'error stop : windad'
   93   continue
        write(6,*) 'error opening wind file :'
        write(6,*) name
        stop 'error stop : windad'
	end

c*************************************************************
c*************************************************************
c routines for unformatted file format
c*************************************************************
c*************************************************************

	function binwin(file,nkn)

c test if wind file is unformatted

	implicit none

	logical binwin
	integer nkn
	character*(*) file

	integer nwin,ios,it,n,na
	real u,v
	integer ifileo

	binwin = .false.

c-----------------------------------------------------------------
c try unformatted
c-----------------------------------------------------------------

	nwin = ifileo(0,file,'unformatted','old')
	if( nwin .le. 0 ) return		!error opening wind file

	read(nwin,iostat=ios) it,n

c	if there was an error we know that it is not a binary wind file
c	if not, we do some more checks, since the linux compiler
c	does not flag an error...

	if( ios .eq. 0 ) then
	    na = abs(n)
	    if( na .eq. nkn .or. na .eq. 1 ) then
		binwin = .true.			!no error -> unformatted
	  	write(6,*) 'the wind file is unformated: ',file
		close( nwin )
		return
	    end if
	end if

	close( nwin )

c-----------------------------------------------------------------
c try formatted
c-----------------------------------------------------------------

	nwin = ifileo(0,file,'formatted','old')
	if( nwin .le. 0 ) return		!error opening wind file

	read(nwin,*,iostat=ios) it,u,v

	if( ios .ne. 0 ) then
	  write(6,*) 'The wind file is not ASCII'
	  write(6,*) 'However there was an error reading it'
	  write(6,*) 'as unformatted.'
	  write(6,*) 'The total number of nodes could be not'
	  write(6,*) 'compatible with the basin file.'
	  write(6,*) 'nkn = ',nkn,'  nwind = ',na
	  stop 'error stop binwin: wind file not compatible'
	end if

	close( nwin )

	write(6,*) 'the wind file is ASCII: ',file

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c**************************************************************

        subroutine wrwin(iunit,it,n,u,v,p)

c writes unformatted wind record

        implicit none

	integer iunit		!unit number
	integer it		!time
	integer n		!dimension (in), negative if pressure is ok
	real u(1),v(1)		!wind speed
	real p(1)		!pressure

        integer i

        write(iunit) it,n

        if(n.gt.0) then
          write(iunit) (u(i),v(i),i=1,n)
        else
          write(iunit) (u(i),v(i),i=1,-n),(p(i),i=1,-n)
        end if

        end

c*************************************************************

        subroutine rdwin(iunit,it,ndim,u,v,p,ierr)

c reads unformatted wind record

        implicit none

	integer iunit		!unit number
	integer it		!time (out)
	integer ndim		!number of data to return
	integer ierr		!error code  0: ok  <0: EOF  >0: error
	real u(1),v(1)		!wind speed [m/s]
	real p(1)		!pressure [Pa]

	real pstd
	parameter ( pstd = 1013.25 )

        integer i,n
	logical bpress,bconst

c-----------------------------------------------------------------
c read first record
c-----------------------------------------------------------------

        read(iunit,iostat=ierr) it,n

	if( ierr .ne. 0 ) return

	bpress = n .lt. 0		!read pressure
	n = abs(n)
	bconst = n .eq. 1		!constant wind field

c-----------------------------------------------------------------
c error handling
c-----------------------------------------------------------------

	if( ierr .ne. 0 ) then		!read error or end of file
	  return
        else if(n.ne.ndim .and. n.ne.1) then	!wrong number of data
          read(iunit,iostat=ierr)       !dummy read
	  write(6,*) '*** rdwin: ',n,ndim
	  write(6,*) '*** not the value I was expecting'
	  write(6,*) 'n = ',n
	  write(6,*) 'expecting either 1 or ',ndim
          ierr=999
	  return
	end if

c-----------------------------------------------------------------
c read data record
c-----------------------------------------------------------------

        if( bpress ) then		!read also pressure
          read(iunit,iostat=ierr) (u(i),v(i),i=1,n),(p(i),i=1,n)
        else				!no pressure -> initialize it
          read(iunit,iostat=ierr) (u(i),v(i),i=1,n)
        end if

	if( ierr .ne. 0 ) return

c-----------------------------------------------------------------
c adjust data
c-----------------------------------------------------------------

        if( .not. bpress ) then
          do i=1,ndim
            p(i)=100.*pstd		!set to standard pressure (in Pa)
          end do
	end if

	if( bconst ) then
          do i=2,ndim
            u(i)=u(1)
            v(i)=v(1)
            p(i)=p(1)
          end do
	end if

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

        end

c**************************************************************
c**************************************************************
c routines for read of STR section $wind
c**************************************************************
c**************************************************************

	subroutine inwnds

	implicit none

	integer nwdim
	parameter(nwdim=100000)

	integer nwdi,nwind,iwtype,iwsect
	common /winpar/ nwdi,nwind,iwtype,iwsect
	real itwin(nwdim)
	common /wintim/ itwin
	real winx(nwdim), winy(nwdim), pres(nwdim)
	common /windat/ winx, winy, pres
	save /winpar/, /windat/, /wintim/

	nwdi = nwdim		!dimension of arrays
	nwind = 0		!total number of records
	iwtype = 0		!type of wind data
	iwsect = 0		!wind read from STR section ?

	end

c**************************************************************

	subroutine rdwnds

	implicit none

	integer nwdi,nwind,iwtype,iwsect
	common /winpar/ nwdi,nwind,iwtype,iwsect

	character*80 line
	integer nrdlin

        do while( nrdlin(line) .gt. 0 )
	   iwsect = iwsect + 1			!how many wind lines read
	   call rdwndl(line)
	end do

	end
	
c**************************************************************

	subroutine ckwnds

c checks iwtype and converts data
c
c speed [m/s] = rnodes * speed [nodes]

	implicit none

	real rnodes
	parameter( rnodes = 1852. / 3600. )	! n. miles / hour
	integer nwdim
	parameter(nwdim=100000)

	integer nwdi,nwind,iwtype,iwsect
	common /winpar/ nwdi,nwind,iwtype,iwsect
	real itwin(nwdim)
	common /wintim/ itwin
	real winx(nwdim), winy(nwdim), pres(nwdim)
	common /windat/ winx, winy, pres

	integer i
	real s,d,u,v
	real sfact

	integer iround
	real getpar

	if( nwdim .ne. nwdi ) stop 'ckwnds: error stop nwdim'

c-----------------------------------------------------------------
c check for correctness of iwtype
c-----------------------------------------------------------------

	iwtype = iround(getpar('iwtype'))

	if( iwtype .lt. 0 ) goto 98
	if( iwtype .eq. 0 ) return	!no wind data needed
	if( iwtype .eq. 1 ) return	!wind data in wx,wy format
	if( iwtype .eq. 2 ) return	!wind data in stress
	if( iwtype .gt. 4 ) goto 98

c-----------------------------------------------------------------
c convert from (s,d) to (wx,wy) - (iwtype=3,4)
c-----------------------------------------------------------------

	sfact = 1.
	if( iwtype .eq. 4 ) sfact = rnodes

	do i=1,nwind
	  s = sfact * winx(i)
	  d = winy(i)
          call convert_wind(s,d,u,v)
	  winx(i) = u
	  winy(i) = v
	end do

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	return
   98	continue
	write(6,*) 'error in wind type : ',iwtype
	stop 'error stop ckwnds'
	end

c**************************************************************

	subroutine prwnds

	implicit none

	integer nwdim
	parameter(nwdim=100000)

	integer nwdi,nwind,iwtype,iwsect
	common /winpar/ nwdi,nwind,iwtype,iwsect
	real itwin(nwdim)
	common /wintim/ itwin
	real winx(nwdim), winy(nwdim), pres(nwdim)
	common /windat/ winx, winy, pres

	integer i,n

	integer iround
	real getpar

	if( nwdim .ne. nwdi ) then
		write(6,*) 'prwnds: ',nwdi,nwind,iwtype
		stop 'prwnds: error stop nwdim'
	end if

	iwtype = iround(getpar('iwtype'))

	if( iwtype .le. 0 ) return	!read wind data from file
	if( iwsect .le. 0 ) return	!read wind data from file

        n = min(100,nwind)

	write(6,*)
	write(6,*) 'Wind data: ',nwdi,nwind,iwtype,iwsect
	do i=1,n
	  write(6,*) itwin(i),winx(i),winy(i),pres(i)
	end do

	end

c**************************************************************

	subroutine tswnds

	implicit none

	integer nwdi,nwind,iwtype,iwsect
	common /winpar/ nwdi,nwind,iwtype,iwsect

	write(6,*) 'tswnds: ',nwdi,nwind,iwtype,iwsect
	call prwnds

	end

c**************************************************************
c**************************************************************
c comodity routines for accessing wind array
c**************************************************************
c**************************************************************

	function nwndget()

c gets number of wind records available

	implicit none

	integer nwndget

	integer nwdi,nwind,iwtype,iwsect
	common /winpar/ nwdi,nwind,iwtype,iwsect

	nwndget = nwind

	end

c**************************************************************

	subroutine nwndset( n )

c sets number of wind records

	implicit none

	integer n

	integer nwdi,nwind,iwtype,iwsect
	common /winpar/ nwdi,nwind,iwtype,iwsect

	nwind = n

	end

c**************************************************************

	subroutine gtwnds(it,u,v,p,nkn,ierr)

c gets wind data from section (or ASCII file already read)

	implicit none

	integer it,nkn,ierr
	real u(1),v(1),p(1)

	integer nwdim
	parameter(nwdim=100000)

	integer nwdi,nwind,iwtype,iwsect
	common /winpar/ nwdi,nwind,iwtype,iwsect
	real itwin(nwdim)
	common /wintim/ itwin
	real winx(nwdim), winy(nwdim), pres(nwdim)
	common /windat/ winx, winy, pres

	integer i
	integer iwind
	save iwind
	data iwind / 0 /

	if( nwdim .ne. nwdi ) stop 'error stop nwdim'

c-----------------------------------------------------------------
c new record
c-----------------------------------------------------------------

	ierr = 0
	iwind = iwind + 1

c-----------------------------------------------------------------
c no more records ?
c-----------------------------------------------------------------

	if( iwind .gt. nwind ) then
		ierr = -1
		return
	end if

c-----------------------------------------------------------------
c copy to arrays
c-----------------------------------------------------------------

	it = itwin(iwind)

        do i=1,nkn
          u(i) = winx(iwind)
          v(i) = winy(iwind)
          p(i) = pres(iwind)
        end do

c-----------------------------------------------------------------
c write to terminal
c-----------------------------------------------------------------

        !write(6,*) 'wind field read (wx,wy,p):',it,u(1),v(1),p(1)

        !write(6,*) 'wind field read :',it
	!write(6,*) 'wx,wy,p : ',u(1),v(1),p(1)
        !write(6,*)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c**************************************************************
c**************************************************************
c reading formatted data
c**************************************************************
c**************************************************************

	subroutine rdwnda( nwin )

c reads ASCII wind data from file

	implicit none

	integer nwin

	integer nwdi,nwind,iwtype,iwsect
	common /winpar/ nwdi,nwind,iwtype,iwsect

	character*80 line

	nwind = 0	!reset

    1	continue
	   read(nwin,'(a)',end=2) line
	   call rdwndl(line)
	goto 1
    2	continue

	call ckwnds

	end

c**************************************************************

	subroutine rdwndl( line )

c reads one line of wind data

	implicit none

	integer nwdim
	parameter(nwdim=100000)
	integer ndim
	parameter(ndim=3)	!how many values to read on line

	integer nwdi,nwind,iwtype,iwsect
	common /winpar/ nwdi,nwind,iwtype,iwsect
	real itwin(nwdim)
	common /wintim/ itwin
	real winx(nwdim), winy(nwdim), pres(nwdim)
	common /windat/ winx, winy, pres

	real pstd
	parameter ( pstd = 1013.25 )

	character*(*) line
	integer ianz
	integer iscanf,iround
	real f(ndim)

	if( nwdim .ne. nwdi ) stop 'rdwndl: error stop nwdim'

	ianz = iscanf(line,f,ndim)

	if( ianz .lt. 0 ) goto 99
	if( ianz .eq. 0 ) return
	if( ianz .lt. 3 ) goto 98

	nwind = nwind + 1
	if( nwind .gt. nwdi ) goto 96
	itwin(nwind) = iround(f(1))
	winx(nwind) = f(2)
	winy(nwind) = f(3)
	pres(nwind) = 100. * pstd	!pressure in Pa

c	write(6,*) 'rdwndl: ',itwin(nwind),winx(nwind),winy(nwind)

	return
   96	continue
        write(6,*) nwdi,nwdim
	write(6,*) 'dimension error nwdim'
	stop 'error stop rdwndl'
   98	continue
	write(6,*) 'less than 3 values read on line'
	write(6,*) line
	stop 'error stop rdwndl'
   99	continue
	write(6,*) 'read error in line:'
	write(6,*) line
	stop 'error stop rdwndl'
	end
	
c**************************************************************
c**************************************************************
c utility routines
c**************************************************************
c**************************************************************

	subroutine wstress(n,wx,wy,tx,ty)

c converts wind data into normalized stress

c (1/rho_0) tau_x  ==  cd * (rho_air/rho_0) * |w| w_x
c
c stress in kg / ( m s**2 )
c
c rho_0		density water
c rho_air	density air
c cd		drag coefficient
c
c wk = cd * (rho_air/rho_0)
c
c wind   ==>          wc = wk * |w|
c stress ==>          wc = 1 / rho_0

	implicit none

	integer n
	real wx(1), wy(1)
	real tx(1), ty(1)

	integer k
	real roluft,rowass
	real aux,wxy,txy,wxymax

	real getpar

	logical bstress			!true if wind data is in stress
	integer iwtype,itdrag
	real dragco,wfact,wsmax
	save bstress,iwtype,itdrag
	save dragco,wfact,wsmax

	integer icall
	save icall
	data icall / 0 /

c-----------------------------------------------------------------
c initialization
c-----------------------------------------------------------------

	if( icall .eq. 0 ) then
	  iwtype = nint(getpar('iwtype'))
	  wsmax = getpar('wsmax')
	  bstress = iwtype .eq. 2
	  itdrag = nint(getpar('itdrag'))
	  dragco = getpar('dragco')
	  roluft = getpar('roluft')
	  rowass = getpar('rowass')
	  if( bstress ) then
	    wfact = 1. / rowass
	  else
	    wfact = roluft / rowass
	  end if
	end if

c-----------------------------------------------------------------
c normal call
c-----------------------------------------------------------------

	icall = icall + 1

c-----------------------------------------------------------------
c compute normalized stress
c-----------------------------------------------------------------

	wxymax = 0.

	if( bstress ) then		!data is stress -> normalize it
	  if( dragco .le. 0 ) dragco = 2.5E-3	!use standard value
	  do k=1,n
	    tx(k) = wfact * wx(k)
	    ty(k) = wfact * wy(k)
	    txy = sqrt( tx(k)**2 + ty(k)**2 )
	    wxy = sqrt(txy/dragco)
	    wxymax = max(wxymax,wxy)
	    wx(k) = tx(k) / (dragco*wxy)
	    wy(k) = ty(k) / (dragco*wxy)
	  end do
	else				!data is wind velocity [m/s]
	  do k=1,n
	    wxy = sqrt( wx(k)**2 + wy(k)**2 )
	    wxymax = max(wxymax,wxy)
	    if( itdrag .gt. 0 ) call get_drag(itdrag,wxy,dragco)
	    tx(k) = wfact * dragco * wxy * wx(k)
	    ty(k) = wfact * dragco * wxy * wy(k)
	  end do
	end if

c-----------------------------------------------------------------
c check wind speed
c-----------------------------------------------------------------

	if( wxymax .gt. wsmax ) then
	  write(6,*) 'maximum wind speed: ',wxymax
	  write(6,*) 'maximum allowed wind speed: ',wsmax
	  write(6,*) 'Are you sure the wind is in the correct format?'
	  write(6,*) 'If no, please set iwtype, else increase wsmax.'
	  stop 'error stop wstress: wind speed too high'
	end if

c-----------------------------------------------------------------
c debug output
c-----------------------------------------------------------------

c	call mima2i('Min/Max of wind stress (x/y) ...',n,tx,ty)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end
	
c**************************************************************

	subroutine get_drag(itdrag,wxy,dragco)

c computes drag coefficient

	implicit none

	integer itdrag		!type of formula
	real wxy		!wind speed
	real dragco		!computed drag coefficient

	if( itdrag .le. 0 ) then
	  !nothing
	else if( itdrag .eq. 1 ) then	!Smith and Banke (1975)
          dragco = 0.001 * (0.63 + 0.066*wxy)
	else if( itdrag .eq. 2 ) then	!Large and Pond (1981)
          if ( wxy .gt. 11. ) then
            dragco = 0.001 * (0.49 + 0.066*wxy)
          else
            dragco = 0.001 * 1.2
	  end if
	else
	  write(6,*) 'erroneous value for itdrag = ',itdrag
	  stop 'error stop get_drag: itdrag'
	end if
	
	end

c**************************************************************

	subroutine get_wind(k,wx,wy)

c helper function -> return wind for node k

	implicit none

	integer k
	real wx,wy

        real wxv(1),wyv(1)
        common /wxv/wxv,/wyv/wyv

	wx = wxv(k)
	wy = wyv(k)

	end

c**************************************************************

        subroutine convert_wind(s,d,u,v)

        implicit none

        real s,d,u,v

        real dir
        real pi,rad
        parameter(pi=3.14159,rad=pi/180.)

        dir = d
        dir = 90. - dir + 180.
        do while( dir .lt. 0. )
          dir = dir + 360.
        end do
        dir = mod(dir,360.)

        u = s*cos(rad*dir)
        v = s*sin(rad*dir)

        end

c**************************************************************

