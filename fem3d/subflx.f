c
c $Id: subflx.f,v 1.15 2007-03-20 13:14:53 georg Exp $
c
c utility routines to read/write FLX file - file type 537
c
c contents :
c
c	subroutine iniflx
c
c	subroutine rfflx	(iunit,nvers
c     +				,nscddi,nfxddi,nlvddi
c     +				,nsect,kfluxm,idtflx,nlmax
c     +				,kflux
c     +				,nlayers
c     +				)
c	subroutine wfflx	(iunit,nvers
c     +				,nsect,kfluxm,idtflx,nlmax
c     +				,kflux
c     +				,nlayers
c     +				)
c	subroutine rdflx(iunit,it,nlvddi,nsect,nlayers,fluxes,ierr)
c	subroutine wrflx(iunit,it,nlvddi,nsect,nlayers,fluxes)
c
c revision log :
c
c 18.10.2011	ggu	created from subnos.f
c 10.10.2015	ggu	new routines for new flx framework
c
c************************************************************

	subroutine iniflx

c sets up initial common block - internal routine

	implicit none

c parameters
	integer ftype,maxvers
	parameter(ftype=537,maxvers=5)
c common
	integer mtype,maxver,vers
	common /flxcom/ mtype,maxver,vers
	save /flxcom/
c local
	integer i,n
c save
	logical, save :: binit = .false.

	if( binit ) return

	binit = .true.

	mtype = ftype
	maxver = maxvers
	vers = 0

	end

c************************************************************
c************************************************************
c************************************************************
c************************************************************
c************************************************************

	subroutine infoflx	(iunit,nvers
     +				,nsect,kfluxm,idtflx,nlmax
     +				)

c reads first info record of FLX file
c
c nvers		on entry maximal version that can be read
c		-> must be an input, used to check the corectness
c		.. of the call parameters
c		on return actual version read

	implicit none

c arguments
	integer iunit,nvers
	integer nsect,kfluxm,idtflx,nlmax
c common
	integer mtype,maxver,vers
	common /flxcom/ mtype,maxver,vers
c local
	integer ntype,irec,i

c initialize

	call iniflx

c control newest version number for call

	if( maxver .ne. nvers ) goto 95

c rewind file

	rewind(iunit,err=96)

c first record - find out what version

	irec = 1
	read(iunit,end=91,err=99) ntype,nvers

c control version number and type of file

	if( ntype .ne. mtype ) goto 97
	if( nvers .le. 0 .or. nvers .gt. maxver ) goto 98
	vers = nvers

c next records

	irec = 2
	if( nvers .le. 3 ) then
	   read(iunit,err=99)	 nsect,kfluxm,idtflx
	   nlmax = 0
	else if( nvers .ge. 4 ) then
	   read(iunit,err=99)	 nsect,kfluxm,idtflx,nlmax
	else
	   stop 'error stop rfflx: internal error (1)'
	end if

	return
   99	continue
	write(6,*) 'rfflx: Error encountered while'
	write(6,*) 'reading record number ',irec
	write(6,*) 'of FLX file header'
	write(6,*) 'nvers = ',nvers
	stop 'error stop rfflx: error 99'
   98	continue
	write(6,*) 'rfflx: Version not recognized : ',nvers
	stop 'error stop rfflx: error 98'
   97	continue
	write(6,*) 'rfflx: Wrong type of file : ',ntype
	write(6,*) 'Expected ',mtype
	stop 'error stop rfflx: error 97'
   96	continue
	write(6,*) 'rfflx: Cannot rewind file for unit : ',iunit
	stop 'error stop rfflx: error 96'
   95	continue
	write(6,*) 'rfflx: Old function call ',nvers
	write(6,*) 'nvers = ',nvers,'   maxver = ',maxver
	write(6,*) 'Please adjust call to rfflx and recompile'
	stop 'error stop rfflx: error 95'
   91	continue
	write(6,*) 'rfflx: File is empty'
	stop 'error stop rfflx: error 91'
	end

c************************************************************
c************************************************************
c************************************************************

	function check_flx_file(file)

	implicit none

        logical check_flx_file
        character*(*) file

        integer nb,nvers,ierr
	integer nsect,kfluxm,nlmax
        integer ifileo

        check_flx_file = .false.

        nb = ifileo(0,file,'unform','old')
        if( nb .le. 0 ) return
	call flx_check_header(nb,nvers,nsect,kfluxm,nlmax,ierr)
        close(nb)

        check_flx_file = ierr == 0

	end

c*********************************************************

	subroutine flx_is_flx_file(iunit,nvers)

	implicit none

	integer iunit,nvers

	integer ierr
	integer nsect,kfluxm,nlmax

	call flx_check_header(iunit,nvers,nsect,kfluxm,nlmax,ierr)

	if( ierr .ne. 0 ) nvers = 0

	rewind(iunit)

	end

c*********************************************************
c*********************************************************
c*********************************************************

	subroutine flx_peek_header(iunit,nvers,nsect,kfluxm,nlmax)

	implicit none

	integer iunit,nvers
	integer nsect,kfluxm,nlmax

	integer ierr

	call flx_check_header(iunit,nvers,nsect,kfluxm,nlmax,ierr)
	if( ierr .ne. 0 ) then
	  stop 'error stop flx_peek_header: error reading header'
	end if

	rewind(iunit)

	end

c*********************************************************

	subroutine flx_peek_record(iunit,nvers,it,ierr)

	implicit none

	integer iunit,nvers,it,ierr

	integer ios

	read(iunit,iostat=ios)  it
	ierr = ios

	backspace(iunit)

	end

c*********************************************************

	subroutine flx_check_header(iunit,nvers,nsect,kfluxm,nlmax,ierr)

c checks version of flx file and returns number of points

	implicit none

	integer iunit,nvers
	integer nsect,kfluxm,nlmax
	integer ierr

	integer mtype,maxver,vers
	common /flxcom/ mtype,maxver,vers

	integer kaux,ios,i,it,j
	integer iaux,idtflx,ntype
	real haux,tt,xaux
	character*80 title

	call iniflx

	nvers = 0
	nsect = 0
	kfluxm = 0
	nlmax = 0
	ierr = -1

	read(iunit,iostat=ierr) ntype,nvers
	if( ierr .ne. 0 ) return
	!write(6,*) ierr,ntype,nvers

	ierr = 97
	if( ntype .ne. mtype ) return
	ierr = 99
	if( nvers .le. 0 .or. nvers .gt. maxver ) return
	vers = nvers
	!write(6,*) ierr,ntype,nvers

c next records

	if( nvers .le. 3 ) then
	   read(iunit,iostat=ierr)	 nsect,kfluxm,idtflx
	   nlmax = 0
	else if( nvers .ge. 4 ) then
	   read(iunit,iostat=ierr)	 nsect,kfluxm,idtflx,nlmax
	else
	   stop 'error stop flx_check_header: internal error (1)'
	end if
	!write(6,*) ierr,nsect,kfluxm,idtflx,nlmax
	if( ierr .ne. 0 ) return

	read(iunit,iostat=ierr)	(iaux,i=1,kfluxm)
	if( ierr .ne. 0 ) return

	if( nvers .le. 2 ) then
	else if( nvers .eq. 3 ) then
	  read(iunit,iostat=ierr)	nlmax,(iaux,i=1,nsect)
	else if( nvers .ge. 4 ) then
	  read(iunit,iostat=ierr)	(iaux,i=1,nsect)
	end if
	if( ierr .ne. 0 ) return

	read(iunit,iostat=ierr)  it
	if( ierr .ne. 0 ) return

	ierr = 0

	end

c*********************************************************
c*********************************************************
c*********************************************************

	subroutine flx_read_header(iunit,nscddi,nfxddi,nlvddi
     +				,nvers
     +				,nsect,kfluxm,idtflx,nlmax
     +                          ,kflux,nlayers
     +				)

	implicit none

	integer iunit
	integer nscddi,nfxddi,nlvddi
	integer nvers
	integer nsect,kfluxm,idtflx,nlmax
	integer kflux(kfluxm)
	integer nlayers(nsect)

	call iniflx

	call rfflx		(iunit,nvers
     +				,nscddi,nfxddi,nlvddi
     +				,nsect,kfluxm,idtflx,nlmax
     +				,kflux
     +				,nlayers
     +				)

	end

c*********************************************************

	subroutine flx_read_record(iunit,nvers,it
     +			,nlvddi,nsect,ivar
     +			,nlayers,fluxes,ierr)

	implicit none

	integer iunit,nvers,it
	integer nlvddi,nsect,ivar
	integer nlayers(nsect)
	real fluxes(0:nlvddi,3,nsect)
	integer ierr

	call rdflx(iunit,it,nlvddi,nsect,ivar,nlayers,fluxes,ierr)

	end

c*********************************************************

	subroutine flx_write_header(iunit
     +				,nvers
     +				,nsect,kfluxm,idtflx,nlmax
     +                          ,kflux,nlayers
     +				)

	implicit none

	integer iunit
	integer nvers
	integer nsect,kfluxm,idtflx,nlmax
	integer kflux(kfluxm)
	integer nlayers(nsect)

	call iniflx

	call wfflx		(iunit,nvers
     +				,nsect,kfluxm,idtflx,nlmax
     +				,kflux
     +				,nlayers
     +				)

	end

c*********************************************************

	subroutine flx_write_record(iunit,nvers,it
     +				,nlvddi,nsect,ivar
     +				,nlayers,fluxes)

	implicit none

	integer iunit,nvers,it
	integer nlvddi,nsect,ivar
	integer nlayers(nsect)
	real fluxes(0:nlvddi,3,nsect)

	call wrflx(iunit,it,nlvddi,nsect,ivar,nlayers,fluxes)

	end

c*********************************************************
c*********************************************************
c*********************************************************

c************************************************************
c************************************************************
c************************************************************

	subroutine rfflx	(iunit,nvers
     +				,nscddi,nfxddi,nlvddi
     +				,nsect,kfluxm,idtflx,nlmax
     +				,kflux
     +				,nlayers
     +				)

c reads first record of FLX file
c
c nvers		on entry maximal version that can be read
c		-> must be an input, used to check the corectness
c		.. of the call parameters
c		on return actual version read

	implicit none

c arguments
	integer iunit,nvers
	integer nscddi,nfxddi,nlvddi
	integer nsect,kfluxm,idtflx,nlmax
	integer kflux(kfluxm)
	integer nlayers(nsect)
c common
	integer mtype,maxver,vers
	common /flxcom/ mtype,maxver,vers
c local
	integer ntype,irec,i

c initialize

	call iniflx

c control newest version number for call

	if( maxver .ne. nvers ) goto 95

c rewind file

	rewind(iunit,err=96)

c first record - find out what version

	irec = 1
	read(iunit,end=91,err=99) ntype,nvers

c control version number and type of file

	if( ntype .ne. mtype ) goto 97
	if( nvers .le. 0 .or. nvers .gt. maxver ) goto 98
	vers = nvers

c next records

	irec = 2
	if( nvers .le. 3 ) then
	   read(iunit,err=99)	 nsect,kfluxm,idtflx
	   nlmax = 0
	else if( nvers .ge. 4 ) then
	   read(iunit,err=99)	 nsect,kfluxm,idtflx,nlmax
	else
	   stop 'error stop rfflx: internal error (1)'
	end if

	if( nsect .gt. nscddi ) then
	  write(6,*) 'nscddi,nsect: ',nscddi,nsect
	  stop 'error stop rfflx: dimension nscddi'
	else if( kfluxm .gt. nfxddi ) then
	  write(6,*) 'nfxddi,kfluxm: ',nfxddi,kfluxm
	  stop 'error stop rfflx: dimension nfxddi'
	end if

	irec = 3
	read(iunit,err=99)	(kflux(i),i=1,kfluxm)

	irec = 4
	if( nvers .le. 2 ) then
	  do i=1,nsect
	    nlayers(i) = 1
	  end do
	else if( nvers .eq. 3 ) then
	  read(iunit,err=99)	nlmax,(nlayers(i),i=1,nsect)
	else if( nvers .ge. 4 ) then
	  read(iunit,err=99)	(nlayers(i),i=1,nsect)
	end if

	if( nlmax .gt. nlvddi ) then
	  write(6,*) 'nlvddi,nlmax: ',nlvddi,nlmax
	  stop 'error stop rfflx: dimension nlvddi'
	end if

	return
   99	continue
	write(6,*) 'rfflx: Error encountered while'
	write(6,*) 'reading record number ',irec
	write(6,*) 'of FLX file header'
	write(6,*) 'nvers = ',nvers
	stop 'error stop rfflx: error 99'
   98	continue
	write(6,*) 'rfflx: Version not recognized : ',nvers
	stop 'error stop rfflx: error 98'
   97	continue
	write(6,*) 'rfflx: Wrong type of file : ',ntype
	write(6,*) 'Expected ',mtype
	stop 'error stop rfflx: error 97'
   96	continue
	write(6,*) 'rfflx: Cannot rewind file for unit : ',iunit
	stop 'error stop rfflx: error 96'
   95	continue
	write(6,*) 'rfflx: Old function call ',nvers
	write(6,*) 'nvers = ',nvers,'   maxver = ',maxver
	write(6,*) 'Please adjust call to rfflx and recompile'
	stop 'error stop rfflx: error 95'
   91	continue
	write(6,*) 'rfflx: File is empty'
	stop 'error stop rfflx: error 91'
	end

c********************************************************************

	subroutine wfflx	(iunit,nvers
     +				,nsect,kfluxm,idtflx,nlmax
     +				,kflux
     +				,nlayers
     +				)

c writes first record of FLX file
c
c nvers		on entry maximal version
c		-> must be an input, used to check the corectness
c		.. of the call parameters

	implicit none

c arguments
	integer iunit,nvers
	integer nsect,kfluxm,idtflx,nlmax
	integer kflux(kfluxm)
	integer nlayers(nsect)
c common
	integer mtype,maxver,vers
	common /flxcom/ mtype,maxver,vers
c local
	integer i

	call iniflx

c control newest version number for call

	if( nvers.ne.maxver ) goto 95

	rewind(iunit)

	write(iunit)		mtype,maxver
	write(iunit)		nsect,kfluxm,idtflx,nlmax
	write(iunit)		(kflux(i),i=1,kfluxm)
	write(iunit)		(nlayers(i),i=1,nsect)

	return
   95	continue
	write(6,*) 'wfflx: old function call'
	write(6,*) 'nvers = ',nvers,'   maxver = ',maxver
	write(6,*) 'Please adjust call to wfflx and recompile'
	stop 'error stop wfflx: old version'
	end

c************************************************************

	subroutine rdflx(iunit,it,nlvddi,nsect,ivar,nlayers,fluxes,ierr)

c reads data record of FLX file

	implicit none

c arguments
	integer iunit,it
	integer nlvddi
	integer nsect
	integer ivar
	integer nlayers(nsect)
	real fluxes(0:nlvddi,3,nsect)
	integer ierr
c common
	integer mtype,maxver,vers
	common /flxcom/ mtype,maxver,vers
c local
	integer i,l,j
	integer nvers

	nvers = vers
	ivar = 0

	if( nvers .le. 2 ) then
          read(iunit,end=88,err=99) it,nsect
     +                  ,(fluxes(0,1,i),i=1,nsect)
     +                  ,(fluxes(0,2,i),fluxes(0,3,i),i=1,nsect)
	else if( nvers .le. 4 ) then
          read(iunit,end=88,err=99) it,nsect
     +                  ,(nlayers(i)
     +                  ,((fluxes(l,j,i),l=0,nlayers(i)),j=1,3)
     +                  ,i=1,nsect)
	else if( nvers .ge. 5 ) then
          read(iunit,end=88,err=99) it,nsect,ivar
     +                  ,(nlayers(i)
     +                  ,((fluxes(l,j,i),l=0,nlayers(i)),j=1,3)
     +                  ,i=1,nsect)
	else
	   stop 'error stop rdflx: internal error (1)'
	end if

	ierr=0

	return
   88	continue
	ierr=-1
	return
   99	continue
	write(6,*) 'rdflx: Error while reading'
	write(6,*) 'data record of FLX file'
	write(6,*) 'it = ',it,'  ivar = ',ivar
	ierr=99
	return
	end

c************************************************************

	subroutine wrflx(iunit,it,nlvddi,nsect,ivar,nlayers,fluxes)

c writes data record of FLX file

	implicit none

c arguments
	integer iunit,it
	integer nlvddi
	integer nsect
	integer ivar
	integer nlayers(nsect)
	real fluxes(0:nlvddi,3,nsect)
c common
	integer mtype,maxver,vers
	common /flxcom/ mtype,maxver,vers
c local
	integer i,l,j

        write(iunit) it,nsect,ivar
     +                  ,(nlayers(i)
     +                  ,((fluxes(l,j,i),l=0,nlayers(i)),j=1,3)
     +                  ,i=1,nsect)

	return
	end

c************************************************************

