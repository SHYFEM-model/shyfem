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
c local
	integer i,n
c save
	logical binit
	save binit
	save /flxcom/
c data
	data binit /.false./

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
	integer kflux(1)
	integer nlayers(1)
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
	integer kflux(1)
	integer nlayers(1)
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
	integer nlayers(1)
	real fluxes(0:nlvddi,3,1)
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
	integer nlayers(1)
	real fluxes(0:nlvddi,3,1)
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

