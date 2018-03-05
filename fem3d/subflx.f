c
c $Id: subflx.f,v 1.15 2007-03-20 13:14:53 georg Exp $
c
c utility routines to read/write FLX file - file type 537
c
c contents :
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
c 03.03.2018	ggu	determine nvar for versions < 6
c
c************************************************************

!==================================================================
        module flxfile
!==================================================================

        implicit none

        integer, parameter :: flx_type = 537
        integer, parameter :: flx_maxvers = 6
        integer, save :: flx_vers = 0		!to be eliminated later

	integer, save :: h1recs(flx_maxvers) = (/3,3,3,2,2,2/)
	integer, save :: h2recs(flx_maxvers) = (/1,1,2,2,2,5/)

!==================================================================
        contains
!==================================================================

!==================================================================
        end module flxfile
!==================================================================

c************************************************************
c************************************************************
c************************************************************

	function check_flx_file(file)

	implicit none

        logical check_flx_file
        character*(*) file

        integer nb,nvers,ierr
	integer nsect,kfluxm,idtflx,nlmax,nvar
        integer ifileo

        check_flx_file = .false.

        nb = ifileo(0,file,'unform','old')
        if( nb .le. 0 ) return
	call flx_check_header(nb,nvers,nsect,kfluxm,idtflx
     +					,nlmax,nvar,ierr)
        close(nb)

        check_flx_file = ( ierr == 0 )

	end

c*********************************************************

	subroutine flx_is_flx_file(iunit,nvers)

	implicit none

	integer iunit,nvers

	integer ierr
	integer nsect,kfluxm,idtflx,nlmax,nvar

	call flx_check_header(iunit,nvers,nsect,kfluxm,idtflx
     +					,nlmax,nvar,ierr)

	if( ierr .ne. 0 ) nvers = 0

	end

c*********************************************************
c*********************************************************
c*********************************************************

	subroutine flx_peek_header(iunit,nvers,nsect,kfluxm,idtflx
     +					,nlmax,nvar,ierr)

	implicit none

	integer iunit,nvers
	integer nsect,kfluxm,idtflx,nlmax,nvar
	integer ierr

	call flx_check_header(iunit,nvers,nsect,kfluxm,idtflx
     +					,nlmax,nvar,ierr)

	end

c*********************************************************

	subroutine flx_peek_record(iunit,nvers,atime,ivar,ierr)

	implicit none

	integer iunit,nvers,ivar,ierr
	double precision atime

	integer nsect,it

	atime = 0.
	it = 0
	ivar = 0

	if( nvers .le. 4 ) then
          read(iunit,iostat=ierr) it,nsect
	else if( nvers .eq. 5 ) then
          read(iunit,iostat=ierr) it,nsect,ivar
	else if( nvers .eq. 6 ) then
          read(iunit,iostat=ierr) atime,nsect,ivar
	else
	   stop 'error stop flx_peek_record: internal error (1)'
	end if

	if( nvers .le. 5 ) atime = it

	backspace(iunit)

	end

c*********************************************************
c*********************************************************
c*********************************************************

	subroutine flx_skip_header(iunit,nvers)

	use flxfile

	implicit none

	integer iunit,nvers

	integer iskip,i

	iskip = h1recs(nvers) +  h2recs(nvers)
	write(6,*) 'iskip = ',iskip
	do i=1,iskip
	  read(iunit)		!empty read - must succeed
	end do
	
	end

c*********************************************************

	subroutine flx_skip_record(iunit,nvers,atime,ivar)

	use flxfile

	implicit none

	integer iunit,nvers
	double precision atime
	integer ivar

	integer ierr

	call flx_peek_record(iunit,nvers,atime,ivar,ierr)
	if( ierr /= 0 ) stop 'error stop: skipping record'

	read(iunit) !peek_record backspaces, we do an empty read
	
	end
	
c*********************************************************
c*********************************************************
c*********************************************************

	subroutine flx_determine_nvar(iunit,nvers,nvar)

	implicit none

	integer iunit,nvers,nvar

	integer ivar,ivar0
	double precision atime

	rewind(iunit)

	call flx_skip_header(iunit,nvers)

	call flx_skip_record(iunit,nvers,atime,ivar0)
	nvar = 1

	do
	  call flx_skip_record(iunit,nvers,atime,ivar)
	  if( ivar == ivar 0 ) exit
	  nvar = nvar + 1
	end do

	rewind(iunit)

	end

c*********************************************************
c*********************************************************
c*********************************************************

	subroutine flx_check_header(iunit,nvers,nsect,kfluxm,idtflx
     +					,nlmax,nvar,ierr)

c checks version of flx file and returns number of points

	use flxfile

	implicit none

	integer iunit,nvers
	integer nsect,kfluxm,idtflx,nlmax,nvar
	integer ierr

	integer ntype

	nvers = 0
	nsect = 0
	kfluxm = 0
	nlmax = 0
	nvar = 0
	ierr = -1

	read(iunit,iostat=ierr) ntype,nvers
	if( ierr .ne. 0 ) goto 99

	ierr = 97
	if( ntype .ne. flx_type ) goto 99
	ierr = 99
	if( nvers .le. 0 .or. nvers .gt. flx_maxvers ) goto 99

c next records

	if( nvers .le. 3 ) then
	   read(iunit,iostat=ierr)	 nsect,kfluxm,idtflx
	   read(iunit,iostat=ierr)	 nlmax	!extra read to get nlmax
	else if( nvers .ge. 4 .and. nvers .le. 5 ) then
	   read(iunit,iostat=ierr)	 nsect,kfluxm,idtflx,nlmax
	else if( nvers .ge. 6 .and. nvers .le. 6 ) then
	   read(iunit,iostat=ierr)	 nsect,kfluxm,idtflx,nlmax,nvar
	else
	   stop 'error stop flx_check_header: internal error (1)'
	end if

	if( ierr .ne. 0 ) goto 99

	ierr = 0

	if( nvers < 6 ) call flx_determine_nvar(iunit,nvers,nvar)

   99	continue
	rewind(iunit)

	end

c*********************************************************
c*********************************************************
c*********************************************************

	subroutine flx_read_header(iunit,nvers,nsect,kfluxm,idtflx
     +					,nlmax,nvar,ierr)

	use flxfile

	implicit none

	integer iunit,nvers
	integer nsect,kfluxm,idtflx,nlmax,nvar,ierr

	integer iskip,i

	call flx_check_header(iunit,nvers,nsect,kfluxm,idtflx
     +					,nlmax,nvar,ierr)
	if( ierr /= 0 ) return

	! the above call rewinds - we will have to skip the header

	iskip = h1recs(nvers)
	do i=1,iskip
	  read(iunit)		!empty read - must succeed
	end do

	end

c*********************************************************

	subroutine flx_read_header2(iunit,nvers,nsect,kfluxm
     +				,kflux,nlayers
     +				,atime0,title,femver,strings,ierr)

	implicit none

	integer iunit,nvers
	integer nsect,kfluxm
	integer kflux(kfluxm)
	integer nlayers(nsect)
	double precision atime0
	character*80 title,femver
	character*80 strings(nsect)
	integer ierr

	integer irec,nlmax

	nlmax = 0
	atime0 = 0.
	title = ' '
	femver = ' '
	strings = ' '

	irec = 3
	read(iunit,err=99)	kflux

	irec = 4
	if( nvers .le. 2 ) then
	  nlayers = 1
	else if( nvers .eq. 3 ) then
	  read(iunit,err=99)	nlmax,nlayers
	else if( nvers .ge. 4 ) then
	  read(iunit,err=99)	nlayers
	end if

	irec = 5
	if( nvers .ge. 6 ) then
	  read(iunit,err=99)	atime0
	  read(iunit,err=99)    title,femver
	  read(iunit,err=99)    strings
	end if

	return
   99	continue
	ierr = irec
	write(6,*) 'read error in second header of FLX file'
	write(6,*) 'nvers,irec: ',nvers,irec
	return
	end

c*********************************************************

	subroutine flx_read_record(iunit,nvers,atime
     +			,nlvddi,nsect,ivar
     +			,nlayers,fluxes,ierr)

	implicit none

	integer iunit,nvers
	double precision atime
	integer nlvddi,nsect,ivar
	integer nlayers(nsect)
	real fluxes(0:nlvddi,3,nsect)
	integer ierr

	integer it,i,j,l

	if( nvers .le. 2 ) then
          read(iunit,end=88,err=99) it,nsect
     +                  ,(fluxes(0,1,i),i=1,nsect)
     +                  ,(fluxes(0,2,i),fluxes(0,3,i),i=1,nsect)
	else if( nvers .le. 4 ) then
          read(iunit,end=88,err=99) it,nsect
     +                  ,(nlayers(i)
     +                  ,((fluxes(l,j,i),l=0,nlayers(i)),j=1,3)
     +                  ,i=1,nsect)
	else if( nvers .eq. 5 ) then
          read(iunit,end=88,err=99) it,nsect,ivar
     +                  ,(nlayers(i)
     +                  ,((fluxes(l,j,i),l=0,nlayers(i)),j=1,3)
     +                  ,i=1,nsect)
	else if( nvers .eq. 6 ) then
          read(iunit,end=88,err=99) atime,nsect,ivar
     +                  ,(nlayers(i)
     +                  ,((fluxes(l,j,i),l=0,nlayers(i)),j=1,3)
     +                  ,i=1,nsect)
	else
	   stop 'error stop rdflx: internal error (1)'
	end if

	if( nvers .le. 5 ) atime = it

	ierr=0

	return
   88	continue
	ierr=-1
	return
   99	continue
	ierr=1
	return
	end

c*********************************************************

	subroutine flx_write_header(iunit,nvers,nsect,kfluxm,idtflx
     +					,nlmax,nvar,ierr)

	use flxfile

	implicit none

	integer iunit,nvers
	integer nsect,kfluxm,idtflx,nlmax,nvar
	integer ierr

	if( nvers /= 0 .and. nvers /= flx_maxvers ) then
          write(6,*) 'cannot write this version for FLX file: ',nvers
          write(6,*) 'please either use 0 or ',flx_maxvers
          ierr = 999
          return
	end if

	rewind(iunit,iostat=ierr)
	if( ierr /= 0 ) return

	write(iunit,iostat=ierr)	flx_type,flx_maxvers
	if( ierr /= 0 ) return
	write(iunit,iostat=ierr)	nsect,kfluxm,idtflx,nlmax,nvar
	if( ierr /= 0 ) return

	end

c*********************************************************

	subroutine flx_write_header2(iunit,nvers,nsect,kfluxm
     +				,kflux,nlayers
     +				,atime0,title,femver,strings,ierr)

	implicit none

	integer iunit,nvers
	integer nsect,kfluxm
	integer kflux(kfluxm)
	integer nlayers(nsect)
	double precision atime0
	character*80 title,femver
	character*80 strings(nsect)
	integer ierr

	integer irec,nlmax

	write(iunit,iostat=ierr)	kflux
	if( ierr /= 0 ) return
	write(iunit,iostat=ierr)	nlayers
	if( ierr /= 0 ) return
	write(iunit,iostat=ierr)	atime0
	if( ierr /= 0 ) return
	write(iunit,iostat=ierr)	title,femver
	if( ierr /= 0 ) return
	write(iunit,iostat=ierr)	strings
	if( ierr /= 0 ) return

	end

c*********************************************************

	subroutine flx_write_record(iunit,nvers,atime
     +			,nlvddi,nsect,ivar
     +			,nlayers,fluxes,ierr)

	implicit none

	integer iunit,nvers
	double precision atime
	integer nlvddi,nsect,ivar
	integer nlayers(nsect)
	real fluxes(0:nlvddi,3,nsect)
	integer ierr

	integer i,j,l

        write(iunit,iostat=ierr) atime,nsect,ivar
     +                  ,(nlayers(i)
     +                  ,((fluxes(l,j,i),l=0,nlayers(i)),j=1,3)
     +                  ,i=1,nsect)

	end

c*********************************************************

	
c*********************************************************
c*********************************************************
c*********************************************************
c old routines ... can be deleted
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

	use flxfile

	implicit none

c arguments
	integer iunit,nvers
	integer nscddi,nfxddi,nlvddi
	integer nsect,kfluxm,idtflx,nlmax
	integer kflux(kfluxm)
	integer nlayers(nsect)
c local
	integer ntype,irec,i

c initialize

c control newest version number for call

	if( flx_maxvers .ne. nvers ) goto 95

c rewind file

	rewind(iunit,err=96)

c first record - find out what version

	irec = 1
	read(iunit,end=91,err=99) ntype,nvers

c control version number and type of file

	if( ntype .ne. flx_type ) goto 97
	if( nvers .le. 0 .or. nvers .gt. flx_maxvers ) goto 98
	flx_vers = nvers

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
	  nlayers = 1
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
	write(6,*) 'Expected ',flx_type
	stop 'error stop rfflx: error 97'
   96	continue
	write(6,*) 'rfflx: Cannot rewind file for unit : ',iunit
	stop 'error stop rfflx: error 96'
   95	continue
	write(6,*) 'rfflx: Old function call ',nvers
	write(6,*) 'nvers = ',nvers,'   maxver = ',flx_maxvers
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

	use flxfile

	implicit none

c arguments
	integer iunit,nvers
	integer nsect,kfluxm,idtflx,nlmax
	integer kflux(kfluxm)
	integer nlayers(nsect)
c local
	integer i

c control newest version number for call

	if( nvers.ne.flx_maxvers ) goto 95

	rewind(iunit)

	write(iunit)		flx_type,flx_maxvers
	write(iunit)		nsect,kfluxm,idtflx,nlmax
	write(iunit)		(kflux(i),i=1,kfluxm)
	write(iunit)		(nlayers(i),i=1,nsect)

	return
   95	continue
	write(6,*) 'wfflx: old function call'
	write(6,*) 'nvers = ',nvers,'   maxver = ',flx_maxvers
	write(6,*) 'Please adjust call to wfflx and recompile'
	stop 'error stop wfflx: old version'
	end

c************************************************************

	subroutine rdflx(iunit,it,nlvddi,nsect,ivar,nlayers,fluxes,ierr)

c reads data record of FLX file

	use flxfile

	implicit none

c arguments
	integer iunit,it
	integer nlvddi
	integer nsect
	integer ivar
	integer nlayers(nsect)
	real fluxes(0:nlvddi,3,nsect)
	integer ierr
c local
	integer i,l,j
	integer nvers

	nvers = flx_vers
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
c local
	integer i,l,j

        write(iunit) it,nsect,ivar
     +                  ,(nlayers(i)
     +                  ,((fluxes(l,j,i),l=0,nlayers(i)),j=1,3)
     +                  ,i=1,nsect)

	return
	end

c************************************************************

