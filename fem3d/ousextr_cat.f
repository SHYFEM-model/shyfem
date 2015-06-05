c
c $Id: ousextr_cat.f,v 1.3 2009-09-14 08:20:58 georg Exp $
c
c concatenates OUS files
c
c revision log :
c
c 02.09.2003	ggu	adapted to new OUS format
c 24.01.2005	ggu	computes maximum velocities for 3D (only first level)
c 23.03.2010	ggu	extracts reocrds
c 26.03.2010	ggu	bug fix: set nkn and nel
c 03.06.2011    ggu     routine adjourned
c 08.06.2011    ggu     routine rewritten
c 31.01.2012    ggu     choice of concatenation mode (mode, itc)
c 29.04.2015    ggu     mode 1 implemented
c
c***************************************************************

	program ousextr_cat

c concatenates two ous files into one

	implicit none

        include 'param.h'

        integer nrdim
        parameter ( nrdim = 2000 )

        integer irec(nrdim)

	character*80 title

	include 'basin.h'

	include 'nlevel.h'
	include 'levels.h'

	include 'depth.h'

	include 'hydro.h'
	include 'hydro_vel.h'

	real ut2v(neldim)
	real vt2v(neldim)
	real u2v(neldim)
	real v2v(neldim)

	integer ilnv(nlvdim,nkndim)

	character*80 name,file
	character*80 file1,file2
	logical ball,bwrite
        integer nvers,nin
        integer itanf,itend,idt,idtous
	integer it,ie,i
	integer itfirst,itsecond
        integer ierr,nread,nextr,nb
	integer nread1,nread2
        integer nknous,nelous,nlvous
        real href,hzoff,hlvmin
	real zmin,zmax
	real umin,umax
	real vmin,vmax
	real xe,ye
	integer k,ke,ivar
	integer lmax,l
	integer mode,itc

	integer iapini,ideffi,ifileo
	logical berror

c-------------------------------------------------------------------
c initialize params
c-------------------------------------------------------------------

        itfirst = -1
        itsecond = -1
	nread=0
	nread1=0
	nread2=0
	nextr=0

c-------------------------------------------------------------------
c get mode of operation
c-------------------------------------------------------------------

	mode = 0
	itc = 0
	write(6,*) 'operation mode:'
	write(6,*) '  0   all of first file, rest of second file'
	write(6,*) '  1   first file until start of second file'
	write(6,*) '  2   first file up to itc (inclusive), then second'
	write(6,*) '  in case of mode 2 itc value is also requested'
	write(6,*) 'Enter mode:'
	read(5,'(i10)') mode
	if( mode .lt. 0 .or. mode .gt. 2 ) stop 'error stop: mode'
	if( mode .eq. 2 ) then
	  write(6,*) 'Enter time up to when to read first file:'
	  read(5,'(i10)') itc
	end if
        write(6,*) 'mode,itc: ',mode,itc

        write(6,*) 'Enter first file:'
        read(5,'(a)') file1
        write(6,*) file1

        write(6,*) 'Enter second file:'
        read(5,'(a)') file2
        write(6,*) file2

	if( mode .eq. 1 ) then
          call ous_get_it_start(file2,itc)
          if( itc .eq. -1 ) then
            write(6,*) 'cannot read starting time of second file'
            stop 'error stop: starting time'
          end if
          itc = itc - 1         !last time for first file
	end if

c--------------------------------------------------------------------
c open first OUS file and read header
c--------------------------------------------------------------------

	call open_ous_file(file1,'old',nin)

	nvers=1
        call rfous(nin
     +			,nvers
     +			,nknous,nelous,nlvous
     +			,href,hzoff
     +			,title
     +			,ierr)
        if(ierr.ne.0) goto 100

        write(6,*)
        write(6,*)   title
        write(6,*)
        write(6,*) ' nvers        : ',nvers
        write(6,*) ' href,hzoff   : ',href,hzoff
        write(6,*) ' nkn,nel      : ',nknous,nelous
        write(6,*) ' nlv          : ',nlvous
        write(6,*)

	nkn=nknous
	nel=nelous
	nlv=nlvous
	call dimous(nin,nkndim,neldim,nlvdim)

	call rsous(nin,ilhv,hlv,hev,ierr)
        if(ierr.ne.0) goto 100

        write(6,*) 'Available levels: ',nlv
        write(6,*) (hlv(l),l=1,nlv)

c-------------------------------------------------------------------
c open OUS output file
c-------------------------------------------------------------------

        call mkname(' ','ous_cat','.ous',file)
        write(6,*) 'writing file ',file(1:50)
        nb = ifileo(55,file,'unform','new')
        if( nb .le. 0 ) goto 98
	call wfous(nb,1,nkn,nel,nlv,href,hzoff,title,ierr)
        if( ierr .ne. 0 ) goto 99
	call wsous(nb,ilhv,hlv,hev,ierr)
        if( ierr .ne. 0 ) goto 99

c-------------------------------------------------------------------
c loop on input records
c-------------------------------------------------------------------

  300   continue

        call rdous(nin,it,nlvdim,ilhv,znv,zenv,utlnv,vtlnv,ierr)

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 100

	
	if( mode .eq. 2 .and. it .gt. itc ) goto 100
	itfirst = it
	nread=nread+1
	nread1=nread1+1
	write(6,*) 'time : ',nread,it

        call wrous(nb,it,nlvdim,ilhv,znv,zenv,utlnv,vtlnv,ierr)
        if( ierr .ne. 0 ) goto 99

	goto 300
  100	continue

	close(nin)
	call delous(nin)

c--------------------------------------------------------------------
c open second OUS file and read header
c--------------------------------------------------------------------

	call open_ous_file(file2,'old',nin)

	nvers=1
        call rfous(nin
     +			,nvers
     +			,nknous,nelous,nlvous
     +			,href,hzoff
     +			,title
     +			,ierr)
        if(ierr.ne.0) goto 100

        write(6,*)
        write(6,*)   title
        write(6,*)
        write(6,*) ' nvers        : ',nvers
        write(6,*) ' href,hzoff   : ',href,hzoff
        write(6,*) ' nkn,nel      : ',nknous,nelous
        write(6,*) ' nlv          : ',nlvous
        write(6,*)

	nkn=nknous
	nel=nelous
	nlv=nlvous
	call dimous(nin,nkndim,neldim,nlvdim)

	call rsous(nin,ilhv,hlv,hev,ierr)
        if(ierr.ne.0) goto 100

        write(6,*) 'Available levels: ',nlv
        write(6,*) (hlv(l),l=1,nlv)

c-------------------------------------------------------------------
c loop on input records
c-------------------------------------------------------------------

  301   continue

        call rdous(nin,it,nlvdim,ilhv,znv,zenv,utlnv,vtlnv,ierr)

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 101
	if( it .le. itfirst ) goto 301

	if( itsecond .eq. -1 ) itsecond = it
	nread=nread+1
	nread2=nread2+1
	write(6,*) 'time : ',nread,it

        call wrous(nb,it,nlvdim,ilhv,znv,zenv,utlnv,vtlnv,ierr)
        if( ierr .ne. 0 ) goto 99

	goto 301
  101	continue

	close(nin)

c-------------------------------------------------------------------
c end of loop
c-------------------------------------------------------------------

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)
        write(6,*) 'nread1/2: ',nread1,nread2
        write(6,*) 'it: ',itfirst,itsecond
        write(6,*)

        if( nextr .le. 0 ) stop 'no file written'

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

	stop
   98   continue
        write(6,*) 'error opening output file'
        stop 'error stop ousextr_records'
   99   continue
        write(6,*) 'error writing file'
        stop 'error stop ousextr_records'
	end

c******************************************************************

