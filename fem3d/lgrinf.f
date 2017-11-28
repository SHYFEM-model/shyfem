c
c $Id: lgrinf.f,v 1.1 2008-07-16 15:41:39 georg Exp $
c
c revision log :
c
c 18.11.1998    ggu     check dimensions with dimnos
c 06.04.1999    ggu     some cosmetic changes
c 03.12.2001    ggu     some extra output -> place of min/max
c 09.12.2003    ggu     check for NaN introduced
c 07.03.2007    ggu     easier call
c 23.04.2015    ggu     for new version 5
c
c**************************************************************

	program lgrinf

c reads lgr file

	implicit none

	integer nread,nin,i,it,nc
	integer mtype,nvers
	integer nbdy,nn,nout
	integer id,ie,ies,lmax,lb
	real x,y,z,xs,ys
	integer id_special
	character*80 file

	integer iapini
	integer ifem_open_file

c--------------------------------------------------------------

	nread=0
	id_special = 4
	id_special = 1
	id_special = 0

c--------------------------------------------------------------
c open simulation
c--------------------------------------------------------------

	nin = 1
	nc = command_argument_count()
        if( nc == 0 ) stop 'no file given'
        call get_command_argument(1,file)
	open(nin,file=file,status='old',form='unformatted')

	read(nin) mtype,nvers
	write(6,*) 'mtype,nvers: ',mtype,nvers
	if( nvers > 4 ) read(nin) lmax
	if( mtype .ne. 367265 ) stop 'error stop: mtype'

c--------------------------------------------------------------
c loop on data
c--------------------------------------------------------------

	do while(.true.)

	   read(nin,end=100) it,nbdy,nn,nout
	   write(6,*) it,nbdy,nn,nout

	   nread = nread + 1

	   do i=1,nn
	     if( nvers < 5 ) then
	       read(nin) id,x,y,z,ie,xs,ys,ies
	     else
	       read(nin) id,x,y,z,ie,lb,xs,ys,ies
	     end if
	     !if( ie .lt. 0 ) write(6,*) -ie,x,y
	     if( id_special < 0 .or. id == id_special ) then
	       write(6,*) id,x,y,z,lb
	     end if
	   end do

	end do	!do while

c--------------------------------------------------------------
c end of loop on data
c--------------------------------------------------------------

  100	continue

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	end

c***************************************************************

