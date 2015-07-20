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

c reads nos file

	implicit none

	include 'param.h'

	integer nread,nin,i,it
	integer mtype,nvers
	integer nbdy,nn,nout
	integer id,ie,ies,lmax,lb
	real x,y,z,xs,ys
	integer id_special

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

	if(iapini(2,0,0,0).eq.0) then
		stop 'error stop : iapini'
	end if

	nin = ifem_open_file('.lgr','old')

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

        subroutine mimar(xx,n,xmin,xmax,rnull)

c computes min/max of vector
c
c xx            vector
c n             dimension of vector
c xmin,xmax     min/max value in vector
c rnull		invalid value

        implicit none

        integer n,i,nmin
        real xx(n)
        real xmin,xmax,x,rnull

	do i=1,n
	  if(xx(i).ne.rnull) goto 1
	end do
    1	continue

	if(i.le.n) then
	  xmax=xx(i)
	  xmin=xx(i)
	else
	  xmax=rnull
	  xmin=rnull
	end if

	nmin=i+1

        do i=nmin,n
          x=xx(i)
	  if(x.ne.rnull) then
            if(x.gt.xmax) xmax=x
            if(x.lt.xmin) xmin=x
	  end if
        end do

        end

c***************************************************************

        subroutine mimar_s(xx,nlvddi,n,xmin,xmax,rnull)

c computes min/max of vector
c
c xx            vector
c n             dimension of vector
c xmin,xmax     min/max value in vector
c rnull		invalid value

        implicit none

        integer n
        integer nlvddi
        real xx(nlvddi,n)
        real xmin,xmax,rnull

        integer k,l
        real x

        do k=1,n
          do l=1,nlvddi
            x=xx(l,k)
	    if(x.ne.rnull) then
              if( x .lt. xmin .or. x .gt. xmax ) then
                write(6,*) l,k,x
              end if
	    end if
          end do
        end do

        end

c***************************************************************

