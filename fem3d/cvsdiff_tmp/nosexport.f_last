c
c $Id: nosexport.f,v 1.1 2005/11/03 17:12:43 georg Exp $
c
c creates files for initialization of NOS fields
c
c revision log :
c
c 05.03.2004    ggu     copied from nosextr
c
c****************************************************************

	program nosexport

c exports nos files to ascii file

	implicit none

	include 'param.h'

c--------------------------------------------------
        character*80 descrr
        common /descrr/descrr
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real xgv(nkndim), ygv(nkndim)
        real hm3v(3,neldim)
        integer nen3v(3,neldim)
        integer ipv(nkndim), ipev(neldim)
        integer iarv(neldim)

        common /xgv/xgv, /ygv/ygv
        common /hm3v/hm3v
        common /nen3v/nen3v
        common /ipv/ipv, /ipev/ipev
        common /iarv/iarv
c--------------------------------------------------

	character*80 title
	real cv(nkndim)
	real cv3(nlvdim,nkndim)

	integer ilhkv(nkndim)
	real hlv(nlvdim)
	real hev(neldim)

	logical berror
	integer nin,nvers
	integer nlv,nvar,ierr
	integer k,l,nread,ivar,it
	double precision xgb,ygb,xgb0,ygb0

	integer iapini,ideffi

c--------------------------------------------------
c Gauss-Boaga
c--------------------------------------------------

	xgb0 = 2330000. - 50000.
	ygb0 = 5000000.

c---------------------------------------------------------------
c open simulation and basin
c---------------------------------------------------------------

	if(iapini(3,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

c---------------------------------------------------------------
c open files and read headers
c---------------------------------------------------------------

	nin=ideffi('datdir','runnam','.nos','unform','old')
	if(nin.le.0) goto 100

        nvers=3
	call rfnos(nin,nvers,nkn,nel,nlv,nvar,title,ierr)
        if(ierr.ne.0) goto 100

        write(6,*) 'nvers    : ',nvers
        write(6,*) 'nkn,nel  : ',nkn,nel
        write(6,*) 'nlv,nvar : ',nlv,nvar
        write(6,*) 'title    : ',title

        call dimnos(nin,nkndim,neldim,nlvdim)

	call rsnos(nin,ilhkv,hlv,hev,ierr)
        if(ierr.ne.0) goto 100

	write(6,*) 'Available levels: ',nlv
	write(6,*) (hlv(l),l=1,nlv)

c---------------------------------------------------------------
c what to extract
c---------------------------------------------------------------

        open(1,file='nosexport.txt',status='unknown',form='formatted')
        write(1,*) nkn,nlv
        write(1,*) (hlv(l),l=1,nlv)                     !we need level struct.

c---------------------------------------------------------------
c time loop
c---------------------------------------------------------------

  300   continue

	call rdnos(nin,it,ivar,nlvdim,ilhkv,cv3,ierr)

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 100

	nread=nread+1
	write(6,*) 'time : ',it,ivar

	write(1,*) it,ivar
	do k=1,nkn
	  xgb = xgb0 + xgv(k)
	  ygb = ygb0 + ygv(k)
          write(1,*) k,xgb,ygb,(cv3(l,k),l=1,nlv)
	end do

	goto 300

c---------------------------------------------------------------
c end of time loop
c---------------------------------------------------------------

  100	continue

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)
	write(6,*) 'data written to file 79'
	write(6,*)

        close(1)

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

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

