c
c $Id: nosinit.f,v 1.2 2005/11/03 16:59:25 georg Exp $
c
c creates files for initialization of NOS fields
c
c revision log :
c
c 05.03.2004    ggu     copied from nosextr
c 17.06.2005    ggu     write level structure only up to nlv
c
c****************************************************************

	program nosinit

c creates files for initialization of NOS fields

	include 'param.h'

c--------------------------------------------------
        character*80 descrr
        common /descrr/descrr
	include 'nbasin.h'

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

c--------------------------------------------------
	integer nudim
	parameter(nudim=300)
	integer iunit(nudim)
	character*50 file

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

        write(6,*) 'Specify it, ivar, ilev'
        write(6,*) '   it     time of map to extract'
        write(6,*) '   ivar   variable (ivar=0 -> all vars)'
        write(6,*) '   ilev   level (ilev=0 -> all levels)'
        write(6,*) 'Enter it ivar ilev : '
        read(5,*)  itact,ivaract,ilevact
        write(6,*) 'it,ivar,ilev: ',itact,ivaract,ilevact

        if( ilevact .gt. nlv ) then
          write(6,*) nlv,ilevact
          stop 'error stop: ilev must be smaller than nlv'
        end if

        ivar = nvar
        if( ivaract .gt. 0 ) ivar = 1
        ilev = nlv
        if( ilevact .gt. 0 ) ilev = 1

        open(1,file='nosinit.dat',status='unknown',form='unformatted')
        write(1) nkn,ilev,ivar

        if( ilevact .eq. 0 .and. nlv .gt. 1 ) then       !write all levels
          write(1) (hlv(l),l=1,nlv)                     !we need level struct.
        end if

        open(2,file='nosinit.nos',status='unknown',form='unformatted')
        call wfnos(2,3,nkn,nel,nlv,1,title,ierr)
        call wsnos(2,ilhkv,hlv,hev,ierr)

c---------------------------------------------------------------
c time loop
c---------------------------------------------------------------

  300   continue

	call rdnos(nin,it,ivar,nlvdim,ilhkv,cv3,ierr)

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 100

	nread=nread+1
	write(6,*) 'time : ',it,ivar

        if( it .eq. itact ) then
          if( ivaract .eq. 0 .or. ivar .eq. ivaract ) then
            if( ilev .eq. 1 ) then
              write(1) (cv3(ilevact,k),k=1,nkn)
            else
              write(1) ((cv3(l,k),l=1,nlv),k=1,nkn)
            end if
            call wrnos(2,it,ivar,nlvdim,ilhkv,cv3,ierr)
          end if
        end if

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
        close(2)

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

