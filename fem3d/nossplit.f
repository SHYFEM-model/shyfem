c
c $Id: nossplit.f,v 1.3 2009-04-07 10:43:57 georg Exp $
c
c 18.11.1998    ggu     check dimensions with dimnos
c
c**********************************************************

	program nossplit

c splits nos file

	include 'param.h'

	parameter ( ndim = 150 )

	integer iu(ndim)
	character*80 name,file

	character*80 title
	real cv(nkndim)
	real cv3(nlvdim,nkndim)

	integer ilhkv(nkndim)
	real hlv(nlvdim)
	real hev(neldim)

c-----------------------------------------------------------

	do i=1,ndim
	  iu(i) = 0
	end do

	nread=0
	rnull=0.

	if(iapini(2,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

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

c loop on input records

  300   continue

	call rdnos(nin,it,ivar,nlvdim,ilhkv,cv3,ierr)

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 100

	nread=nread+1
	write(6,*) 'time : ',it,ivar

	if( ivar .gt. ndim ) goto 95

	if( iu(ivar) .le. 0 ) then	!open file
	  write(name,'(i4)') ivar
	  call mkname(' ',name,'.nos',file)
	  write(6,*) 'opening file ',file(1:50)
	  nb = ifileo(55,file,'unform','new')
	  if( nb .le. 0 ) goto 98
	  call wfnos(nb,3,nkn,nel,nlv,1,title,ierr)
	  if( ierr .ne. 0 ) goto 99
	  call wsnos(nb,ilhkv,hlv,hev,ierr)
	  if( ierr .ne. 0 ) goto 99
	  iu(ivar) = nb
	end if

	nb = iu(ivar)
	call wrnos(nb,it,ivar,nlvdim,ilhkv,cv3,ierr)
	if( ierr .ne. 0 ) goto 99

	goto 300

  100	continue

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)

	stop
   95	continue
	write(6,*) 'error ndim : ',ivar,ndim
	stop 'error stop nossplit'
   98	continue
	write(6,*) 'error opening file'
	stop 'error stop nossplit'
   99	continue
	write(6,*) 'error writing file'
	stop 'error stop nossplit'
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

