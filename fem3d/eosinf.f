c
c $Id: eosinf.f,v 1.8 2008-11-20 10:51:34 georg Exp $
c
c info on eos files
c
c revision log :
c
c 18.11.1998    ggu     check dimensions with dimnos
c 06.04.1999    ggu     some cosmetic changes
c 03.12.2001    ggu     some extra output -> place of min/max
c 09.12.2003    ggu     check for NaN introduced
c 07.03.2007    ggu     easier call
c 08.11.2008    ggu     do not compute min/max in non-existing layers
c 07.12.2010    ggu     write statistics on depth distribution (depth_stats)
c 23.01.2011    ggu     write statistics on depth distribution (depth_stats)
c
c**************************************************************

	program eosinf

c reads nos file

	implicit none

	include 'param.h'

	real cv(nkndim)
	real cv3(nlvdim,nkndim)

	integer ilhkv(nkndim)
	integer ilhv(neldim)
	real hlv(nlvdim)
	real hev(neldim)

	integer nread,nin
	integer nvers
	integer nkn,nel,nlv,nvar
	integer ierr
	integer it,ivar
	integer l,k,ie
	character*80 title
	real rnull
	real cmin,cmax

	integer iapini
	integer ifem_open_file

c--------------------------------------------------------------

	nread=0
	rnull=0.
	rnull=-1.

c--------------------------------------------------------------
c open basin and simulation
c--------------------------------------------------------------

	if(iapini(2,0,0,0).eq.0) then
		stop 'error stop : iapini'
	end if

	nin = ifem_open_file('.eos','old')

        nvers=3
	call rheos(nin,nvers,nkndim,neldim,nlvdim,nkn,nel,nlv,nvar
     +				,ilhv,hlv,hev,title)

	call depth_stats(nel,ilhv)

c--------------------------------------------------------------
c loop on data
c--------------------------------------------------------------

	do while(.true.)

	   call rdeos(nin,it,ivar,nlvdim,ilhv,cv3,ierr)

           if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
           if(ierr.ne.0) goto 100

	   nread=nread+1
	   write(6,*) 'time : ',it,'   ivar : ',ivar

	   do l=1,nlv
	     do ie=1,nel
	       cv(ie)=cv3(l,ie)
	       if( l .gt. ilhv(ie) ) cv(ie) = rnull
	     end do
	     call mimar(cv,nel,cmin,cmax,rnull)
             call check1Dr(nel,cv,0.,-1.,"NaN check","cv")
	     write(6,*) 'level,cmin,cmax : ',l,cmin,cmax
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

	subroutine depth_stats(nkn,ilhkv)

c	computes statistics on levels

	implicit none

	include 'param.h'

	integer nkn
	integer ilhkv(1)

	integer count(nlvdim)
	integer ccount(nlvdim)

	integer nlv,lmax,l,k,nc,ll

	nlv = 0
	do l=1,nlvdim
	  count(l) = 0
	  ccount(l) = 0
	end do

	do k=1,nkn
	  lmax = ilhkv(k)
		if( lmax .gt. nlvdim ) stop 'error stop depth_stats: lmax'
		count(lmax) = count(lmax) + 1
		nlv = max(nlv,lmax)
	end do

	do l=nlv,1,-1
		nc = count(l)
	  do ll=1,l
			ccount(ll) = ccount(ll) + nc
		end do
	end do

	nc = 0
	write(6,*) 'statistics for depth: ',nlv
	do l=1,nlv
	  write(6,*) l,count(l),ccount(l)
		nc = nc + count(l)
	end do
	write(6,*) 'total count: ',nc

	end

c***************************************************************

