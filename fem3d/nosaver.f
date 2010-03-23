c
c $Id: nosaver.f,v 1.4 2008-04-11 16:05:37 georg Exp $
c
c 18.11.1998    ggu     check dimensions with dimnos
c 07.03.2007    ggu     easier calls
c 10.04.2008    ggu     new name for file: aver.nos -> aver_amms.nos
c
c**********************************************************

	program nosaver

c averages records from nos file
c
c creates one file containing 4 records: aver, min, max, sum
c averages over whole file

        implicit none

	include 'param.h'

        integer ndim
	parameter ( ndim = 150 )

	integer iu(ndim)
	character*80 name,file

	character*80 title
	real cv(nkndim)
	real cv3(nlvdim,nkndim)

        double precision accum(nlvdim,nkndim)
	real amin(nlvdim,nkndim)
	real amax(nlvdim,nkndim)
	real asum(nlvdim,nkndim)

	integer ilhkv(nkndim)
	real hlv(nlvdim)
	real hev(neldim)

        integer nread,ivarold
        integer l,k,nin,nb
        integer nkn,nel,nlv,nvar
        integer it,ivar
        integer ierr
        integer nvers
        real r,rnull
	real conz,high

        integer iapini,ideffi,ifileo,ifem_open_file

c-----------------------------------------------------------

	nread=0
	rnull=0.
        ivarold = 0
	high = 1.e+30

        do l=1,nlvdim
          do k=1,nkndim
            accum(l,k) = 0.
            amin(l,k) = high
            amax(l,k) = -high
          end do
        end do

c-----------------------------------------------------------------
c open files
c-----------------------------------------------------------------

	if(iapini(2,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

        nin = ifem_open_file('.nos','old')

        nvers=3
        call rhnos(nin,nvers,nkndim,neldim,nlvdim,nkn,nel,nlv,nvar
     +                          ,ilhkv,hlv,hev,title)

c-----------------------------------------------------------------
c loop on input records
c-----------------------------------------------------------------

	do while(.true.)

	call rdnos(nin,it,ivar,nlvdim,ilhkv,cv3,ierr)

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 100
        if( ivarold .eq. 0 ) ivarold = ivar
        if( ivar .ne. ivarold ) goto 91

	nread=nread+1
	write(6,*) 'time : ',it,ivar

        do l=1,nlv
          do k=1,nkn
	    conz = cv3(l,k)
            accum(l,k) = accum(l,k) + conz
	    amin(l,k) = min(amin(l,k),conz)
	    amax(l,k) = max(amax(l,k),conz)
          end do
        end do

	end do

c-----------------------------------------------------------------
c end of loop on input records
c-----------------------------------------------------------------

  100	continue

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)

        if( nread .le. 0 ) stop 'no file written'

        r = 1./nread

        do l=1,nlv
          do k=1,nkn
            asum(l,k) = accum(l,k)
            cv3(l,k) = accum(l,k) * r
          end do
        end do

c-----------------------------------------------------------------
c write averaged file
c-----------------------------------------------------------------

        it = 0
        ivar = ivarold

	call mkname(' ','aver_amms','.nos',file)
	write(6,*) 'writing file ',file(1:50)
	nb = ifileo(55,file,'unform','new')
	if( nb .le. 0 ) goto 98

	call wfnos(nb,3,nkn,nel,nlv,1,title,ierr)
	if( ierr .ne. 0 ) goto 99
	call wsnos(nb,ilhkv,hlv,hev,ierr)
	if( ierr .ne. 0 ) goto 99
	call wrnos(nb,it,ivar,nlvdim,ilhkv,cv3,ierr)	!aver
	if( ierr .ne. 0 ) goto 99
	call wrnos(nb,it,ivar,nlvdim,ilhkv,amin,ierr)	!min
	if( ierr .ne. 0 ) goto 99
	call wrnos(nb,it,ivar,nlvdim,ilhkv,amax,ierr)	!max
	if( ierr .ne. 0 ) goto 99
	call wrnos(nb,it,ivar,nlvdim,ilhkv,asum,ierr)	!sum
	if( ierr .ne. 0 ) goto 99

	stop
   91	continue
	write(6,*) 'error ivar : ',ivar,ivarold
	stop 'error stop nosaver'
   95	continue
	write(6,*) 'error ndim : ',ivar,ndim
	stop 'error stop nosaver'
   98	continue
	write(6,*) 'error opening file'
	stop 'error stop nosaver'
   99	continue
	write(6,*) 'error writing file'
	stop 'error stop nosaver'
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

