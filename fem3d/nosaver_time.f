c
c $Id: nosaverf.f,v 1.1 2008-04-11 16:05:37 georg Exp $
c
c revision log :
c
c 18.11.1998    ggu     check dimensions with dimnos
c 07.03.2007    ggu     easier calls
c 10.04.2008    ggu     copied from nosaver -> frequency introduced
c 29.04.2010    ggu     some changes
c
c**********************************************************

	program nosaver_time

c averages records from nos file (with frequence given)
c
c creates 4 files: aver, min, max, sum
c averages with given frequency (number of records to use)
c if frequency is 0, averages the whole file
c the last written record might be averaged over less records

        implicit none

	include 'param.h'

	character*80 title
	real cv3(nlvdim,nkndim)

        double precision accum(nlvdim,nkndim)
	real amin(nlvdim,nkndim)
	real amax(nlvdim,nkndim)
	real asum(nlvdim,nkndim)

	integer ilhkv(nkndim)
	real hlv(nlvdim)
	real hev(neldim)

        integer nread,ivarold,nwrite
	integer nfreq,naccum
        integer l,k,nin
	integer nbaver,nbmin,nbmax,nbsum
        integer nkn,nel,nlv,nvar
        integer it,ivar
        integer ierr
        integer nvers
        real r
	real conz,high

        integer iapini,ideffi,ifileo,ifem_open_file

c-----------------------------------------------------------------
c initialize params
c-----------------------------------------------------------------

	nread = 0
	nwrite = 0
	nfreq = 30
        ivarold = 0
	high = 1.e+30

	naccum = 0
        do l=1,nlvdim
          do k=1,nkndim
            accum(l,k) = 0.
            amin(l,k) = high
            amax(l,k) = -high
          end do
        end do

c-----------------------------------------------------------------
c open input file
c-----------------------------------------------------------------

	if(iapini(2,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

        nin = ifem_open_file('.nos','old')

        nvers=3
        call rhnos(nin,nvers,nkndim,neldim,nlvdim,nkn,nel,nlv,nvar
     +                          ,ilhkv,hlv,hev,title)

c-----------------------------------------------------------------
c get frequency
c-----------------------------------------------------------------

	write(6,*) 'Enter frequency of averaging (0 for whole file): '
	read(5,'(i10)') nfreq
	write(6,*) 'averaging with frequency ',nfreq

c-----------------------------------------------------------------
c open output files
c-----------------------------------------------------------------

	call open_nos_file('aver','new',nbaver)
	call whnos(nbaver,3,nkn,nel,nlv,1,ilhkv,hlv,hev,title)

	call open_nos_file('min','new',nbmin)
	call whnos(nbmin,3,nkn,nel,nlv,1,ilhkv,hlv,hev,title)

	call open_nos_file('max','new',nbmax)
	call whnos(nbmax,3,nkn,nel,nlv,1,ilhkv,hlv,hev,title)

	call open_nos_file('sum','new',nbsum)
	call whnos(nbsum,3,nkn,nel,nlv,1,ilhkv,hlv,hev,title)

c-----------------------------------------------------------------
c loop on input records
c-----------------------------------------------------------------

	do while(.true.)

c	----------------------------------------------------------
c	read next record
c	----------------------------------------------------------

	call rdnos(nin,it,ivar,nlvdim,ilhkv,cv3,ierr)

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 100
        if( ivarold .eq. 0 ) ivarold = ivar
        if( ivar .ne. ivarold ) goto 91

	nread=nread+1
	write(6,*) 'time : ',it,ivar

c	----------------------------------------------------------
c	accumulate results
c	----------------------------------------------------------

	naccum = naccum + 1
        do l=1,nlv
          do k=1,nkn
	    conz = cv3(l,k)
            accum(l,k) = accum(l,k) + conz
	    amin(l,k) = min(amin(l,k),conz)
	    amax(l,k) = max(amax(l,k),conz)
          end do
        end do

c	----------------------------------------------------------
c	if it is time -> write output
c	----------------------------------------------------------

	if( nfreq .gt. 0 .and. mod(nread,nfreq) .eq. 0 ) then

	  nwrite = nwrite + 1
	  write(6,*) 'writing output files: ',nread,naccum,nwrite

          r = 1./naccum
          do l=1,nlv
            do k=1,nkn
              asum(l,k) = accum(l,k)
              cv3(l,k) = accum(l,k) * r
            end do
          end do

	  call wrnos(nbaver,it,ivar,nlvdim,ilhkv,cv3,ierr)	!aver
	  call wrnos(nbmin,it,ivar,nlvdim,ilhkv,amin,ierr)	!min
	  call wrnos(nbmax,it,ivar,nlvdim,ilhkv,amax,ierr)	!max
	  call wrnos(nbsum,it,ivar,nlvdim,ilhkv,asum,ierr)	!sum

	  naccum = 0
          do l=1,nlvdim
            do k=1,nkndim
              accum(l,k) = 0.
              amin(l,k) = high
              amax(l,k) = -high
            end do
          end do

	end if

	end do

c-----------------------------------------------------------------
c end of loop on input records
c-----------------------------------------------------------------

  100	continue

c-----------------------------------------------------------------
c elaborate last records if there are any
c-----------------------------------------------------------------

	if( naccum .gt. 0 ) then

	  nwrite = nwrite + 1
	  write(6,*) 'writing output files: ',nread,naccum,nwrite

          r = 1./naccum
          do l=1,nlv
            do k=1,nkn
              asum(l,k) = accum(l,k)
              cv3(l,k) = accum(l,k) * r
            end do
          end do

	  call wrnos(nbaver,it,ivar,nlvdim,ilhkv,cv3,ierr)	!aver
	  call wrnos(nbmin,it,ivar,nlvdim,ilhkv,amin,ierr)	!min
	  call wrnos(nbmax,it,ivar,nlvdim,ilhkv,amax,ierr)	!max
	  call wrnos(nbsum,it,ivar,nlvdim,ilhkv,asum,ierr)	!sum

	end if

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	write(6,*)
	write(6,*) 'frequency used',nfreq
	write(6,*)
	write(6,*) nread,' records read in total'
	write(6,*) nwrite,' records written in total'
	write(6,*)

	stop
   91	continue
	write(6,*) 'error ivar : ',ivar,ivarold
	write(6,*) 'can handle only one variable type in file'
	write(6,*) 'please use nossplit to extract the variable first'
	stop 'error stop nosaver'
	end

c***************************************************************

