
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

c routines for random number generation
c
c contents :
c
c function grand(idum)
c function ran(idum)
c
c function ran0(idum)		!shuffling random results
c function ran1(idum)		!best random number generator
c function ran2(idum)		!good enough random number generator
c function ran3(idum)		!different method (Knuth)
c function ran4(idum)		!using des
c
c revision log :
c
c 19.11.1998    ggu     routines copied from numerical receipies
c 23.11.1998    ggu     test routines and comments
c 30.03.2000    ggu     pause statement eliminated, variables saved (bug)
c 11.12.2007	ggu	new routine ggrand
c
c******************************************************

	function grand(idum)

	grand = ran2(idum)

	end

c******************************************************

	function ggrand(init)

c works for initialization and normal call

	integer icall,idum
	save icall,idum
	data icall,idum / 0 , 99 /

	if( icall .eq. 0 .or. init .lt. 0 ) then
	  icall = 1
	  if( init .ne. 0 ) idum = abs(init)
	end if

	ggrand = ran2(idum)

	end

c******************************************************

	function ran00(idum)

	ran00 = ran1(idum)

	end

c******************************************************

      function ran0(idum)
      dimension v(97)
      save y,v
      data iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        iseed=abs(idum)
        idum=1
        do 11 j=1,97
          dum=ran(iseed)
11      continue
        do 12 j=1,97
          v(j)=ran(iseed)
12      continue
        y=ran(iseed)
      endif
      j=1+int(97.*y)
      if(j.gt.97.or.j.lt.1) stop 'error stop ran0: internal error'
      y=v(j)
      ran0=y
      v(j)=ran(iseed)
      return
      end

c******************************************************

      function ran1(idum)
      dimension r(97)
      parameter (m1=259200,ia1=7141,ic1=54773,rm1=3.8580247e-6)
      parameter (m2=134456,ia2=8121,ic2=28411,rm2=7.4373773e-6)
      parameter (m3=243000,ia3=4561,ic3=51349)
      save ix1,ix2,ix3,r
      data iff /0/
      if (idum.lt.0.or.iff.eq.0) then
        iff=1
        ix1=mod(ic1-idum,m1)
        ix1=mod(ia1*ix1+ic1,m1)
        ix2=mod(ix1,m2)
        ix1=mod(ia1*ix1+ic1,m1)
        ix3=mod(ix1,m3)
        do 11 j=1,97
          ix1=mod(ia1*ix1+ic1,m1)
          ix2=mod(ia2*ix2+ic2,m2)
          r(j)=(float(ix1)+float(ix2)*rm2)*rm1
11      continue
        idum=1
      endif
      ix1=mod(ia1*ix1+ic1,m1)
      ix2=mod(ia2*ix2+ic2,m2)
      ix3=mod(ia3*ix3+ic3,m3)
      j=1+(97*ix3)/m3
      if(j.gt.97.or.j.lt.1) stop 'error stop ran1: internal error'
      ran1=r(j)
      r(j)=(float(ix1)+float(ix2)*rm2)*rm1
      return
      end

c******************************************************

      function ran2(idum)
      parameter (m=714025,ia=1366,ic=150889,rm=1.4005112e-6)
      dimension ir(97)
      save iy,ir
      data iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        idum=mod(ic-idum,m)
        do 11 j=1,97
          idum=mod(ia*idum+ic,m)
          ir(j)=idum
11      continue
        idum=mod(ia*idum+ic,m)
        iy=idum
      endif
      j=1+(97*iy)/m
      if(j.gt.97.or.j.lt.1) stop 'error stop ran2: internal error'
      iy=ir(j)
      ran2=iy*rm
      idum=mod(ia*idum+ic,m)
      ir(j)=idum
      return
      end

c******************************************************

      function ran3(idum)
c         implicit real*4(m)
c         parameter (mbig=4000000.,mseed=1618033.,mz=0.,fac=2.5e-7)
      parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1.e-9)
      dimension ma(55)
      data iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=mseed-iabs(idum)
        mj=mod(mj,mbig)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.mz)mk=mk+mbig
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.mz)ma(i)=ma(i)+mbig
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz)mj=mj+mbig
      ma(inext)=mj
      ran3=mj*fac
      return
      end

c******************************************************

      function ran4(idum)
      parameter (im=11979,ia=430,ic=2531,nacc=24)
      dimension inp(64),jot(64),key(64),pow(65)
      data iff/0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        idum=mod(idum,im)
        pow(1)=0.5
        do 11 j=1,64
          idum=mod(idum*ia+ic,im)
          key(j)=(2*idum)/im
          inp(j)=mod((4*idum)/im,2)
          pow(j+1)=0.5*pow(j)
11      continue
        newkey=1
      endif
      isav=inp(64)
      if(isav.ne.0)then
        inp(4)=1-inp(4)
        inp(3)=1-inp(3)
        inp(1)=1-inp(1)
      endif
      do 12 j=64,2,-1
        inp(j)=inp(j-1)
12    continue
      inp(1)=isav
      jot = 0
      call des(inp,key,newkey,0,jot)
      ran4=0.0
      do 13 j=1,nacc
        if(jot(j).ne.0)ran4=ran4+pow(j)
13    continue
      return
      end

c******************************************************
c******************************************************
c            dummy and test routines
c******************************************************
c******************************************************

	subroutine des(inp,key,newkey,i,jot)
        integer inp(1),jot(1),key(1)
	end

c******************************************************

	subroutine rndtest(n,ndim,icount)

c tests random routine

	implicit none

	integer n
	integer ndim
	integer icount(ndim)

	integer i,j
	real grand

	integer iseed
	save iseed
	data iseed / -346723 /

	do i=1,ndim
	  icount(i) = 0
	end do

	do i=1,n
	  j = 1 + int(ndim*grand(iseed))
	  icount(j) = icount(j) + 1
	end do

	call rndprnt(n,ndim,icount)

	end

c******************************************************

	subroutine rndprnt(n,ndim,icount)

c prints result of random test

	implicit none

	integer n
	integer ndim
	integer icount(ndim)

	integer i
	integer medium,ic,im

	medium = n / ndim

	do i=1,ndim
	  ic = icount(i)
	  im = ic - medium
	  write(6,'(3i10,f10.4)') i,ic,im,abs((100.*im)/medium)
	end do

	end

c******************************************************

	subroutine rnddrv(n,k)

c driver for random test
c
c n	total number of iterations / test
c k	total number of tests

	implicit none

	integer n,k

        integer ndim
        parameter(ndim=20)
        integer icount(ndim)
        integer itot(ndim)

	integer i,l

	do i=1,ndim
	  itot(i) = 0
	end do

	do l=1,k
	  write(6,*) 'Iteration : ',l,' of ',k
	  call rndtest(n,ndim,icount)
	  do i=1,ndim
	    itot(i) = itot(i) + icount(i)
	  end do
	end do

	write(6,*) 'Final result : '
	call rndprnt(n*k,ndim,itot)

	end

c******************************************************
c
c uncomment next lines to run test on random number generator
c
c	program test
c	call rnddrv(1000000,10)
c	end
c
c******************************************************
	

