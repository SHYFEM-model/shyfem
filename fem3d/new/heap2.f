
!--------------------------------------------------------------------------
!
!    Copyright (C) 2010,2019  Georg Umgiesser
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

c**************************************************************

! revision log :
!
! 23.03.2010	ggu	changed v6.1.1
! 16.02.2019	ggu	changed VERS_7_5_60

c**************************************************************

      SUBROUTINE SORT(N,RA)

c heapsort

      implicit none

      integer n
      real RA(N)
c      integer ra(n)

      integer i,ir,j,l
      real rra
c      integer rra

      if(n.lt.2) return

      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
      GO TO 10
      END

c******************************************************************

      SUBROUTINE iSORT(N,RA,index)

c heapsort

      implicit none

      integer n
      integer RA(N),index(n)

      integer i,ir,j,l,irra
      integer rra

      do i=1,n
	index(i)=i
      end do

      if(n.lt.2) return

      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
	  irra=index(l)
          RRA=RA(irra)
        ELSE
	  irra=index(ir)
          RRA=RA(irra)
          index(ir)=index(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            index(1)=irra
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(index(j)).LT.RA(index(J+1)))J=J+1
          ENDIF
          IF(RRA.LT.RA(index(j)))THEN
            index(i)=index(j)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        index(i)=irra
      GO TO 10
      END

c**********************************************************

      integer function locate(n,ixx,index,ix)

c locates index for ix in array ixx ( ixx(index(i)) is sorted )

      implicit none

      integer n
      integer ixx(n)
      integer index(n)
      integer ix

      integer jl,ju,jm

      jl=1
      ju=n

      do while(ju.gt.jl)
        jm=(ju+jl)/2
        if(ix.gt.ixx(index(jm)))then
          jl=jm+1
        else
          ju=jm
        endif
      end do

      if(ju.le.0) then                    !empty list
        locate=0
      else if(ix.eq.ixx(index(ju))) then  !found
        locate=index(ju)
      else                                !not found
        locate=0
      end if

      return
      end

c******************************************************************

	subroutine heap_init2(n,ix)

	implicit none

	integer n
	integer ix(2,1)

	integer i

	do i=1,n
	  ix(1,i) = i
	  ix(2,i) = i
	end do

	end

c******************************************************************

	subroutine heap_sort2(n,ra,ix)

	implicit none

	integer n
	real ra(1)
	integer ix(2,1)

	call heap_init2(n,ix)
	call heap_make2(n,ra,ix)
	call heap_retire2(n,ra,ix)

	end

c******************************************************************

	subroutine heap_make2(n,ra,ix)

	implicit none

	integer n
	real ra(1)
	integer ix(2,1)

	integer l

	if( n .lt. 2 ) return

	do l=n/2,1,-1
	  call heap_demote2(l,n,ra,ix)
	end do

	end

c******************************************************************

	subroutine heap_print2(n,ra,ix,text)

	implicit none

	integer n
	real ra(1)
	integer ix(2,1)
	character*(*) text

	integer i

	write(6,*) 'heap print: ',text

	do i=1,n
	  write(6,*) i,n,ra(i),ix(1,i),ix(2,i),ra(ix(1,i))
	end do

	end

c******************************************************************

	subroutine heap_get_index(i,ix,index)

	implicit none

	integer i,index
	integer ix(2,1)

	index = ix(1,i)

	end

c******************************************************************

	subroutine heap_get_val(i,ra,ix,value)

	implicit none

	integer i
	real ra(1)
	integer ix(2,1)
	real value

	value = ra(ix(1,i))

	end

c******************************************************************

	subroutine heap_swap2(i1,i2,ix)

	implicit none

	integer i1,i2
	integer ix(2,1)

	integer ind1,ind2,iaux

	!logical bdebug
	!bdebug = i1 .eq. 3 .and. i2 .eq. 1
	!bdebug = .false.

	ind1 = ix(1,i1)
	ind2 = ix(1,i2)

	iaux = ix(2,ind1)
	ix(2,ind1) = ix(2,ind2)
	ix(2,ind2) = iaux

	ix(1,i1) = ind2
	ix(1,i2) = ind1

	end

c******************************************************************

	subroutine heap_print_index2(n,ix)

	implicit none

	integer n
	integer ix(2,1)

	integer i

	do i=1,n
	  write(6,*) i,ix(1,i),ix(2,i),ix(1,ix(2,i))
	end do

	end

c******************************************************************

	subroutine heap_check2(n,ra,ix)

	implicit none

	integer n
	real ra(1)
	integer ix(2,1)

	call heap_check_index2(n,ix)
	call heap_check_value2(n,ra,ix)

	end

c******************************************************************

	subroutine heap_check_index2(n,ix)

	implicit none

	integer n
	integer ix(2,1)

	integer i,j

	do i=1,n
	  j = ix(1,ix(2,i))
	  if( i .ne. j ) goto 98
	  j = ix(2,ix(1,i))
	  if( i .ne. j ) goto 98
	end do

	return
   98   continue
        write(6,*) n,i,(ix(j,i),j=1,2),ix(1,ix(2,i)),ix(2,ix(1,i))
        stop 'error stop heap_check_index2: index inconsistency'
	end

c******************************************************************

	subroutine heap_check_value2(n,ra,ix)

	implicit none

	integer n
	real ra(1)
	integer ix(2,1)

	integer i,j

	real getval
	getval(i) = ra(ix(1,i))

	do i=1,n/2
	  j = i+i
	  if( getval(i) .lt. getval(j) ) goto 99
	  j = j + 1
	  if( j .le. n .and. getval(i) .lt. getval(j) ) goto 99
	end do

	return
   99   continue
        write(6,*) n,i,j,getval(i),getval(j)
        stop 'error stop heap_check_value2: heap property violated'
	end

c******************************************************************

	subroutine heap_insert2(n,ra,ix,raa)

	implicit none

	integer n
	real ra(1)
	integer ix(2,1)
	real raa

	n = n + 1
	ra(n) = raa
	ix(1,n) = n
	ix(2,n) = n

	call heap_promote2(n,n,ra,ix)
	!call heap_check2(n,ra,ix)

	end

c******************************************************************

	subroutine heap_retire2(n,ra,ix)

	implicit none

	integer n
	real ra(1)
	integer ix(2,1)

	integer ir

	do ir=n,2,-1
	  call heap_swap2(1,ir,ix)
	  call heap_demote2(1,ir-1,ra,ix)
	  !call heap_print2(n,ra,ix,'retire')
	end do

	end

c******************************************************************

	subroutine heap_adjust2(l,n,ra,ix)

c adjusts entry l to right place

	implicit none

	integer l,n
	real ra(1)
	integer ix(2,1)

	call heap_promote2(l,n,ra,ix)
	call heap_demote2(l,n,ra,ix)

	end

c******************************************************************

	subroutine heap_promote2(l,n,ra,ix)

c promotes entry l entry to right place

	implicit none

	integer l,n
	real ra(1)
	integer ix(2,1)

	integer i,j
	real rra,rraj

	real getval
	getval(i) = ra(ix(1,i))

	i = l			!adjust this index
	j = l/2			!node to compare with
	rra = getval(i)		!value to promote

	do while( j .ge. 1 )

	  rraj = getval(j)
	  if( rra .gt. rraj ) then	!promote rra
	    call heap_swap2(i,j,ix)
	    i = j
	    j = j/2
	  else				!finished
	    j = 0
	  end if

	end do

	end

c******************************************************************

	subroutine heap_demote2(l,n,ra,ix)

c demotes entry l to right place

	implicit none

	integer l,n
	real ra(1)
	integer ix(2,1)

	integer i,j
	real rra,rraj,rraj1

	real getval
	getval(i) = ra(ix(1,i))

	i = l			!adjust this index
	j = l+l			!node to compare with
	rra = getval(i)		!value to demote

	do while( j .le. n )

	  rraj = getval(j)
	  if( j .lt. n ) then
	    rraj1 = getval(j+1)
	    if( rraj .lt. rraj1 ) then
		j = j + 1
		rraj = rraj1
	    end if
	  end if

	  if( rra .lt. rraj ) then	!demote rra
	    call heap_swap2(i,j,ix)
	    i = j
	    j = j + j
	  else				!finished
	    j = n + 1
	  end if

	end do

	end

c******************************************************************

      function heap_random()

      implicit none

      real heap_random

      real heap_ran0
      real heap_ran

      integer idum
      save idum
      data idum /-356/

      heap_random = heap_ran(idum)

      end

c******************************************************************

      function heap_ran0(idum)
      dimension v(97)
      save y,v
      data iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        iseed=abs(idum)
        idum=iseed
        do 11 j=1,97
          dum=heap_ran(iseed)
11      continue
        do 12 j=1,97
          v(j)=heap_ran(iseed)
12      continue
        y=heap_ran(iseed)
      endif
      j=1+int(97.*y)
      if(j.gt.97.or.j.lt.1) stop 'error stop ran0: internal error'
      y=v(j)
      ran0=y
      v(j)=heap_ran(iseed)
      return
      end

      function heap_ran(idum)
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
      heap_ran=r(j)
      r(j)=(float(ix1)+float(ix2)*rm2)*rm1
      return
      end

c******************************************************************

	subroutine heap_test

	implicit none

	integer ndim
	parameter (ndim=200000)
	!parameter (ndim=20)
	!parameter (ndim=200)
	!parameter (ndim=80000)

	real ra(ndim)
	real ra1(ndim)
	real ra2(ndim)
	real ra3(ndim)
	integer ix(2,ndim)
	integer ix1(2,ndim)
	integer ix2(2,ndim)
	integer ix3(2,ndim)
	real raa,raa1
	integer i,n,j
	logical berror

	real heap_random

	j = 0
	n = ndim
	berror = .false.

	do i=1,n
	  raa = heap_random()
	  ra(i) = raa
	  ra1(i) = raa
	  ra2(i) = raa
	  !write(6,*) i,ra(i)
	  call heap_insert2(j,ra3,ix3,raa)
	  !call heap_check2(j,ra3,ix3)
	end do

	call heap_sort2(n,ra1,ix1)
	call sort(n,ra2)
	call heap_retire2(n,ra3,ix3)

	!write(6,*) 1,ra1(1),ra2(1)
	do i=2,n
	  call heap_get_val(i,ra1,ix1,raa)
	  call heap_get_val(i-1,ra1,ix1,raa1)
	  if( raa .lt. raa1 ) then
	    !write(6,*) 'not greater... ',i,raa,raa1
	    berror = .true.
	  end if
	  raa1 = ra2(i)
	  if( raa .ne. raa1 ) then
	    !write(6,*) 'not equal... ',i,raa,raa1
	    berror = .true.
	  end if
	  call heap_get_val(i,ra3,ix3,raa1)
	  if( raa .ne. raa1 ) then
	    !write(6,*) 'not equal... ',i,raa,raa1
	    berror = .true.
	  end if
	  !write(6,*) i,raa,ra2(i),raa1
	end do

	if( berror ) then
	  write(6,*) 'there have been errors...'
	else
	  write(6,*) 'test passed.'
	end if

	end
	  
c******************************************************************
	call heap_test
	end
c******************************************************************
c******************************************************************

