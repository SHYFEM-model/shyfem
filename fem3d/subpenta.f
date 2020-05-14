
!--------------------------------------------------------------------------
!
!    Copyright (C) 2018-2019  Georg Umgiesser
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

! revision log :
!
! 21.12.2018	ggu	changed VERS_7_5_53
! 18.01.2019	ggu	changed VERS_7_5_55
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!**************************************************************************

! for PENTA see the following information:
! https://www.hindawi.com/journals/mpe/2015/232456/
! Mathematical Problems in Engineering
! Volume 2015, Article ID 232456, 9 pages
! http://dx.doi.org/10.1155/2015/232456
! http://www.math.uakron.edu/~kreider/anpde/penta.f

! for penta_fact and penta_solve see the following information:
! http://www.academia.edu/24859989/
!		On_the_inverse_of_a_general_pentadiagonal_matrix

! please note that the subdiagonals are defined differently for 
! both types of subroutines
!
! the code for tri-diagonal solvers is taken from common text books

!**********************************************************************
!**********************************************************************
!**********************************************************************
! penta-diagonal solvers
!**********************************************************************
!**********************************************************************
!**********************************************************************

      SUBROUTINE PENTA(N,E,A,D,C,F,B,X)

! solves directly pentadiagonal matrix

!   RESULTS:  matrix has 5 bands, EADCF, with D being the main diagonal,
!   E and A are the lower diagonals, and C and F are the upper diagonals.

!     E is defined for rows i = 3:N, but is defined as E(1) to E(N-2)
!     A is defined for rows i = 2:N, but is defined as A(1) to A(N-1)
!     D is defined for rows i = 1:N
!     C is defined for rows i = 1:N-1, but the last element isn't used
!     F is defined for rows i = 1:N-2, but the last 2 elements aren't used

!   B is the right-hand side
!   X is the solution vector

      implicit none

      integer i,n
      double precision E(N),A(N),D(N),C(N),F(N),B(N),X(N),XMULT

      DO I = 2,N-1
        XMULT = A(I-1)/D(I-1)
        D(I) = D(I) - XMULT*C(I-1)
        C(I) = C(I) - XMULT*F(I-1)
        B(I) = B(I) - XMULT*B(I-1)
        XMULT = E(I-1)/D(I-1)
        A(I) = A(I) - XMULT*C(I-1)
        D(I+1) = D(I+1) - XMULT*F(I-1)
        B(I+1) = B(I+1) - XMULT*B(I-1)
      END DO

      XMULT = A(N-1)/D(N-1)
      D(N) = D(N) - XMULT*C(N-1)
      X(N) = (B(N) - XMULT*B(N-1))/D(N)
      X(N-1) = (B(N-1) - C(N-1)*X(N))/D(N-1)

      DO I = N-2,1,-1
        X(I) = (B(I) - F(I)*X(I+2) - C(I)*X(I+1))/D(I)
      END DO

      END

!**********************************************************************

	subroutine penta_solve(n,val,r,s)

! solves already factorized pentadiagonal system

	implicit none

	integer n
	double precision val(-2:2,n)	!matrix
	double precision r(n)		!right hand side
	double precision s(n)		!solution

	integer i
	double precision g(n)
	double precision, parameter :: zero = 0.0d+0
	integer, parameter :: ia = 0
	integer, parameter :: ix = 0
	integer, parameter :: ib = 1
	integer, parameter :: iy = 1
	integer, parameter :: ic = 2
	integer, parameter :: id = -1
	integer, parameter :: iz = -1
	integer, parameter :: ie = -2

	g(1) = r(1)
	g(2) = r(2) - val(iz,2)*g(1)
	do i=3,n
	  g(i) = r(i) - val(iz,i)*g(i-1) - val(ie,i)*g(i-2)
	end do

	s(n) = g(n) / val(ix,n)
	s(n-1) = (g(n-1)-val(iy,n-1)*s(n))/val(ix,n-1)
	do i=n-2,1,-1
	  s(i) = (g(i)-val(iy,i)*s(i+1)-val(ic,i)*s(i+2))/val(ix,i)
	end do

	end

!**********************************************************************

	subroutine penta_factf(n,val)

! factorization of pentadiagonal matrix (faster version)

	implicit none

	integer n
	double precision val(-2:2,n)

	integer i,n1,n2
	double precision xx
	double precision, parameter :: zero = 0.0d+0
	integer, parameter :: ia = 0
	integer, parameter :: ix = 0
	integer, parameter :: ib = 1
	integer, parameter :: iy = 1
	integer, parameter :: ic = 2
	integer, parameter :: id = -1
	integer, parameter :: iz = -1
	integer, parameter :: ie = -2

	if( n < 4 ) then
	  write(6,*) 'penta n = ',n
	  write(6,*) val
	  stop 'error stop penta_factf: n<4'
	end if

	i = 1
	xx = val(ia,1)
	if( xx == zero ) goto 99
	val(iz,2) = val(id,2) / xx
	val(ie,3) = val(ie,3) / xx

	i = 2
	xx = val(ia,2) - val(iy,1)*val(iz,2)
	if( xx == zero ) goto 99
	val(ix,2) = xx
	val(iy,2) = val(ib,2) - val(iz,2)*val(ic,1)
	val(iz,3) = (val(id,3)-val(ie,3)*val(iy,1)) / xx
	val(ie,4) = val(ie,4) / xx

	do i=3,n-2
	  xx = val(ia,i) - val(iy,i-1)*val(iz,i) - val(ie,i)*val(ic,i-2)
	  if( xx == zero ) goto 99
	  val(ix,i) = xx
	  val(iy,i) = val(ib,i) - val(iz,i)*val(ic,i-1)
	  val(iz,i+1) = (val(id,i+1)-val(ie,i+1)*val(iy,i-1)) / xx
	  val(ie,i+2) = val(ie,i+2) / xx
	end do

	n1 = n - 1
	n2 = n - 2
	i = n1
	xx = val(ia,n1) - val(iy,n2)*val(iz,n1) - val(ie,n1)*val(ic,n-3)
	if( xx == zero ) goto 99
	val(ix,n1) = xx
	val(iy,n1) = val(ib,n1) - val(iz,n1)*val(ic,n2)
	val(iz,n) = (val(id,n)-val(ie,n)*val(iy,n2)) / xx

	i = n
	xx = val(ia,n) - val(iy,n1)*val(iz,n) - val(ie,n)*val(ic,n2)
	if( xx == zero ) goto 99
	val(ix,n) = xx

	return
   99	continue
	write(6,*) i,xx
	stop 'error stop penta_factf: matrix singular'
	end

!**********************************************************************

	subroutine penta_fact(n,val)

! factorization of pentadiagonal matrix

	implicit none

	integer n
	double precision val(-2:2,n)

	integer i
	double precision x(n),y(n),z(n)
	double precision xx
	double precision, parameter :: zero = 0.0d+0
	integer, parameter :: ia = 0
	integer, parameter :: ib = 1
	integer, parameter :: ic = 2
	integer, parameter :: id = -1
	integer, parameter :: ie = -2

	if( n < 4 ) then
	  write(6,*) 'penta n = ',n
	  write(6,*) val
	  stop 'error stop penta_fact: n<4'
	end if

	z(1) = zero

	i = 1
	xx = val(ia,i)
	if( xx == zero ) goto 99
	x(i) = xx
	y(i) = val(ib,i)
	z(i+1) = val(id,i+1) / xx
	val(ie,i+2) = val(ie,i+2) / xx

	i = 2
	xx = val(ia,i) - y(i-1)*z(i)
	if( xx == zero ) goto 99
	x(i) = xx
	y(i) = val(ib,i) - z(i)*val(ic,i-1)
	z(i+1) = (val(id,i+1)-val(ie,i+1)*y(i-1)) / xx
	val(ie,i+2) = val(ie,i+2) / xx

	do i=3,n-2
	  xx = val(ia,i) - y(i-1)*z(i) - val(ie,i)*val(ic,i-2)
	  if( xx == zero ) goto 99
	  x(i) = xx
	  y(i) = val(ib,i) - z(i)*val(ic,i-1)
	  z(i+1) = (val(id,i+1)-val(ie,i+1)*y(i-1)) / xx
	  val(ie,i+2) = val(ie,i+2) / xx
	end do

	i = n - 1
	xx = val(ia,i) - y(i-1)*z(i) - val(ie,i)*val(ic,i-2)
	if( xx == zero ) goto 99
	x(i) = xx
	y(i) = val(ib,i) - z(i)*val(ic,i-1)
	z(i+1) = (val(id,i+1)-val(ie,i+1)*y(i-1)) / xx

	i = n
	xx = val(ia,i) - y(i-1)*z(i) - val(ie,i)*val(ic,i-2)
	if( xx == zero ) goto 99
	x(i) = xx
	y(i) = zero

	val(ia,:) = x
	val(ib,:) = y
	val(id,:) = z

	return
   99	continue
	write(6,*) i,xx
	stop 'error stop penta_fact: matrix singular'
	end

!**********************************************************************
!**********************************************************************
!**********************************************************************
! tri-diagonal solvers
!**********************************************************************
!**********************************************************************
!**********************************************************************

	subroutine tria(n,v,r,s)

! solves tridiagonal linear system

	implicit none 

	integer n
	double precision v(3,n)	!lower, central, and upper diagonal values
	double precision r(n)	!right hand side
	double precision s(n)	!solution

	double precision aux
	integer i

! forward elimination

	do i=2,n
  	  aux=v(1,i)/v(2,i-1)
  	  v(2,i)=v(2,i)-aux*v(3,i-1)
  	  r(i)=r(i)-aux*r(i-1)
	end do

! back substitution

	s(n) = r(n)/v(2,n)

	do i=n-1,1,-1
   	  s(i) = (r(i)- v(3,i)*s(i+1))/v(2,i)
	end do

	end

!**********************************************************************

	subroutine tria_multi(n,m,v,r,s)

! solves multiple tridiagonal linear systems

	implicit none 

	integer n		!size of system
	integer m		!how many systems to solve
	double precision v(3,n)	!lower, central, and upper diagonal values
	double precision r(n,m)	!right hand sides
	double precision s(n,m)	!solutions

	double precision aux
	integer i,j

! forward elimination

	do i=2,n
  	  aux=v(1,i)/v(2,i-1)
  	  v(2,i)=v(2,i)-aux*v(3,i-1)
	  do j=1,m
  	    r(i,j)=r(i,j)-aux*r(i-1,j)
	  end do
	end do

! back substitution

	do j=1,m
	  s(n,j) = r(n,j)/v(2,n)
	  do i=n-1,1,-1
   	    s(i,j) = (r(i,j)- v(3,i)*s(i+1,j))/v(2,i)
	  end do
	end do

	end

!**********************************************************************
!**********************************************************************
!**********************************************************************
! testing
!**********************************************************************
!**********************************************************************
!**********************************************************************

	subroutine penta_example_2

	implicit none

	integer, parameter :: n = 10
	
	double precision :: d(n) = (/1, 2, 3, -4, 5, 6, 7, -1, 1, 8/)
	double precision :: u1(n) = (/2, 2, 1, 5, -7, 3, -1, 4, 5, 0/)
	double precision :: u2(n) = (/1, 5, -2, 1, 5, 2, 4, -3, 0, 0/)
	double precision :: l1(n) = (/0, 3, 2, 1, 2, 1, 2, 1, -2, 4/)
	double precision :: l2(n) = (/0, 0, 1, 3, 1, 5, 2, 2, 2, -1/)
	!double precision :: l1(n) = (/ 3, 2, 1, 2, 1, 2, 1, -2, 4, 0/)
	!double precision :: l2(n) = (/ 1, 3, 1, 5, 2, 2, 2, -1, 0, 0/)
	double precision :: r(n) = 
     +			(/8, 33, 8, 24, 29, 98, 99, 17, 57, 108/)
	double precision :: t(n) = 
     +			(/1, 2, 3, 4, 5, 6, 7, 8, 9, 10/)
	double precision :: s(n)
	double precision :: val(-2:2,n)
	double precision :: val1(-2:2,n)
	double precision :: ss

	val(-2,:) = l2
	val(-1,:) = l1
	val(0,:) = d
	val(1,:) = u1
	val(2,:) = u2

	val1 = val

	call penta_fact(n,val)
	call penta_solve(n,val,r,s)

	!write(6,*) s
	call band_check(n,s,t,ss)
	write(6,*) n,ss

	call penta_factf(n,val1)
	call penta_solve(n,val1,r,s)

	!write(6,*) s
	call band_check(n,s,t,ss)
	write(6,*) n,ss

	end

!**********************************************************************

	subroutine penta_example_1

	implicit none

	integer, parameter :: n = 10
	
	double precision :: d(n) = (/1, 2, 3, -4, 5, 6, 7, -1, 1, 8/)
	double precision :: u1(n) = (/2, 2, 1, 5, -7, 3, -1, 4, 5, 0/)
	double precision :: u2(n) = (/1, 5, -2, 1, 5, 2, 4, -3, 0, 0/)
	!double precision :: l1(n) = (/0, 3, 2, 1, 2, 1, 2, 1, -2, 4/)
	!double precision :: l2(n) = (/0, 0, 1, 3, 1, 5, 2, 2, 2, -1/)
	double precision :: l1(n) = (/ 3, 2, 1, 2, 1, 2, 1, -2, 4, 0/)
	double precision :: l2(n) = (/ 1, 3, 1, 5, 2, 2, 2, -1, 0, 0/)
	double precision :: r(n) = 
     +			(/8, 33, 8, 24, 29, 98, 99, 17, 57, 108/)
	double precision :: t(n) = 
     +			(/1, 2, 3, 4, 5, 6, 7, 8, 9, 10/)
	double precision :: x(n)
	double precision :: ss

        call PENTA(N,l2,l1,d,u1,u2,r,x)

	!write(6,*) x
	call band_check(n,x,t,ss)
	write(6,*) n,ss

	end

!**********************************************************************

	subroutine tri_example_1

	implicit none

	integer, parameter :: n = 2
	
	double precision :: l1(n) = (/0, 3/)
	double precision :: d(n) = (/1, 4/)
	double precision :: u1(n) = (/2, 0/)
	double precision :: r(n) = (/5,11/)
	double precision :: t(n) = (/1,2/)

	double precision :: s(n)
	double precision :: val(-1:1,n)
	double precision :: ss

	val(-1,:) = l1
	val(0,:) = d
	val(1,:) = u1

	call tria(n,val,r,s)

	!write(6,*) s
	call band_check(n,s,t,ss)
	write(6,*) n,ss

	end

!**********************************************************************

	subroutine band_check(n,s,t,ss)

	implicit none

	integer n
	double precision s(n)
	double precision t(n)	!known and expected solution

	integer i
	double precision ss

	ss = 0.

	do i=1,n
	  ss = ss + (t(i)-s(i))**2
	end do

	ss = sqrt( ss / n )
	!write(6,*) n,sqrt(ss)

	end

!**********************************************************************

	subroutine test_random(nloop,ndim,m)

! executes random tests for banded matrices (penta or tria)

	implicit none

	integer :: nloop		!total number of loops
	integer :: ndim			!max size of system
	integer :: m			!upper or lower diagonals
					!2 for penta, 1 for tria

	double precision :: val(-m:m,ndim)
	double precision :: s(ndim)
	double precision :: r(ndim)
	double precision :: x(ndim)

	integer :: l,i,j,n
	double precision :: rt,rr
	logical :: bint = .false.

	integer rand_int, rand_sign
	double precision rand_double
	double precision :: ss
	double precision :: ssmax = 0.

	do l=1,nloop
	  n = rand_int(2*m,ndim)		!ensure minimum number of diags
	  do i=1,n
	    rt = 0.
	    do j=-m,m
	      rr = rand_double(-1,+1)
	      if( bint ) rr = int(10*rr)
	      val(j,i) = rr
	      rt = rt + abs(rr)
	    end do
	    val(0,i) = rt * rand_sign()		!ensure diagonally dominance
	    s(i) = rand_double(-1,+1)
	    if( bint ) s(i) = int(10*rand_double(-1,+1))
	  end do
	  call band_multiply(n,m,val,s,r)

	  if( bint ) then
	    write(6,*) n,m
	    do i=1,n
	      write(6,1000) nint(val(:,i))
     +				,nint(s(i)),nint(r(i))
	    end do
 1000	    format(9i8)
	  end if

	  if( m == 2 ) then
	    !call penta_fact(n,val)
	    call penta_factf(n,val)
	    call penta_solve(n,val,r,x)
	  else if( m == 1 ) then
	    !call tria(n,val,r,x)
	    call tria_multi(n,1,val,r,x)
	  else
	    stop 'error stop test_random: m must be 1 or 2'
	  end if

	  if( bint ) write(6,1000) nint(x(1:n))

	  call band_check(n,s,x,ss)
	  ssmax = max(ssmax,ss)
	  !write(6,*) l,n,ss

	end do 

	write(6,*) 'm = ',2*m+1,'  nloop = ',nloop,'  ssmax = ',ssmax

	end

!**********************************************************************

	function rand_double(min,max)

	implicit none

	double precision rand_double
	integer min,max

	double precision r

	call random_number(r)

	rand_double = min + (max-min)*r

	end

!**********************************************************************

	function rand_int(min,max)

	implicit none

	integer rand_int
	integer min,max

	double precision r

	call random_number(r)

	rand_int = min + (1+max-min)*r

	end

!**********************************************************************

	function rand_sign()

	implicit none

	integer rand_sign

	double precision r

	call random_number(r)

	rand_sign = +1
	if( r < 0.5 ) rand_sign = -1

	end

!**********************************************************************

	subroutine band_multiply(n,m,val,s,r)

! multiplies band matrix with vector

	implicit none

	integer n,m
	double precision :: val(-m:m,n)
	double precision :: s(n)
	double precision :: r(n)

	integer i,j,jmin,jmax
	double precision acu

	do i=1,n
	  jmin = max(-m,1-i)
	  jmax = min(m,n-i)
	  acu = 0.
	  do j=jmin,jmax
	    acu = acu + val(j,i)*s(i+j)
	  end do
	  r(i) = acu
	end do

	end

!**********************************************************************

	subroutine test_penta
	call penta_example_1
	call penta_example_2
	call tri_example_1
	call test_random(1000,20,2)
	call test_random(1000,20,1)
	!call test_random(1,5,1)
	end

!**********************************************************************
!	program test_penta_main
!	call test_penta
!	end
!**********************************************************************

