
!--------------------------------------------------------------------------
!
!    Copyright (C) 2007-2008,2019  Georg Umgiesser
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

c matrix inversion routines (non symmetric band matrix) (Gauss inversion)
c
c contents :
c
c lp_init_system	initializes matrix and vector
c lp_solve_system	factors and solves for solution
c lp_subst_system	solves for solution
c lp_mult_band		multiplies matrix with vector
c
c dlp_init_system	initializes matrix and vector
c dlp_solve_system	factors and solves for solution
c dlp_subst_system	solves for solution
c dlp_mult_band		multiplies matrix with vector
c
c loclp			finds position in matrix
c
c revision log :
c
c 02.04.2007	ggu	assembled from lapack
c 06.06.2007	ggu	new routines for back substitution and initialization
c 22.04.2008	ggu	new SAXPY_NEW for parallelization trial
c 04.01.2019	ggu	linpack and blas routines transfered
c 18.01.2019	ggu	changed VERS_7_5_55
c 16.02.2019	ggu	changed VERS_7_5_60
c
c*************************************************************************
c
c access of unsymmetric banded matrix in linpack
c
c linear matrix must have the following dimension:
c		ndim = ( 3*m + 1 ) * n
c
c*************************************************************************

	subroutine lp_init_system(n,m,abd,b)

c initializes band matrix

	implicit none

	integer n,m
	real abd(1)		!band matrix
	real b(1)		!right hand side [n], at return x

	integer nb,i

	nb = ( 3*m + 1 ) * n

	do i=1,nb
	  abd(i) = 0.
	end do

	do i=1,n
	  b(i) = 0.
	end do

	end

c*************************************************************************

	subroutine lp_solve_system(n,m,abd,b,ipvt,z)

c solves system a*x=b

	implicit none

	integer n,m
	real abd(1)		!band matrix
	real b(1)		!right hand side [n], at return x
	integer ipvt(1)		!pivot information [n]
	real z(1)		!aux vector [n]

	integer lda,ml,mu,job,info
	real rcond

	info = 0
	rcond = 1.

	ml = m
	mu = m
	lda = 3 * m + 1
	job = 0			!we solve direct problem, not transpose

	!call sgbco(abd,lda,n,ml,mu,ipvt,rcond,z)	!rcond should be > 0
        call sgbfa(abd,lda,n,ml,mu,ipvt,info)		!info should be 0
	call sgbsl(abd,lda,n,ml,mu,ipvt,b,job)

	if( info .ne. 0 .or. rcond .eq. 0. ) then
	  write(6,*) 'condition number: ',rcond,info
	  stop 'error stop lp_solve_system: info'
	end if

	end

c*************************************************************************

	subroutine lp_subst_system(n,m,abd,b,ipvt)

c solves a*x=b with a already factored (by lp_solve_system)

	integer n,m
	real abd(1)		!band matrix (already factored)
	real b(1)		!right hand side [n], at return x
	integer ipvt(1)		!pivot information [n]

	integer lda,ml,mu,job

	ml = m
	mu = m
	lda = 3 * m + 1
	job = 0			!we solve direct problem, not transpose

	call sgbsl(abd,lda,n,ml,mu,ipvt,b,job)

	end

c*************************************************************************

	subroutine lp_mult_band(n,m,abd,b,res)

c multiplies band matrix a with vector b and returns result in res

	implicit none

	integer n,m
	real abd(1)		!band matrix
	real b(1)		!right hand side [n]
	real res(1)		!result [n]

	integer i,j,k,mm,j1,j2
	real acu

	integer loclp

	mm = 2 * m + 1

	do i=1,n
	  j1 = max(1,i-m)
	  j2 = min(n,i+m)
	  acu = 0.
	  do j=j1,j2
            k = loclp(i,j,n,m)
	    !write(6,*) i,j,k,j1,j2,abd(k),b(j)
	    acu = acu + abd(k) * b(j)
	  end do
	  res(i) = acu
	end do

	end

c*************************************************************************

        subroutine dlp_init_system(n,m,abd,b)

c initializes band matrix

        implicit none

        integer n,m
        double precision abd(1)      !band matrix
        double precision b(1)        !right hand side [n], at return x

        integer nb,i
	double precision zero

        nb = ( 3*m + 1 ) * n
	zero = 0.

        do i=1,nb
          abd(i) = zero
        end do

        do i=1,n
          b(i) = zero
        end do

        end

c*************************************************************************

	subroutine dlp_solve_system(n,m,abd,b,ipvt,z)

c solves system a*x=b

	implicit none

	integer n,m
	double precision abd(1)		!band matrix
	double precision b(1)		!right hand side [n], at return x
	integer ipvt(1)			!pivot information [n]
	double precision z(1)		!aux vector [n]

	integer lda,ml,mu,job,info
	double precision rcond

        info = 0
        rcond = 1.

	ml = m
	mu = m
	lda = 3 * m + 1
	job = 0			!we solve direct problem, not transpose

	call dgbco(abd,lda,n,ml,mu,ipvt,rcond,z)	!rcond should be > 0
        !call dgbfa(abd,lda,n,ml,mu,ipvt,info)          !info should be 0
	call dgbsl(abd,lda,n,ml,mu,ipvt,b,job)

	if( info .ne. 0 .or. rcond .eq. 0. ) then
	  write(6,*) 'condition number: ',rcond,info
	  stop 'error stop dlp_solve_system: info'
	end if

	end

c*************************************************************************

	subroutine dlp_subst_system(n,m,abd,b,ipvt)

c solves a*x=b with a already factored (by dlp_solve_system)

	integer n,m
	double precision abd(1)		!band matrix (already factored)
	double precision b(1)		!right hand side [n], at return x
	integer ipvt(1)			!pivot information [n]

	integer lda,ml,mu,job

	ml = m
	mu = m
	lda = 3 * m + 1
	job = 0			!we solve direct problem, not transpose

	call dgbsl(abd,lda,n,ml,mu,ipvt,b,job)

	end

c*************************************************************************

	subroutine dlp_mult_band(n,m,abd,b,res)

c multiplies band matrix a with vector b and returns result in res

	implicit none

	integer n,m
	double precision abd(1)		!band matrix
	double precision b(1)		!right hand side [n]
	double precision res(1)		!result [n]

	integer i,j,k,mm,j1,j2
	double precision acu

	integer loclp

	mm = 2 * m + 1

	do i=1,n
	  j1 = max(1,i-m)
	  j2 = min(n,i+m)
	  acu = 0.
	  do j=j1,j2
            k = loclp(i,j,n,m)
	    !write(6,*) i,j,k,j1,j2,abd(k),b(j)
	    acu = acu + abd(k) * b(j)
	  end do
	  res(i) = acu
	end do

	end

c*************************************************************************

        function loclp(i,j,n,m)

c access linpack routines (unsymmetric banded matrix)
c
c (i,j)   position of element in square matrix (row,column)
c n       dimension of square matrix
c m       band width of square matrix
c loclp   position of element in band matrix

        implicit none

        integer loclp
        integer i,j,n,m

	integer lda,k

	loclp = 0
	if( abs(i-j) .gt. m ) return

	lda = 3*m + 1
	k = i - j + lda - m

	loclp = k + (j-1) * lda

        end

c*************************************************************************
c*************************************************************************
c*************************************************************************
c*************************************************************************
c*************************************************************************

	subroutine loclp_test

c use test of linpack routine

	implicit none

	integer n,m,lda,lnmax
	parameter (n=6,m=2,lda=3*m+1,lnmax=lda)

	integer i,j,k
	integer i1,i2
	real val
	real a(n,n)
	real b(lnmax*lnmax)
	real c(lda,n)

	real br(lnmax)
	real z(lnmax)
	integer ipvt(lnmax)

	integer loclp

	write(6,*) 'n,m,lda,lnmax: ',n,m,lda,lnmax
	do j=1,n
	  do i=1,n
	    if( j-i .gt. m ) then
		val = 0.
	    else if( i-j .gt. m-1 ) then
		val = 0.
	    else
		val = 10.*i + j
	    end if
	    a(i,j) = val
	  end do
	end do

	write(6,*) 'matrix a'
	call loclp_print(n,n,a)
	    
	call loclp_init(lnmax,lnmax,b)
	write(6,*) 'matrix b init'
	call loclp_print(lnmax,lnmax,b)

	do j=1,n
	  do i=1,n
	    k = loclp(i,j,n,m)
	    b(k) = a(i,j)
	  end do
	end do

	write(6,*) 'matrix b'
	call loclp_print(lda,n,b)
	    
	write(6,*) 'matrix b through loclp'
	do i=1,n 
	  do j=1,n
	    k = loclp(i,j,n,m)
	    val = 0.
	    if( k .gt. 0 ) val = b(k)
	    br(j) = val
	  end do
	  write(6,*) (br(j),j=1,n)
	end do

	call loclp_init(lda,n,c)

	do j=1,n
	  i1 = max(1,j-m)
	  i2 = min(n,j+m)
	  do i=i1,i2
	    k = i - j + 2*m + 1
	    c(k,j) = a(i,j)
	  end do
	end do

	write(6,*) 'matrix c'
	call loclp_print(lda,n,c)
	    
	do i=1,n
	  z(i) = i
	end do

	call lp_mult_band(n,m,b,z,br)

	write(6,*) 'right hand side'
	do i=1,n
	  write(6,*) i,br(i)
	end do

	call lp_solve_system(n,m,b,br,ipvt,z)

	write(6,*) 'solution system'
	do i=1,n
	  write(6,*) i,br(i)
	end do

	end

c*************************************************************************

	subroutine dloclp_test

c use test of linpack routine

	implicit none

	integer n,m,lda,lnmax
	parameter (n=6,m=2,lda=3*m+1,lnmax=lda)

	integer i,j,k
	integer i1,i2
	double precision val
	double precision a(n,n)
	double precision b(lnmax*lnmax)
	double precision c(lda,n)

	double precision br(lnmax)
	double precision z(lnmax)
	integer ipvt(lnmax)

	integer loclp

	write(6,*) 'n,m,lda,lnmax: ',n,m,lda,lnmax
	do j=1,n
	  do i=1,n
	    if( j-i .gt. m ) then
		val = 0.
	    else if( i-j .gt. m-1 ) then
		val = 0.
	    else
		val = 10.*i + j
	    end if
	    a(i,j) = val
	  end do
	end do

	write(6,*) 'matrix a'
	call dloclp_print(n,n,a)
	    
	call dloclp_init(lnmax,lnmax,b)
	write(6,*) 'matrix b init'
	call dloclp_print(lnmax,lnmax,b)

	do j=1,n
	  do i=1,n
	    k = loclp(i,j,n,m)
	    b(k) = a(i,j)
	  end do
	end do

	write(6,*) 'matrix b'
	call dloclp_print(lda,n,b)
	    
	write(6,*) 'matrix b through loclp'
	do i=1,n 
	  do j=1,n
	    k = loclp(i,j,n,m)
	    val = 0.
	    if( k .gt. 0 ) val = b(k)
	    br(j) = val
	  end do
	  write(6,*) (br(j),j=1,n)
	end do

	call dloclp_init(lda,n,c)

	do j=1,n
	  i1 = max(1,j-m)
	  i2 = min(n,j+m)
	  do i=i1,i2
	    k = i - j + 2*m + 1
	    c(k,j) = a(i,j)
	  end do
	end do

	write(6,*) 'matrix c'
	call dloclp_print(lda,n,c)
	    
	do i=1,n
	  z(i) = i
	end do

	call dlp_mult_band(n,m,b,z,br)

	write(6,*) 'right hand side'
	do i=1,n
	  write(6,*) i,br(i)
	end do

	call dlp_solve_system(n,m,b,br,ipvt,z)

	write(6,*) 'solution system'
	do i=1,n
	  write(6,*) i,br(i)
	end do

	end

c*************************************************************************
c*************************************************************************
c*************************************************************************
c*************************************************************************
c*************************************************************************

	subroutine loclp_init(l,n,a)

	implicit none

	integer l,n
	real a(l,n)

	integer i,j

	do j=1,n
	  do i=1,l
	    a(i,j) = 0.
	  end do
	end do

	end

c*************************************************************************

	subroutine loclp_print(l,n,a)

	implicit none

	integer l,n
	real a(l,n)

	integer i,j

	do i=1,l
	  write(6,*) (a(i,j),j=1,n)
	end do

	end

c*************************************************************************

	subroutine dloclp_init(l,n,a)

	implicit none

	integer l,n
	double precision a(l,n)

	integer i,j

	do j=1,n
	  do i=1,l
	    a(i,j) = 0.
	  end do
	end do

	end

c*************************************************************************

	subroutine dloclp_print(l,n,a)

	implicit none

	integer l,n
	double precision a(l,n)

	integer i,j

	do i=1,l
	  write(6,*) (a(i,j),j=1,n)
	end do

	end

c*************************************************************************
c*************************************************************************
c*************************************************************************
c*************************************************************************
c*************************************************************************

c	program loclp_main
c	call loclp_test
c	call dloclp_test
c	end

c*************************************************************************

