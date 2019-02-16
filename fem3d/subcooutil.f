
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

!*******************************************************************

	subroutine csr_show(n,nnz,icsr,jcsr,ccsr)

	implicit none

	integer n			!row dimension
	integer nnz			!number of non zeros in coo matrix
	integer icsr(n+1)		!start of row of elements
	integer jcsr(nnz)		!column of elements
	double precision ccsr(nnz)	!value of elements

	integer k,i,j,ip,ipe
	double precision v
	character*1 c
	character*132 laux
	character*132 line(132)

	if( n > 132 ) stop 'error stop coo_show: dimension n'

	call coo_init_line(laux,'.')
	line = laux

	do i=1,n
	  ip = icsr(i)
	  ipe = icsr(i+1) - 1
	  do k=ip,ipe
	    j = jcsr(k)
	    v = ccsr(k)
	    call coo_convert_char(v,c)
	    call coo_insert_char(i,j,c,line)
	  end do
	end do

	call coo_print_lines(n,line)

	end

!*******************************************************************

	subroutine coo_print(n,nnz,icoo,jcoo,ccoo,rhs)

	implicit none

	integer n			!row dimension
	integer nnz			!number of non zeros in coo matrix
	integer icoo(nnz)		!row of elements
	integer jcoo(nnz)		!column of elements
	double precision ccoo(nnz)	!value of elements
	double precision rhs(n)		!value of elements

	integer k,i,j
	double precision v
	double precision mat(n,n)
	double precision vec(n)
	double precision aux(n)
	integer iaux(n)

	mat = 0.
	vec = rhs

	do k=1,nnz
	  i = icoo(k)
	  j = jcoo(k)
	  v = ccoo(k)
	  mat(i,j) = v
	end do

	write(6,*) '----------------------'
	do i=1,n
	  write(6,'(7g11.3)') (mat(i,j),j=1,n),vec(i)
	end do
	write(6,*) '----------------------'

	call dmatinv(mat,iaux,aux,n,n)

	write(6,*) '----------------------'
	do i=1,n
	  write(6,'(7g11.3)') (mat(i,j),j=1,n),vec(i)
	end do
	write(6,*) '----------------------'

	do i=1,n
	  v = 0
	  do j=1,n
	    v = v + mat(i,j) * vec(j)
	  end do
	  write(6,*) i,v
	end do

	end

!*******************************************************************

	subroutine coo_show(n,nnz,icoo,jcoo,ccoo)

	implicit none

	integer n			!row dimension
	integer nnz			!number of non zeros in coo matrix
	integer icoo(nnz)		!row of elements
	integer jcoo(nnz)		!column of elements
	double precision ccoo(nnz)	!value of elements

	integer k,i,j
	double precision v
	character*1 c
	character*132 laux
	character*132 line(132)

	if( n > 132 ) stop 'error stop coo_show: dimension n'

	call coo_init_line(laux,'.')
	line = laux

	do k=1,nnz
	  i = icoo(k)
	  j = jcoo(k)
	  v = ccoo(k)
	  call coo_convert_char(v,c)
	  call coo_insert_char(i,j,c,line)
	end do

	call coo_print_lines(n,line)

	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine coo_convert_char(v,c)

	implicit none

	double precision v
	character*1 c

	integer ival,ic

	  if( v == 0. ) then
	    c = '0'
	  else if( v >= 1. .and. v <= 9. ) then
	    ival = nint(v)
	    ic = ival + ichar('0')
	    c = char(ic)
	  else
	    c = '*'
	  end if

	end

!*******************************************************************

	subroutine coo_insert_char(i,j,c,line)

	implicit none

	integer i,j
	character*1 c
	character*132 line(1:132)

	character*132 laux

	laux = line(i)
	laux(j:j) = c
	line(i) = laux

	end

!*******************************************************************

	subroutine coo_init_line(line,c)

	implicit none

	character*132 line
	character*1 c

	integer i

	do i=1,132
	  line(i:i) = c
	end do

	end

!*******************************************************************

	subroutine coo_print_lines(n,line)

	implicit none

	integer n
	character*132 line(1:132)

	integer i
	character*132 laux

	call coo_init_line(laux,'=')
	write(6,*) laux(1:n)

	do i=1,n
	  laux = line(i)
	  write(6,*) laux(1:n)
	end do

	call coo_init_line(laux,'=')
	write(6,*) laux(1:n)

	end

!*******************************************************************

	subroutine coo_test

	implicit none

	integer, parameter :: nnz = 8
	integer, parameter :: n = 5
	integer icoo(nnz)
	integer jcoo(nnz)
	double precision ccoo(nnz)
	integer icsr(n+1)
	integer jcsr(nnz)
	double precision ccsr(nnz)

	integer i

	do i=1,n
	  icoo(i) = i
	  jcoo(i) = i
	  ccoo(i) = i
	end do

	i = 6
	icoo(i) = 2
	jcoo(i) = 3
	ccoo(i) = i

	i = 7
	icoo(i) = 1
	jcoo(i) = 5
	ccoo(i) = i

	i = 8
	icoo(i) = 4
	jcoo(i) = 2
	ccoo(i) = i

	call coo_show(n,nnz,icoo,jcoo,ccoo)
	call coocsr(n,nnz,ccoo,icoo,jcoo,ccsr,jcsr,icsr)
	call coo_show(n,nnz,icoo,jcoo,ccoo)
	call csr_show(n,nnz,icsr,jcsr,ccsr)

	end

!*******************************************************************

!	program coo_test_main
!	call coo_test
!	end

!*******************************************************************

