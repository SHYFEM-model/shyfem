
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

c revision log :
c
c 13.04.1999	ggu	written from scratch
c 23.03.2010    ggu     changed v6.1.1
c 22.11.2011    ggu     changed VERS_6_1_37
c 14.02.2019    ggu     changed VERS_7_5_56
c
c*************************************************************

	programm taylor_test

c tests taylor diagram

	implicit none

c parameters
	integer ndim
	parameter (ndim=40)

	integer n,j
        real sigma_r,sigma_f,sigma_n,r_corr
        real a(6,ndim)

	n = 1

        call read_data('DATA_C1.dat',ndim,n,a)

c-----------------------
	call qopen
c-----------------------

        call taylor_start
	call taylor_set_param(9.,5.5,2.,-1)
	call taylor_set_config(.true.,.true.,.true.)
        call taylor_full
        call taylor_data(8.,0.9,'A')
        call taylor_data(6.,0.6,'B')
        call taylor_point(8.,0.4,'C')
        call taylor_end

        call taylor_start
	call taylor_set_param(9.,5.5,2.,-1)
	call taylor_set_config(.false.,.false.,.false.)
        call taylor_quart
        call taylor_data(8.,0.9,'A')
        call taylor_data(6.,0.6,'B')
        call taylor_point(8.,0.4,'C')
        call taylor_end

        call taylor_start
	call taylor_set_param(1.3,1.0,0.2,1)  !normalized
	call taylor_set_config(.false.,.false.,.false.)
        call taylor_full
        do j=1,n
          sigma_r = a(2,j)
          sigma_f = a(4,j)
          sigma_n = sigma_f / sigma_r
          r_corr  = a(5,j)
          call taylor_data(sigma_n,r_corr,'X')
        end do
        call taylor_end

        call taylor_start
	call taylor_set_param(1.3,1.0,0.2,1)  !normalized
	call taylor_set_config(.true.,.true.,.true.)
        call taylor_quart
        do j=1,n
          sigma_r = a(2,j)
          sigma_f = a(4,j)
          sigma_n = sigma_f / sigma_r
          r_corr  = a(5,j)
          call taylor_data(sigma_n,r_corr,'X')
        end do
        call taylor_end

c-----------------------
	call qclose
c-----------------------

	end

c*****************************************************************

        subroutine read_data(file,ndim,n,a)

        implicit none

        character*(*) file
        integer ndim,n
        real a(6,ndim)
        character*10 line

        integer i
        real aux(6)

        n = 0
        open(1,file=file)
    3   continue
        read(1,'(6(f8.4,2x),a)',end=4) (aux(i),i=1,6),line
        n = n + 1
        if( n .gt. ndim ) stop 'ndim.....'
        write(6,*) aux
        do i=1,6
          a(i,n) = aux(i)
        end do
        goto 3
    4   continue

        end

c*****************************************************************

