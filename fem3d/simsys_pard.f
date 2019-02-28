
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

c system routines for Pardiso solver
c
c revision log :
c
c 12.01.2009	ggu	new file for system routines
c 31.03.2009    ggu     call renamed to pard_*
c 15.12.2015    dbf     adjusted for new 3d framework
c
c******************************************************************

        subroutine system_initialize

	use mod_system
	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        write(6,*) '----------------------------------------'
        write(6,*) 'initializing matrix inversion routines'
        write(6,*) 'using Pardiso routines'
        write(6,*) '----------------------------------------'

        call mod_system_init(nkn,nel,ngr,mbw,nlv)

        end

c******************************************************************

	subroutine system_init

	use mod_system

	implicit none

	call pard_init_system

	end

c******************************************************************

        subroutine system_set_explicit

        use mod_system

        implicit none

        bsysexpl = .true.

        end

c******************************************************************

	subroutine system_solve_z(n,z)

	use mod_system 

	implicit none

        integer n
        real z(n)

	call pard_solve_system(.false.,n2max,n,z)

	end

!******************************************************************

	subroutine system_solve_3d(n,nlvdi,nlv,z)

! solves system - z is used for initial guess

	use mod_system

	implicit none

        integer n,nlvdi,nlv
	real z(nlvdi,n)

	integer i,k,l
        real p((nlv+2)*n)
        integer nn 

        i = 0

	do k=1,n
	  i = i + 1
          p(i) = 0.
	  do l=1,nlv
            i = i + 1
            p(i) = z(l,k)
	  end do
	  i = i + 1
	  p(i) = 0.
        end do

	nn = n*(nlv + 2)
	call pard_solve_system(.true.,n3max,nn,p)

	end

c******************************************************************

	subroutine system_assemble(ie,kn,mass,rhs)

	use mod_system

	implicit none

	integer ie,n,m
	integer kn(3)
	real mass(3,3)
	real rhs(3)

	integer i,j,kk

        do i=1,3
          do j=1,3
	    kk=ijp_ie(i,j,ie)
	    if(kk.gt.0) c2coo(kk) = c2coo(kk) + mass(i,j)
          end do
          rvec2d(kn(i)) = rvec2d(kn(i)) + rhs(i)
        end do

	end

c******************************************************************

	subroutine system_assemble_3d(ie,l,nlv,kn,mass,rhs)

! assembles element matrix into system matrix

	use mod_system

	implicit none

        integer ie,l,nlv
        integer kn(3)
        real mass(-1:1,3,3)
        real rhs(3)

        integer i,j,kk

        integer loccoo3d
        external loccoo3d

        do i=1,3
	   do j=1,3
	     kk = loccoo3d(i,j,kn,l,ie)
	     if(kk.gt.0) then
	        c3coo(kk-1) = c3coo(kk-1) + mass(-1,i,j)
		c3coo(kk) = c3coo(kk) + mass(0,i,j)
		c3coo(kk+1) = c3coo(kk+1) + mass(+1,i,j)
	     endif
             if(j.eq.1) rvec3d(i3coo(kk)) = rvec3d(i3coo(kk)) + rhs(i)
	   end do
	enddo
	
        end

c******************************************************************

        subroutine system_adjust_z(n,z)

	use mod_system

        implicit none

	integer n
	real z(n)

        integer k

        z = real(rvec2d)

        end

c******************************************************************

	subroutine system_adjust_3d(n,nlvdi,nlv,z)

! copies solution back to z

	use mod_system

	implicit none

	integer n,nlvdi,nlv
	real z(nlvdi,n)

	integer i,k,l

	i = 0
        do k=1,n
          i = i + 1
          do l=1,nlv
            i = i + 1
            z(l,k) = real(rvec3d(i))
          end do
          i = i + 1
        end do

	end

c******************************************************************

        subroutine system_add_rhs(dt,n,array)

	use mod_system

        implicit none

        real dt
	integer n
        real array(n)

        integer k

        do k=1,n
          rvec2d(k) = rvec2d(k) + dt * array(k)
        end do

        end

c******************************************************************

        subroutine system_adjust_matrix_3d

	implicit none

	call coo_adjust

        end

c******************************************************************

