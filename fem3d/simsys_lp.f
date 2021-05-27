
!--------------------------------------------------------------------------
!
!    Copyright (C) 2009-2010,2014-2015,2017,2019  Georg Umgiesser
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

c matrix inversion administration (Gaussian elimination - traditional)
c
c revision log :
c
c 12.01.2009	ggu	new file for system routines
c 31.01.2009	ggu	prepared for double precision, new aux vect vs1v,...
c 31.01.2009	ggu	prepared for double precision, new aux vect vs1v,...
c 23.03.2010	ggu	changed v6.1.1
c 28.01.2014	ggu	changed VERS_6_1_71
c 05.06.2015	ggu	changed VERS_7_1_12
c 10.07.2015	ggu	changed VERS_7_1_50
c 17.07.2015	ggu	changed VERS_7_1_52
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 30.07.2015	ggu	changed VERS_7_1_83
c 15.12.2015	ggu	added dummy subroutines for 3d case
c 05.12.2017	ggu	changed VERS_7_5_39
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2021	ggu	added routine system_finalize()
c
c******************************************************************
c
c to change from real to double precision
c change amat here and vs*v in common.h
c
c******************************************************************

	subroutine system_initialize

	use mod_system
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	write(6,*) '----------------------------------------'
	write(6,*) 'initializing matrix inversion routines'
	write(6,*) 'using Gaussian elimination'
	write(6,*) '----------------------------------------'

	call mod_system_init(nkn,nel,ngr,mbw,nlv)
	call mod_system_amat_init(nkn,mbw)

	end

c******************************************************************

	subroutine system_init

	use mod_system
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	!call lp_init_system(nkn,mbw,amat,vs1v)
	call dlp_init_system(nkn,mbw,amat,vs1v)

	end

c******************************************************************

        subroutine system_set_explicit

        use mod_system

        implicit none

        bsysexpl = .true.

        end

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine system_solve_z(n,z)

	use mod_system
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer n
	real z(n)

	!call lp_solve_system(nkn,mbw,amat,vs1v,is2v,vs3v)
	call dlp_solve_system(nkn,mbw,amat,vs1v,is2v,vs3v)

	end

c******************************************************************

	subroutine system_assemble(ie,n,m,kn,mass,rhs)

	use mod_system

	implicit none

	integer ie,n,m
	integer kn(3)
	real mass(3,3)
	real rhs(3)

	integer i,j,kk

	integer loclp,loccoo
	external loclp,loccoo

	!write(6,*) 'amat: ',size(amat)
        do i=1,3
          do j=1,3
            kk=loclp(kn(i),kn(j),n,m)
            if(kk.gt.0) amat(kk) = amat(kk) + mass(i,j)
          end do
          vs1v(kn(i)) = vs1v(kn(i)) + rhs(i)
        end do

	end

c******************************************************************

        subroutine system_adjust_z(n,z)

	use mod_system

        implicit none

	integer n
	real z(n)

        integer k

        do k=1,n
          z(k) = vs1v(k)
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
          vs1v(k) = vs1v(k) + dt * array(k)
        end do

        end

c******************************************************************
c******************************************************************
c******************************************************************
c 3d case ... not ready
c******************************************************************
c******************************************************************
c******************************************************************

	subroutine system_solve_3d(n,nlvdi,nlv,z)

        use mod_system

        implicit none

        integer n,nlvdi,nlv
        real z(nlvdi,n)

	stop 'error stop: gauss elimination 3d not ready'

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

	stop 'error stop: gauss elimination 3d not ready'

	end

c******************************************************************

	subroutine system_adjust_3d(n,nlvdi,nlv,z)

! copies solution back to z

        use mod_system

        implicit none

        integer n,nlvdi,nlv
        real z(nlvdi,n)

	stop 'error stop: gauss elimination 3d not ready'

	end

c******************************************************************

	subroutine system_adjust_matrix_3d

        implicit none

	stop 'error stop: gauss elimination 3d not ready'

	end

c******************************************************************

        subroutine system_finalize

        implicit none

        end

c******************************************************************
c******************************************************************
c******************************************************************

