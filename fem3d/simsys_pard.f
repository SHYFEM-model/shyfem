c
c revision log :
c
c 12.01.2009	ggu	new file for system routines
c 31.03.2009    ggu     call renamed to pard_*
c
c******************************************************************

        subroutine system_initialize

	use mod_system
	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        include 'param.h'

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

        include 'param.h'

	call pard_init_system

	end

c******************************************************************

	subroutine system_solve_z(n,z)

	!use mod_system

	implicit none

        integer n
        real z(n)

	call pard_solve_system(n)

	end

c******************************************************************

	subroutine system_assemble(ie,n,m,kn,mass,rhs)

	use mod_system

	implicit none

	integer ie,n,m
	integer kn(3)
	real mass(3,3)
	real rhs(3)

        include 'param.h'

	integer i,j,kk

	integer loclp,loccoo
	external loclp,loccoo

        do i=1,3
          do j=1,3
            kk=loccoo(kn(i),kn(j),n,m)
            if(kk.gt.0) coo(kk) = coo(kk) + mass(i,j)
          end do
          rvec(kn(i)) = rvec(kn(i)) + rhs(i)
        end do

	end

c******************************************************************

        subroutine system_adjust_z(n,z)

	use mod_system

        implicit none

	integer n
	real z(n)

        include 'param.h'

        integer k

        do k=1,n
          z(k) = rvec(k)
        end do

        end

c******************************************************************

        subroutine system_add_rhs(dt,n,array)

	use mod_system

        implicit none

        real dt
	integer n
        real array(n)

        include 'param.h'

        integer k

        do k=1,n
          rvec(k) = rvec(k) + dt * array(k)
        end do

        end

c******************************************************************

