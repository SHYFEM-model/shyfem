c
c revision log :
c
c 12.01.2009	ggu	new file for system routines
c 31.03.2009    ggu     call renamed to pard_*
c
c******************************************************************

        subroutine system_initialize

        implicit none

        include 'common.h'

        write(6,*) '----------------------------------------'
        write(6,*) 'initializing matrix inversion routines'
        write(6,*) 'using Pardiso routines'
        write(6,*) '----------------------------------------'

        end

c******************************************************************

	subroutine system_init

	implicit none

	include 'common.h'

	call pard_init_system

	end

c******************************************************************

	subroutine system_solve_z

	implicit none

	include 'common.h'

	call pard_solve_system

	end

c******************************************************************

	subroutine system_assemble(n,m,kn,mass,rhs)

	implicit none

	integer n,m
	integer kn(3)
	real mass(3,3)
	real rhs(3)

	include 'common.h'

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

        subroutine system_adjust_z

        implicit none

	include 'common.h'

        integer k

        do k=1,nkn
          znv(k) = rvec(k)
        end do

        end

c******************************************************************

        subroutine system_add_rhs(dt,array)

        implicit none

        real dt
        real array(1)

	include 'common.h'

        integer k

        do k=1,nkn
          rvec(k) = rvec(k) + dt * array(k)
        end do

        end

c******************************************************************

