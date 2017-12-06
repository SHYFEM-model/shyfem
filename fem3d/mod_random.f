
!============================================================
	module mod_random
!============================================================

	implicit none

	logical, save :: binit = .false.

	integer, parameter :: int64 = 8

        INTERFACE irand
        MODULE PROCEDURE         irand_max,irand_min_max
        END INTERFACE

!============================================================
	contains
!============================================================

	    subroutine init_random_seed()
            !use iso_fortran_env, only: int64
            implicit none
            integer, allocatable :: seed(:)
            integer :: i, n, un, istat, dt(8), pid
            integer(int64) :: t
            integer getpid
          
	    if( binit ) return
	    binit = .true.

            call random_seed(size = n)
            allocate(seed(n))

            ! First try if the OS provides a random number generator

	    un = 0
	    istat = -1
!            open(newunit=un, file="/dev/urandom", access="stream" 
!     +                 	,form="unformatted", action="read"
!     +			,status="old", iostat=istat)

            if (istat == 0) then
               read(un) seed
               close(un)
            else
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
               call system_clock(t)
               if (t == 0) then
                  call date_and_time(values=dt)
                  t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 
     +                       + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 
     +                       + dt(3) * 24_int64 * 60 * 60 * 1000 
     +                       + dt(5) * 60 * 60 * 1000 
     +                       + dt(6) * 60 * 1000 + dt(7) * 1000 
     +                       + dt(8)
               end if
               pid = getpid()
               t = ieor(t, int(pid, kind(t)))
               do i = 1, n
                  seed(i) = lcg(t)
               end do
            end if
            call random_seed(put=seed)

          end subroutine init_random_seed

!***************************************************************

            ! This simple PRNG might not be good enough for real work, but is
            ! sufficient for seeding a better PRNG.
            function lcg(s)
              !use iso_fortran_env, only: int64
              integer :: lcg
              integer(int64) :: s
              if (s == 0) then
                 s = 104729
              else
                 s = mod(s, 4294967296_int64)
              end if
              s = mod(s * 279470273_int64, 4294967291_int64)
              lcg = int(mod(s, int(huge(0), int64)), kind(0))
            end function lcg

!***************************************************************

        integer function irand_max(max)
        implicit none
        integer max
        real r
        call random_number(r)
        irand_max = 1 + floor(r*max)
        end

        integer function irand_min_max(min,max)
        implicit none
        integer min,max
        real r
        call random_number(r)
        irand_min_max = min + floor(r*(max-min+1))
        end

        logical function rand_prob(ptrue)
        implicit none
        real ptrue
        real r
        call random_number(r)
        rand_prob = ( r < ptrue )
        end

!============================================================
	end module mod_random
!============================================================

	subroutine test_mod_rand

	use mod_random

	implicit none

	integer i

	call init_random_seed

	do i=1,20
	  write(6,*) irand_max(10)
	end do

	end

!***************************************************************

!	program test_mod_rand_main
!	call test_mod_rand
!	end

!***************************************************************

