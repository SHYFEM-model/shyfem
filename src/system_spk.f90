!
! revision log :
!
! 12.01.2009	ggu	new file for system routines
! 31.03.2009	ggu	call renamed to spk_*
! 25.05.2015	ggu	some calls changed (pass array in)
! 09.12.2015	ggu	adapted to new pointers and 3d matrix
! 15.12.2015	ggu&deb	finsihed and validated
!
!******************************************************************
!----------------------------------------------------------------------
        module system
!----------------------------------------------------------------------
        contains
!----------------------------------------------------------------------

        subroutine system_initialize

! allocates data structure

	use system_matrix
	use levels
	use basin, only : nkn,nel,ngr,mbw
        use mpi_common_struct

        implicit none

        write(6,*) '----------------------------------------'
        write(6,*) 'initializing matrix inversion routines'
        write(6,*) 'using Sparskit routines'
        write(6,*) '----------------------------------------'

#ifdef DEBUGON
        call mod_system_init(nkn_local,nel_local,ngr,mbw,nlv,n_threads)
#else
        call mod_system_init(nkn,nel,ngr,mbw,nlv,n_threads)
#endif

        end

!******************************************************************

	subroutine system_init

! on first call initializes pointers - sets arrays to zero
!
! must be called before every assembly of matrix

	use system_matrix
        use sparskit_admin

	implicit none

	if( bsysexpl ) then
	  rvec2d = 0.d0
	  raux2d = 0.d0
	else
	  call spk_init_system
	end if

	end

!******************************************************************

	subroutine system_set_explicit

	use system_matrix

	implicit none

	bsysexpl = .true.

	end

!******************************************************************
!******************************************************************
!******************************************************************
!
!	subroutine system_show_internal
!
!	use basin
!	use system_matrix
!	use shympi
!
!        write(my_unit,*) 'system internal: ',my_id,bsysexpl
!
!        do k=1,nkn
!          write(my_unit,1234) k,ipv(k),id_node(k),xgv(k),ygv(k)
!     +                          ,rvec2d(k),raux2d(k)
! 1234     format(3i5,4f12.4)
!        end do
!
!	end
!
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine system_solve_z(n,z)

! solves system - z is used for initial guess

	use system_matrix
	use shympi
        use petsc_admin
        use basin
        use mpi_io_admin
        use sparskit_admin

	implicit none

        include 'femtime.h'
	integer n
	double precision z(n)
        integer i
               

	if( bsysexpl ) then
#ifndef DEBUGON
          call shympi_exchange_and_sum_2D_nodes(raux2d)
#endif
	  rvec2d = rvec2d / raux2d	!GGUEXPL

	else
          if(b_use_mpi) then    !! using shympi program

#ifdef DEBUGON
            if(n_threads.gt.1) then
               call sol_petsc_debug
            else
               call sol_petsc
            end if
#else
            call sol_petsc
#endif
                 
          else                  !! using shyfem program
	    call spk_solve_system(.false.,n2max,n,z)
          end if
	end if

	end

!******************************************************************

	subroutine system_solve_3d(n,nlvdi,nlv,z)

! solves system - z is used for initial guess

	use system_matrix
        use sparskit_admin

	implicit none

	integer n,nlvdi,nlv
	double precision z(nlvdi,n)

	integer i,k,l
	double precision p((nlv+2)*n)
	integer nn !DEB

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

	nn = n*(nlv + 2) !DEB
	!call spk_solve_system(.true.,n3max,n,p)
	call spk_solve_system(.true.,n3max,nn,p) !DEB

	end

!******************************************************************

	subroutine system_assemble(ie,kn,mass,rhs)

! assembles element matrix into system matrix

	use system_matrix
        use shympi

	implicit none

	integer ie
	integer kn(3)
	double precision mass(3,3)
	double precision rhs(3)
        character*40 formato,fmt2,fmt3

	integer i,j,kk,ieglob

        formato='(A,E24.17,E24.17,I8)'
        fmt2='(A,E24.17,E24.17,I8,I8,I8,I8)'
        fmt3='(A,I8,E24.17,E24.17)'

	if( bsysexpl ) then
          do i=1,3
            raux2d(kn(i)) = raux2d(kn(i)) + mass(i,i)	!GGUEXPL
            rvec2d(kn(i)) = rvec2d(kn(i)) + rhs(i)
	    do j=1,3
	      if( i /= j .and. mass(i,j) /= 0. ) then
	        write(6,*) ie,kn(i),i,j,mass(i,j)
		stop 'error stop system_assemble: non diag elems /= 0'
	      end if
	    end do
	  end do
	else
         do i=1,3
          do j=1,3
            kk=ijp_ie(i,j,ie)			!COOGGU
            if(kk.gt.0) then
              c2coo(kk) = c2coo(kk) + mass(i,j)
              if (i2coo(kk).eq.j2coo(kk)) then
              index_coo(kk) = kn(i)
              end if    
            end if
          end do
          rvec2d(kn(i)) = rvec2d(kn(i)) + rhs(i)
         end do
	end if

	end

!******************************************************************

	subroutine system_assemble_3d(ie,l,nlv,kn,mass,rhs)

! assembles element matrix into system matrix

	use system_matrix
        use coo_matrix

	implicit none

	integer ie,l,nlv
	integer kn(3)
	double precision mass(-1:1,3,3)
	double precision rhs(3)

	integer i,j,kk

        do i=1,3
          do j=1,3
	    kk = loccoo3d(i,j,kn,l,ie)
            if(kk.gt.0) then
	       c3coo(kk-1) = c3coo(kk-1) + mass(-1,i,j)
	       c3coo(kk) = c3coo(kk) + mass(0,i,j)
	       c3coo(kk+1) = c3coo(kk+1) + mass(+1,i,j)
	    end if
          end do
	  kk = (nlv+2)*(kn(i)-1) + l + 1
          rvec3d(kk) = rvec3d(kk) + rhs(i) !DEB
        end do
	
	end

!******************************************************************

        subroutine system_adjust_z(n,z)

! copies solution back to z

	use system_matrix

        implicit none

	integer n
	double precision z(n)

        integer k

	!z = double precision(rvec2d)
        z = rvec2d

        end

!******************************************************************

        subroutine system_adjust_3d(n,nlvdi,nlv,z)

! copies solution back to z

	use system_matrix

        implicit none

	integer n,nlvdi,nlv
	double precision z(nlvdi,n)

	integer i,k,l

	i = 0

	do k=1,n
	  i = i + 1
	  do l=1,nlv
	    i = i + 1
	    z(l,k) = rvec3d(i) !DEB
	  end do
	  i = i + 1
	end do

        end

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine system_add_rhs(dt,n,array)

! adds right hand side to system array

	use system_matrix
        use shympi

        implicit none

        double precision dt
	integer n
        double precision array(n)

        integer k

#ifndef DEBUGON
        call shympi_exchange_and_sum_2D_nodes(rvec2d)
#endif

        do k=1,n
          rvec2d(k) = rvec2d(k) + dt * array(k)
        end do

        end

!******************************************************************

        subroutine system_adjust_matrix_3d

        use coo_matrix

        implicit none

        call coo_adjust

        end

!******************************************************************

!----------------------------------------------------------------------
        end module system
!----------------------------------------------------------------------
