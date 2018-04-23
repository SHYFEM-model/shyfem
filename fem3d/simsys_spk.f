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

        subroutine system_initialize

! allocates data structure

	use mod_system
	use levels
	use basin
	use shympi

        implicit none

        write(6,*) '----------------------------------------'
        write(6,*) 'initializing matrix inversion routines'
        write(6,*) 'using Sparskit routines'
        write(6,*) '----------------------------------------'

        call mod_system_init(nkn,nel,ngr,mbw,nlv,l_matrix)
	call mod_system_insert_elem_index(nel,nen3v,l_matrix)

	call mod_system_set_local

! next we first have to set bstsexpl - not yet done... !FIXME

	if( bmpi ) then	!only needed if not explicit !FIXME
	 if( .not. bsysexpl ) then	!only needed if not explicit !FIXME
          call mod_system_init(nkn_global,nel_global,ngr_global
     +				,mbw,nlv,g_matrix)
	  stop 'error stop system_initialize: not ready'
	  !call mod_system_insert_elem_index(nel,nen3v,g_matrix)
	 end if
	end if

        end

!******************************************************************

	subroutine system_init

! on first call initializes pointers - sets arrays to zero
!
! must be called before every assembly of matrix

	use mod_system
	use shympi

	implicit none

	if( bsysexpl ) then
	  call mod_system_set_local
	  a_matrix%rvec2d = 0.
	  a_matrix%raux2d = 0.
	else
	  if( bmpi ) then
	    call mod_system_set_global
	    call spk_init_system
	  end if
	  call mod_system_set_local
	  call spk_init_system
	end if

	end

!******************************************************************

	subroutine system_set_explicit

	use mod_system

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
!	use mod_system
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

	subroutine system_solve_global_z(n,z)

! this solves the global system

	use mod_system
	use shympi

	implicit none

	integer n
	real z(n)

	stop 'error stop system_solve_global_z: not ready'

	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine system_solve_z(n,z)

! solves system - z is used for initial guess

	use mod_system
	use shympi

	implicit none

	integer n
	real z(n)

	integer n2max
	type(smatrix), pointer :: m

	m => a_matrix

	if( bsysexpl ) then
	  !write(6,*) 'solving explicitly...'
          call shympi_exchange_and_sum_2d_nodes(m%rvec2d)
          call shympi_exchange_and_sum_2d_nodes(m%raux2d)
          !call shympi_comment('shympi_elem: exchange rvec2d, raux2d')
	  m%rvec2d = m%rvec2d / m%raux2d	!GGUEXPL
	else if( bmpi ) then
	  call system_solve_global_z(n,z)
	else			!solve directly locally
	  n2max = m%n2max
	  call spk_solve_system(.false.,n2max,n,z)
	end if

	end

!******************************************************************

	subroutine system_solve_3d(n,nlvdi,nlv,z)

! solves system - z is used for initial guess

	use mod_system

	implicit none

	integer n,nlvdi,nlv
	real z(nlvdi,n)

	integer i,k,l
	integer n3max
	real p((nlv+2)*n)
	integer nn !DEB

	i = 0
	n3max = a_matrix%n3max

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

	use mod_system

	implicit none

	integer ie
	integer kn(3)
	real mass(3,3)
	real rhs(3)

	integer i,j,kk
	type(smatrix), pointer :: m

	m => a_matrix

	if( bsysexpl ) then
          do i=1,3
            m%raux2d(kn(i)) = m%raux2d(kn(i)) + mass(i,i)	!GGUEXPL
            m%rvec2d(kn(i)) = m%rvec2d(kn(i)) + rhs(i)
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
            kk=m%ijp_ie(i,j,ie)			!COOGGU
            if(kk.gt.0) m%c2coo(kk) = m%c2coo(kk) + mass(i,j)
          end do
          m%rvec2d(kn(i)) = m%rvec2d(kn(i)) + rhs(i)
         end do
	end if

	end

!******************************************************************

	subroutine system_assemble_3d(ie,l,nlv,kn,mass,rhs)

! assembles element matrix into system matrix

	use mod_system

	implicit none

	integer ie,l,nlv
	integer kn(3)
	real mass(-1:1,3,3)
	real rhs(3)

	integer i,j,kk
	type(smatrix), pointer :: m

	integer loccoo3d
	external loccoo3d

	m => a_matrix

        do i=1,3
          do j=1,3
	    kk = loccoo3d(i,j,kn,l,ie)
            if(kk.gt.0) then
	       m%c3coo(kk-1) = m%c3coo(kk-1) + mass(-1,i,j)
	       m%c3coo(kk) = m%c3coo(kk) + mass(0,i,j)
	       m%c3coo(kk+1) = m%c3coo(kk+1) + mass(+1,i,j)
	    end if
          end do
	  kk = (nlv+2)*(kn(i)-1) + l + 1
          m%rvec3d(kk) = m%rvec3d(kk) + rhs(i) !DEB
        end do
	
	end

!******************************************************************

        subroutine system_adjust_z(n,z)

! copies solution back to z

	use mod_system

        implicit none

	integer n
	real z(n)

        integer k

	z = real(a_matrix%rvec2d)

        end

!******************************************************************

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
	    z(l,k) = real(a_matrix%rvec3d(i)) !DEB
	  end do
	  i = i + 1
	end do

        end

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine system_add_rhs(dt,n,array)

! adds right hand side to system array

	use mod_system

        implicit none

        real dt
	integer n
        real array(n)

        integer k

        do k=1,n
          a_matrix%rvec2d(k) = a_matrix%rvec2d(k) + dt * array(k)
        end do

        end

!******************************************************************

        subroutine system_adjust_matrix_3d

        implicit none

        call coo_adjust

        end

!******************************************************************

