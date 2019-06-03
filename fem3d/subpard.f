
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
c 31.03.2009	ggu	call renamed to pard_*
c 23.03.2010	ggu	changed v6.1.1
c 10.07.2015	ggu	changed VERS_7_1_50
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 12.10.2015	ggu	changed VERS_7_3_4
c 15.12.2015	dbf	adjusted for new 3d framework
c 18.01.2019	ggu	changed VERS_7_5_55
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2019	ggu	changed VERS_7_5_61
c
c*************************************************************************

      subroutine pard_init_system

! Initialize vector and matrix      

	use mod_system
        use basin

      implicit none

      integer, save :: icall_coo = 0

      if (icall_coo.eq.0) then                ! only first time

	call coo_init_new 

        if( n2zero > n2max ) then
	      stop 'error stop spk_init_system: non zero 2d max'
        end if
        if( n3zero > n3max ) then
	      stop 'error stop spk_init_system: non zero 3d max'
        end if

        write(6,*) 'SOLVER: Pardiso'
        print*, 'coo-matrix initialisation...'
        print*, 'Number of non-zeros 2d:',n2zero,n2max
        print*, 'Number of non-zeros 3d:',n3zero,n3max

        icall_coo=1

      endif

      rvec2d = 0.
      raux2d = 0.
      rvec3d = 0.
      raux3d = 0.

      c2coo = 0.
      c3coo = 0.

      end

c*************************************************************************

	subroutine pard_solve_system(buse3d,nndim,n,z)

	use mod_system

	implicit none

	logical buse3d
        integer nndim           !dimension of non zeros in system
	integer n               !dimension of system(unknowns)
	real z(n)               !first guess

	integer k,i

        real*8, allocatable :: csr(:)
        real*8, allocatable :: rvec(:)
        real*8, allocatable :: raux(:)
        real*8, allocatable :: ddum(:)

        integer, allocatable :: icsr(:),jcsr(:)
        integer, allocatable :: iwork(:)
			      
	integer precision
	integer nth		!number of threads for pardiso
	logical bdirect		!iterative solver
	integer ngl,nnzero


	integer, save :: icall_coo = 0

	allocate(csr(nndim),icsr(n+1),jcsr(nndim),iwork(2*nndim))!DEB
        allocate(rvec(nndim),raux(nndim))!DEB

        ngl = n
	nth = nthpard
        rvec = 0.
        raux = 0.

	precision = iprec	!precision - see common.h

	bdirect = precision .eq. 0

!-----------------------------------------------------------------
! coo to csr conversion and sorting	
!-----------------------------------------------------------------
      
      if( buse3d ) then
           nnzero = n3zero
	   write(6,*)'nnzero',nnzero,'ngl',ngl
           call coocsr(ngl,nnzero,c3coo,i3coo,j3coo,csr,jcsr,icsr)!COOGGU
      else
	   nnzero = n2zero
	   write(6,*)'2D nnzero',nnzero,'n',n
	   call coocsr(ngl,nnzero,c2coo,i2coo,j2coo,csr,jcsr,icsr)!COOGGU
      endif
      if( nnzero .gt. nndim .or. ngl+1 .gt. 2*nndim ) goto 99

      !call csort (ngl,csr,jcsr,icsr,iwork,.true.)
      call csort (ngl,csr,jcsr,icsr,.true.)     !new calling modus spk2

!-----------------------------------------------------------------
! initialization
!-----------------------------------------------------------------

	if (bdirect .or. icall_coo.eq.0) then 
          if(icall_coo.eq.0) then
		  print*, 'Pardiso initialisation'
		  if( bdirect ) then
		  print*, 'Pardiso direct solution: ',precision
		  else
		  print*, 'Pardiso iterative solution: ',precision
		  end if
          endif
          call pardiso_solve(0,nth,ngl,precision,csr,icsr,jcsr,
     +          ddum,ddum)
	  icall_coo = 1
        end if

!-----------------------------------------------------------------
! solving	
!-----------------------------------------------------------------

        if( buse3d ) then
          rvec = rvec3d
          raux = raux3d
	else
          rvec = rvec2d
          raux = raux2d
        end if

        call pardiso_solve(1,nth,ngl,precision,csr,icsr,jcsr,
     .		rvec,raux)

	if( bdirect ) then
          call pardiso_solve(3,nth,ngl,precision,csr,icsr,jcsr,ddum,
     +         ddum)
	end if

        if( buse3d ) then
	       rvec3d(1:ngl) = raux(1:ngl)
	else
	       rvec2d(1:ngl) = raux(1:ngl)
        endif


        deallocate(csr,icsr,jcsr,iwork)
        deallocate(raux,rvec)

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

	return
   99	continue
	write(6,*) nnzero,ngl+1,nndim
	stop 'error stop pard_solve_system: dimension iwork'
	end

c*************************************************************************

      subroutine pardiso_solve(pcall,nth,nrow,precision,aa,iaa,jaa,b,x)

! Solve matrix with pardiso routines

      implicit none
!      external pardiso

! arguments

      integer pcall			!what to do
      integer nth			!number of threads
      integer nrow
      integer precision
      real*8 aa(*)
      integer iaa(nrow+1),jaa(*)
      real*8 b(nrow)			!right hand side
      real*8 x(nrow)			!result x

! Pardiso vars
      integer*8 pt(64)
      integer iparm(64),mtype
      integer maxfct,mnum,phase,nrhs,msglvl,error
      integer perm(nrow) !permutation vector or specifies elements used 
                         !for computing a partial solution
      data nrhs /1/, maxfct /1/, mnum /1/
      integer i
      
! Variables to save
      save pt,iparm
      
      logical pdefault
      
      ! Set the number of threads for the MKL library
      call mkl_set_num_threads(nth)

      ! Choose if default values or custom
      pdefault = .true.
      mtype = 11 ! real unsymmetric

      error = 0 ! initialize error flag
      msglvl = 0 ! print statistical information

      if (pcall.eq.0) then  !Initialization

        if( pdefault ) then
            ! Set default values for nonsymmetric matrix
            call pardisoinit(pt, 11, iparm)
        else
            ! Set the input parameters
            do i = 1, 64
               iparm(i) = 0
            end do
            iparm(1) = 1 ! no solver default
            iparm(2) = 3
            !iparm(2) = 2 ! fill-in reordering from METIS (suggested for symmetric)
            iparm(4) = precision	!use 0 for direct, else 31,61,91 etc..
            iparm(5) = 0 ! no user fill-in reducing permutation
            iparm(6) = 0 ! =0 solution on the first n compoments of x
            !iparm(7) = -1 ! number of iterative refinement steps (out)
            iparm(8) = 0 ! max numbers of iterative refinement steps
            iparm(10) = 13 ! perturbe the pivot elements with 1E-13
            iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
            iparm(12) = 0 ! Solve with transposed or conjugate transposed matrix A
            iparm(13) = 1 ! improved accuracy using nonsymmetric matchings
            !iparm(14) = 0 ! number of perturbed pivots (out)
            !iparm(15) = 0 ! Peak memory on symbolic factorization (out)
            !iparm(16) = 0 ! Permanent  memory on symbolic factorization (out)
            !iparm(17) = 0 ! As 15 added with solution (out)
            iparm(18) = -1 ! Report the number of nonzeros in the factor LU (in/out)
            iparm(19) = 0 ! Report Mflops that are necessary to factor the matrix A (use 0 to speedup)
            !iparm(20) = 0 ! Numbers of CG Iterations
            iparm(21) = 1 ! pivoting
            iparm(24) = 0 ! Parallel factorization control (use 1 with many threads > 8)
            iparm(25) = 0 ! Parallel forward/backward solve control.
            iparm(27) = 0 ! Matrix checker
            iparm(28) = 0 ! Single or double precision of PARDISO (0 = double)
            iparm(31) = 0 ! Partial solve and computing selected components of the solution vectors
            iparm(34) = 0 ! Optimal number of OpenMP threads for conditional numerical reproducibility (CNR) mode
            iparm(35) = 0 ! 1-based/0-based input data indexing (0 = fortran style)
            iparm(36) = 0 ! Schur complement matrix computation control
            iparm(56) = 0 ! Diagonal and pivoting control
            iparm(60) = 0 ! PARDISO mode (2 holds the matrix factors in files)
            
            !Initiliaze the internal solver memory pointer.
            do i = 1, 64
               pt(i) = 0
            end do
            
        end if

         phase = 11 
         call pardiso (pt,maxfct,mnum,mtype,phase,nrow,aa,iaa,jaa,
     &                 perm,nrhs,iparm,msglvl,b,x,error)
	 return
     
      elseif (pcall.eq.1) then

         !.. Numerical factorization, Solve, Iterative refinement
         phase = 23 
         call pardiso (pt,maxfct,mnum,mtype,phase,nrow,aa,iaa,jaa,
     &                 perm,nrhs,iparm,msglvl,b,x,error)
         !print*, 'Number of CG iterations: ',iparm(20)
         !print*, 'Number of perturbed pivots: ',iparm(7)
         !print*, 'Mflops for LU factorisation: ',iparm(19)
         !print*, 'Number of refinement steps: ',iparm(14)
         !print*, 'Number of non-zero: ',iparm(18)
         if (error .ne. 0) then
            write(*,*) '1 The following ERROR was detected: ', error
            stop
         end if
	 return

      elseif (pcall.eq.2) then !Solve system

         !.. Solve and iterative refinement
         phase = 33
         call pardiso (pt,maxfct,mnum,mtype,phase,nrow,aa,iaa,jaa,
     &                 perm,nrhs,iparm,msglvl,b,x,error)
         !write(*,*) 'Solve completed ... '
         if (error .ne. 0) then
            write(*,*) '2 The following ERROR was detected: ', error
            stop
         end if
	 return

      elseif (pcall.eq.3) then ! Release memory

         phase = -1 ! release internal memory
         call pardiso (pt,maxfct,mnum,mtype,phase,nrow,aa,iaa,jaa,
     &                 perm,nrhs,iparm,msglvl,b,x,error)
         !write(*,*) 'Release completed ... '
	 return

      else

	 stop 'Pardiso call: Error.'

      end if

      end

c*************************************************************************
