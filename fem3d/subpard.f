
c*************************************************************************

      subroutine pard_init_system

! Initialize vector and matrix      

	use mod_system
        use basin

      implicit none
      include 'param.h'

      integer n

      integer icall
      data icall /0/
      save icall

      do n=1,nkn
         rvec(n) = 0.
      end do

      if (icall.eq.0) then		! only first time
	 write(6,*) 'SOLVER: Pardiso'
         call coo_init(nel,nkn,mbw,nen3v,csrdim,nnzero,ijp,icoo,jcoo)
         print*, 'coo-matrix initialisation...'
         print*, 'Number of non-zeros: ',nnzero
         icall=1
      end if

      do n=1,nnzero
         coo(n) = 0.
      end do

      end

c*************************************************************************

	subroutine pard_solve_system(n)

	use mod_system

	implicit none

        integer n !dimension of x and b

        include 'param.h'
	integer k

        real*8 csr(csrdim)
        integer icsr(n+1),jcsr(csrdim)
	integer iwork(2*csrdim)		!aux for sorting routine
        real*8 ddum(n)
        !integer indu(n),iwk(n+1) !clean-csr vectors

	integer precision
	logical bdirect		!iterative solver

        integer icall
        data icall /0/
        save icall

        integer nkn

        nkn = n

	precision = iprec	!precision - see common.h

	bdirect = precision .eq. 0

!-----------------------------------------------------------------
! coo to csr conversion and sorting	
!-----------------------------------------------------------------

	call coocsr(nkn,nnzero,coo,icoo,jcoo,csr,jcsr,icsr)
	!nnzero = icsr(nkn+1)-1		!already in common.h
	if( nnzero .gt. csrdim .or. nkn+1 .gt. 2*csrdim ) goto 99

        call csort (nkn,csr,jcsr,icsr,iwork,.true.)
	!call clncsr(3,1,nkn,csr,jcsr,icsr,indu,iwk)

!-----------------------------------------------------------------
! initialization
!-----------------------------------------------------------------

        if ( bdirect .or. icall .eq. 0 ) then
           if( icall .eq. 0 ) then
	     print*, 'Pardiso initialisation'
	     if( bdirect ) then
	       print*, 'Pardiso direct solution: ',precision
	     else
	       print*, 'Pardiso iterative solution: ',precision
	     end if
	   end if
           call pardiso_solve(0,nkn,precision,csr,icsr,jcsr,ddum,ddum)
           icall=1
        end if
	
!-----------------------------------------------------------------
! solving	
!-----------------------------------------------------------------

        call pardiso_solve(1,nkn,precision,csr,icsr,jcsr,rvec,raux)

	if( bdirect ) then
          call pardiso_solve(3,nkn,precision,csr,icsr,jcsr,ddum,ddum)
	end if

	do k=1,nkn
	  rvec(k) = raux(k)
	end do

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

	return
   99	continue
	write(6,*) nnzero,nkn+1,csrdim
	stop 'error stop pard_solve_system: dimension iwork'
	end

c*************************************************************************

      subroutine pardiso_solve(pcall,nrow,precision,aa,iaa,jaa,b,x)

! Solve matrix with pardiso routines

!$      use omp_lib
      implicit none
!      external pardiso

! arguments

      integer pcall			!what to do
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

      mtype = 11 ! real unsymmetric
      error = 0 ! initialize error flag
      msglvl = 0 ! print statistical information

      if (pcall.eq.0) then  !Initialization

         ! Set the input parameters if not using the default values
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
         
         !.. Initiliaze the internal solver memory pointer. This is only
         !   necessary for the FIRST call of the PARDISO solver.
         do i = 1, 64
            pt(i) = 0
         end do

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

