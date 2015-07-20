
c*************************************************************************

      subroutine pard_init_system

! Initialize vector and matrix      

	use mod_system

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

	subroutine pard_solve_system

	use mod_system

	implicit none
        include 'param.h'
	integer k

        real*8 csr(csrdim)
        integer icsr(nkn+1),jcsr(csrdim)
	integer iwork(2*csrdim)		!aux for sorting routine
        real*8 ddum
        !integer indu(nkn),iwk(nkn+1) !clean-csr vectors

	integer precision
	logical bdirect		!iterative solver

        integer icall
        data icall /0/
        save icall

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

      use omp_lib
      implicit none
!      external pardiso

! arguments

      integer pcall			!what to do
      integer nrow
      integer precision
      real*8 aa(*)
      integer iaa(nrow+1),jaa(*)
      real*8 b(*)			!right hand side
      real*8 x(*)			!result x

! Pardiso vars
      integer*8 pt(64)
      integer iparm(64),mtype
      integer maxfct,mnum,phase,nrhs,msglvl,error
      real*8 ddum  !Double prec dummy
      integer idum !Integer dummy
      data nrhs /1/, maxfct /1/, mnum /1/

      integer i

! Variables to save
      save pt,iparm

      mtype = 11 ! real unsymmetric
      error = 0 ! initialize error flag
      msglvl = 0 ! print statistical information

      if (pcall.eq.0) then  !Initialization

         do i = 1, 64
            iparm(i) = 0
         end do
         iparm(1) = 1 ! no solver default
	 iparm(2) = 0
	 if( precision .gt. 0 ) iparm(2) = 2
         !iparm(2) = 2 ! fill-in reordering from METIS (suggested for symmetric)
         !iparm(3) = 2 ! numbers of processors
	 !call mkl_set_num_threads(1)
	 !call omp_set_num_threads(1)
	 iparm(3)= omp_get_max_threads() !OMP_NUM_THREADS envirom. var.
         iparm(4) = precision	!use 0 for direct, else 31,61,91 etc..
         iparm(5) = 0 ! no user fill-in reducing permutation
         iparm(6) = 0 ! =0 solution on the first n compoments of x
         iparm(7) = -1 ! number of iterative refinement steps
         iparm(8) = 0 ! max numbers of iterative refinement steps
         iparm(9) = 0 ! not in use
         iparm(10) = 13 ! perturbe the pivot elements with 1E-13
         iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
         iparm(12) = 0 ! not in use
         iparm(13) = 1 ! improved accuracy using nonsymmetric matchings
         iparm(14) = 0 ! Output: number of perturbed pivots
         iparm(15) = 0 ! not in use
         iparm(16) = 0 ! not in use
         iparm(17) = 0 ! not in use
         iparm(18) = -1 ! Output: number of nonzeros in the factor LU
         iparm(19) = -1 ! Output: Mflops for LU factorization
         iparm(20) = 0 ! Output: Numbers of CG Iterations
         iparm(21) = 1 ! pivoting
         !iparm(23) = 1 ! 
         !iparm(24) = 1 ! 

         !.. Initiliaze the internal solver memory pointer. This is only
         !   necessary for the FIRST call of the PARDISO solver.
         do i = 1, 64
            pt(i) = 0
         end do

         phase = 11 
         call pardiso (pt,maxfct,mnum,mtype,phase,nrow,aa,iaa,jaa,
     &                 idum,nrhs,iparm,msglvl,b,x,error)
	 return
     
      elseif (pcall.eq.1) then

         !.. Numerical factorization, Solve, Iterative refinement
         phase = 23 
         call pardiso (pt,maxfct,mnum,mtype,phase,nrow,aa,iaa,jaa,
     &                 idum,nrhs,iparm,msglvl,b,x,error)
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
     &                 idum,nrhs,iparm,msglvl,b,x,error)
         !write(*,*) 'Solve completed ... '
         if (error .ne. 0) then
            write(*,*) '2 The following ERROR was detected: ', error
            stop
         end if
	 return

      elseif (pcall.eq.3) then ! Release memory

         phase = -1 ! release internal memory
         call pardiso (pt,maxfct,mnum,mtype,phase,nrow,ddum,idum,idum,
     &                 idum,nrhs,iparm,msglvl,ddum,ddum,error)
         !write(*,*) 'Release completed ... '
	 return

      else

	 stop 'Pardiso call: Error.'

      end if

      end

c*************************************************************************

