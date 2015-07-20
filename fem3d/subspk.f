c
c $Id: subspk.f,v 1.4 2010-03-22 15:29:31 georg Exp $
c
c Sparskit solver routines
c
c revision log :
c
c 05.06.2009    ggu     some routines cleaned
c 22.03.2010    ggu     change of some parameters
c 29.03.2012    ggu     introduce zero and one as double (bug, was int before)
c
c*************************************************************************

      subroutine spk_init_system

! Initialize vector and matrix      

	use mod_system
	use basin

      implicit none
      include 'param.h'

      integer n

      integer icall_coo
      data icall_coo /0/
      save icall_coo

      do n=1,nkn
         rvec(n) = 0.
         raux(n) = 0.
      end do

      if (icall_coo.eq.0) then		! only first time
	 write(6,*) 'SOLVER: Sparskit'
         call coo_init(nel,nkn,mbw,nen3v,csrdim,nnzero,ijp,icoo,jcoo)
         print*, 'coo-matrix initialisation...'
         print*, 'Number of non-zeros: ',nnzero
         icall_coo=1
      end if

      do n=1,nnzero
         coo(n) = 0.
      end do

      end

c*************************************************************************

      subroutine spk_solve_system(n,z)

! Solver routine with Sparskit iterative methods.

	use mod_system

      implicit none

	integer n
	real z(n)	!first guess

      include 'param.h'
      !include 'nbasin.h'
      integer k

      real*8 csr(csrdim)
      integer icsr(n+1),jcsr(csrdim)
      integer iwork(2*csrdim)

      external bcg
      external bcgstab
      external tfqmr
      external fom
      external gmres
      external dqgmres
      external fgmres
      external dbcg

      integer  itermax !Maximum number of iterations
      parameter (itermax=1000)
      real*8   alu(csrdim*2), fpar(16), wilut(n+1)
      real*8   wksp(8*n)
      real*8   guess(n)
	integer nkn
      integer  ju(n), jlu(2*csrdim), iperm(2*n)
      integer  iw(n), jw(2*n), ipar(16)
      integer  ierr, lfil, iwk, itsol, itpre, mbloc
      double precision   droptol, permtol
      double precision   zero,one

      !ITPACK
      integer IPARIT(12),IWKSP(3*n), INW
      real*8  FPARIT(12),WKSPIT(6*n+4*itermax), INIU(n)

	nkn = n

!--------------------------------------------------
! CONVERSION AND SORTING
      call coocsr(nkn,nnzero,coo,icoo,jcoo,csr,jcsr,icsr)
      
      if( nnzero .gt. csrdim .or. nkn+1 .gt. 2*csrdim ) goto 99

      call csort (nkn,csr,jcsr,icsr,iwork,.true.)
!--------------------------------------------------

      lfil    = 5	! Number of largest elements, in absolute 
      			! value, of each row or column to keep. The
			! other ones are dropped.
      mbloc   = nkn
      droptol = 0.01d0
      permtol = 0.01d0
      iwk     = 2 * nnzero
      itpre   = 2	! either 2 or 3
      itsol   = 3	! 3 is probably best

      zero = 0.
      one = 1.

      if( itsol .gt. 9 ) stop 'error stop spk_solve_system: itsol'	!ggu

!--------------------------------------------------
! PRECONDITIONERS
      if (itpre .eq. 1) then
        call ilu0(nkn, csr, jcsr, icsr, alu, jlu, ju, iw, ierr)
      else if (itpre .eq. 2) then
        call ilut(nkn,csr,jcsr,icsr,lfil,droptol,alu,
     +            jlu,ju,iwk,wilut,jw,ierr)
      else if (itpre .eq. 3) then
        call ilutp(nkn,csr,jcsr,icsr,lfil,droptol,permtol,mbloc,alu,
     +             jlu,ju,iwk,wilut,jw,iperm,ierr)
      else if (itpre .eq. 4) then
        call ilud(nkn,csr,jcsr,icsr,zero,droptol,alu,
     +            jlu,ju,iwk,wilut,jw,ierr)
      else if (itpre .eq. 5) then
        call iludp(nkn,csr,jcsr,icsr,one,droptol,permtol,mbloc,alu,
     +             jlu,ju,iwk,wilut,jw,iperm,ierr)
      end if
!--------------------------------------------------

!--------------------------------------------------
! SOLVERS
      ipar(1) = 0      		! always 0 to start an iterative solver
      ipar(2) = 1  		! right preconditioning
      ipar(3) = 1      		! use convergence test scheme 1
      ipar(4) = 8*nkn		! size of the work array
      !ipar(5) = 10		! size of the Krylov subspace 
				!(used by GMRES and its variants)
      ipar(5) = 1
      ipar(6) = itermax		! use at most 100 matvec's
      !next values needed in order to run test kreis
      fpar(1) = 1.0E-8		! relative tolerance 1.0E-6
      fpar(2) = 1.0E-8		! absolute tolerance 1.0E-1

	do k=1,nkn
	  !guess(k) = 0.
	  guess(k) = z(k)	!ggu
	end do

       if (itsol .eq. 1) then
         call runrc(nkn,rvec,raux,ipar,fpar,wksp,guess,csr,jcsr,icsr,
     +     alu,jlu,ju,dbcg)
       else if (itsol .eq. 2) then
         call runrc(nkn,rvec,raux,ipar,fpar,wksp,guess,csr,jcsr,icsr,
     +     alu,jlu,ju,bcg)
       else if (itsol .eq. 3) then
         call runrc(nkn,rvec,raux,ipar,fpar,wksp,guess,csr,jcsr,icsr,
     +     alu,jlu,ju,bcgstab)
       else if (itsol .eq. 4) then
         call runrc(nkn,rvec,raux,ipar,fpar,wksp,guess,csr,jcsr,icsr,
     +     alu,jlu,ju,fom)
       else if (itsol .eq. 5) then
         call runrc(nkn,rvec,raux,ipar,fpar,wksp,guess,csr,jcsr,icsr,
     +     alu,jlu,ju,tfqmr)
       else if (itsol .eq. 6) then
         call runrc(nkn,rvec,raux,ipar,fpar,wksp,guess,csr,jcsr,icsr,
     +     alu,jlu,ju,gmres)
       else if (itsol .eq. 7) then
         call runrc(nkn,rvec,raux,ipar,fpar,wksp,guess,csr,jcsr,icsr,
     +     alu,jlu,ju,fgmres)
       else if (itsol .eq. 8) then
         call runrc(nkn,rvec,raux,ipar,fpar,wksp,guess,csr,jcsr,icsr,
     +     alu,jlu,ju,dqgmres)
       end if
       
!-----------------------------------------------------------------

      do k=1,nkn
        rvec(k) = raux(k)
      end do

	return
   99	continue
	write(6,*) nnzero,nkn+1,csrdim
	stop 'error stop spk_solve_system: dimension iwork'
	end

c*************************************************************************

