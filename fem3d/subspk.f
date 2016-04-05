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
c 15.12.2015    ggu&deb adapted to new 3d framework
c 13.01.2016    ggu&ivn bug in allocation of rvec and raux (n instead nndim)
c 04.04.2016    ggu	make some big arrays allocatable and not on stack
c
c*************************************************************************

	subroutine spk_init_system

! Initialize vector and matrix      

	use mod_system
	use basin

	implicit none

	integer n,i
	integer nozero
	integer nozero_max

	integer, save :: icall_coo = 0

	if (icall_coo.eq.0) then		! only first time

	  call coo_init_new

	  if( n2zero > n2max ) then
	    stop 'error stop spk_init_system: non zero 2d max'
	  end if
	  if( n3zero > n3max ) then
	    stop 'error stop spk_init_system: non zero 3d max'
	  end if

	  write(6,*) 'SOLVER: Sparskit'
          print*, 'coo-matrix initialisation...'
          print*, 'Number of non-zeros 2d: ',n2zero,n2max
          print*, 'Number of non-zeros 3d: ',n3zero,n3max

          icall_coo=1
	end if

	rvec2d = 0.
	raux2d = 0.
	rvec3d = 0.
	raux3d = 0.
	c2coo = 0.
	c3coo = 0.

      end

c*************************************************************************

      subroutine spk_solve_system(buse3d,nndim,n,z)

! Solver routine with Sparskit iterative methods.

	use mod_system

      implicit none

	logical buse3d
	integer nndim		!dimension of non zeros in system
	integer n		!dimension of system (unknowns)
	real z(n)		!first guess

      integer k,i

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

      real*8   fpar(16)
      real*8   wilut(n+1)
      real*8   guess(n)
      !real*8   wksp(8*n)
      !real*8   alu(nndim*2)
      !integer  jlu(2*nndim)
      real*8, allocatable ::   wksp(:)
      real*8, allocatable ::   alu(:)
      integer, allocatable ::  jlu(:)
      integer  ju(n), iperm(2*n)
      integer  iw(n), jw(2*n)

      integer  ipar(16)
      integer ngl,nnzero
      integer  ierr, lfil, iwk, itsol, itpre, mbloc

      double precision   droptol, permtol
      double precision   zero,one

      !ITPACK
      integer IPARIT(12), INW
      !integer IWKSP(3*n)
      integer, allocatable :: IWKSP(:)
      real*8  FPARIT(12)
      real*8  INIU(n)
      !real*8  WKSPIT(6*n+4*itermax)
      real*8, allocatable ::  WKSPIT(:)

      real*8, allocatable :: csr(:)
      real*8, allocatable :: rvec(:)
      real*8, allocatable :: raux(:)

      integer, allocatable :: icsr(:),jcsr(:)
      integer, allocatable :: iwork(:)
	integer ii,nn

	allocate(csr(nndim),icsr(n+1),jcsr(nndim),iwork(2*nndim))
	allocate(rvec(n),raux(n))
	allocate(IWKSP(3*n),WKSPIT(6*n+4*itermax))
	allocate(wksp(8*n),alu(nndim*2),jlu(2*nndim))

	ngl = n

	ipar = 0	!ggu
	fpar = 0.

	rvec = 0.
	raux = 0.

!--------------------------------------------------
! CONVERSION AND SORTING
!--------------------------------------------------

	if( buse3d ) then
	  nnzero = n3zero
          call coocsr(ngl,nnzero,c3coo,i3coo,j3coo,csr,jcsr,icsr)
          !write(6,*)'3D nnzero',nnzero,ngl
	  !call coo_show(ngl,nnzero,i3coo,j3coo,c3coo)
	  !call csr_show(ngl,nnzero,icsr,jcsr,csr)
	  !call coo_print(ngl,nnzero,i3coo,j3coo,c3coo,rvec3d)
	else
	  nnzero = n2zero
          call coocsr(ngl,nnzero,c2coo,i2coo,j2coo,csr,jcsr,icsr)
          !write(6,*)'2D nnzero',nnzero,ngl
	  !call coo_show(ngl,nnzero,i2coo,j2coo,c2coo)
	  !call csr_show(ngl,nnzero,icsr,jcsr,csr)
	  !call coo_print(ngl,nnzero,i2coo,j2coo,c2coo,rvec)
	end if
      
      if( nnzero .gt. nndim .or. ngl+1 .gt. 2*nndim ) goto 99

      call csort (ngl,csr,jcsr,icsr,iwork,.true.)

!--------------------------------------------------

      lfil    = 5	! Number of largest elements, in absolute 
      			! value, of each row or column to keep. The
			! other ones are dropped.
      mbloc   = ngl
      droptol = 0.01d0
      permtol = 0.01d0
      iwk     = 2 * nnzero
      itpre   = 2	! either 2 or 3
      itsol   = 3	! 3 is probably best

      !itpre   = 3	! either 2 or 3
      !itsol   = 2	! 3 is probably best

      zero = 0.
      one = 1.

      if( itsol .gt. 9 ) stop 'error stop spk_solve_system: itsol'	!ggu

!--------------------------------------------------
! PRECONDITIONERS
!--------------------------------------------------
      if (itpre .eq. 1) then
        call ilu0(ngl, csr, jcsr, icsr, alu, jlu, ju, iw, ierr)
      else if (itpre .eq. 2) then
        call ilut(ngl,csr,jcsr,icsr,lfil,droptol,alu,
     +            jlu,ju,iwk,wilut,jw,ierr)
      else if (itpre .eq. 3) then
        call ilutp(ngl,csr,jcsr,icsr,lfil,droptol,permtol,mbloc,alu,
     +             jlu,ju,iwk,wilut,jw,iperm,ierr)
      else if (itpre .eq. 4) then
        call ilud(ngl,csr,jcsr,icsr,zero,droptol,alu,
     +            jlu,ju,iwk,wilut,jw,ierr)
      else if (itpre .eq. 5) then
        call iludp(ngl,csr,jcsr,icsr,one,droptol,permtol,mbloc,alu,
     +             jlu,ju,iwk,wilut,jw,iperm,ierr)
      end if
!--------------------------------------------------

!--------------------------------------------------
! SOLVERS
!--------------------------------------------------
      ipar(1) = 0      		! always 0 to start an iterative solver
      ipar(2) = 1  		! right preconditioning
      ipar(3) = 1      		! use convergence test scheme 1
      ipar(4) = 8*ngl		! size of the work array
      !ipar(5) = 10		! size of the Krylov subspace 
				!(used by GMRES and its variants)
      ipar(5) = 1
      ipar(6) = itermax		! use at most 100 matvec's
      !next values needed in order to run test kreis
      fpar(1) = 1.0E-8		! relative tolerance 1.0E-6
      fpar(2) = 1.0E-8		! absolute tolerance 1.0E-1

	guess = z
	!guess = 0.
	
      if( buse3d ) then
	rvec = rvec3d
	raux = raux3d
      else
	rvec = rvec2d
	raux = raux2d
      end if


       if (itsol .eq. 1) then
         call runrc(ngl,rvec,raux,ipar,fpar,wksp,guess,csr,jcsr,icsr,
     +     alu,jlu,ju,dbcg)
       else if (itsol .eq. 2) then
         call runrc(ngl,rvec,raux,ipar,fpar,wksp,guess,csr,jcsr,icsr,
     +     alu,jlu,ju,bcg)
       else if (itsol .eq. 3) then
         call runrc(ngl,rvec,raux,ipar,fpar,wksp,guess,csr,jcsr,icsr,
     +     alu,jlu,ju,bcgstab)
       else if (itsol .eq. 4) then
         call runrc(ngl,rvec,raux,ipar,fpar,wksp,guess,csr,jcsr,icsr,
     +     alu,jlu,ju,fom)
       else if (itsol .eq. 5) then
         call runrc(ngl,rvec,raux,ipar,fpar,wksp,guess,csr,jcsr,icsr,
     +     alu,jlu,ju,tfqmr)
       else if (itsol .eq. 6) then
         call runrc(ngl,rvec,raux,ipar,fpar,wksp,guess,csr,jcsr,icsr,
     +     alu,jlu,ju,gmres)
       else if (itsol .eq. 7) then
         call runrc(ngl,rvec,raux,ipar,fpar,wksp,guess,csr,jcsr,icsr,
     +     alu,jlu,ju,fgmres)
       else if (itsol .eq. 8) then
         call runrc(ngl,rvec,raux,ipar,fpar,wksp,guess,csr,jcsr,icsr,
     +     alu,jlu,ju,dqgmres)
       end if
       
!-----------------------------------------------------------------

        if( buse3d ) then
	  rvec3d(1:ngl) = raux(1:ngl)
	else
	  rvec2d(1:ngl) = raux(1:ngl)
	endif

	deallocate(csr,icsr,jcsr,iwork)
	deallocate(raux,rvec)
	deallocate(IWKSP,WKSPIT)
	deallocate(wksp,alu,jlu)

	return
   99	continue
	write(6,*) nnzero,ngl+1,nndim
	stop 'error stop spk_solve_system: dimension iwork'
	end

c*************************************************************************

