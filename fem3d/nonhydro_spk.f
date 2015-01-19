
	subroutine nonhydro_solve_matrix

! sparskit solver 
	
	implicit none
	include 'param.h'
        include 'nohydlinks.h'
	include 'nlevel.h'
	include 'levels.h'
	include 'basin.h'
        integer k

        real*8 csr(csrdimnh)
        integer icsr(matdimmax+1),jcsr(csrdimnh)
        integer iwork(2*csrdimnh)

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
        real*8   alu(csrdimnh*2), fpar(16), wilut(matdimmax+1)
        real*8   wksp(8*matdimmax)
        real*8   guess(matdimmax)
        integer  ju(matdimmax), jlu(2*csrdimnh), iperm(2*matdimmax)
        integer  iw(matdimmax), jw(2*matdimmax), ipar(16)
        integer  ierr, lfil, iwk, itsol, itpre, mbloc
        double precision   droptol, permtol
        double precision   zero,one
	integer l,lmax,nn,nnn
	include 'nohyd.h'
	integer i
        !ITPACK
        integer IPARIT(12),IWKSP(3*nkndim), INW
        real*8  FPARIT(12),WKSPIT(6*nkndim+4*itermax), INIU(nkndim)

!--------------------------------------------------
! CONVERSION AND SORTING
      call coocsr(matdimmax,nnzeronh,conh,ioii1,iojj1,csr,jcsr,icsr)
      !write(653,*)(csr(i),i=1,nnzeronh) 
      if( nnzeronh.gt.csrdimnh.or. matdimmax+1 .gt. 2*csrdimnh ) goto 99
 
      call csort (matdimmax,csr,jcsr,icsr,iwork,.true.)

!--------------------------------------------------

      !lfil    = 5	! Number of largest elements, in absolute 
      lfil    = ngr + 3	! Number of largest elements, in absolute 
      			! value, of each row or column to keep. The
			! other ones are dropped.
      !mbloc   = nkn*nlv
      mbloc   = matdimmax
      droptol = 0.01d0
      permtol = 0.01d0
      iwk     = 2 * nnzeronh
      itpre   = 2	! either 2 or 3
      itsol   = 3	! 3 is probably best

      zero = 0.
      one = 1.

      if( itsol .gt. 9 ) stop 'error stop spk_solve_system: itsol'	!ggu

!--------------------------------------------------
! PRECONDITIONERS
      if (itpre .eq. 1) then
        call ilu0(matdimmax, csr, jcsr, icsr, alu, jlu, ju, iw, ierr)
      else if (itpre .eq. 2) then
        call ilut(matdimmax,csr,jcsr,icsr,lfil,droptol,alu,
     +            jlu,ju,iwk,wilut,jw,ierr)
       !	write(6,*)'ilut ierr',ierr
      else if (itpre .eq. 3) then
        call ilutp(matdimmax,csr,jcsr,icsr,lfil,droptol,permtol,mbloc,
     +		   alu,jlu,ju,iwk,wilut,jw,iperm,ierr)
      else if (itpre .eq. 4) then
        call ilud(matdimmax,csr,jcsr,icsr,zero,droptol,alu,
     +            jlu,ju,iwk,wilut,jw,ierr)
      else if (itpre .eq. 5) then
        call iludp(matdimmax,csr,jcsr,icsr,one,droptol,permtol,mbloc,
     +		alu,jlu,ju,iwk,wilut,jw,iperm,ierr)
      end if
!--------------------------------------------------

!--------------------------------------------------
! SOLVERS
      ipar(1) = 0      		! always 0 to start an iterative solver
      ipar(2) = 1  		! right preconditioning
      ipar(3) = 1      		! use convergence test scheme 1
      !ipar(4) = 8*nkndim*nlvdim	! size of the work array
      ipar(4) = 8*matdimmax	! size of the work array
      !ipar(5) = 10		! size of the Krylov subspace 
				!(used by GMRES and its variants)
      ipar(5) = 1
      ipar(6) = itermax		! use at most 100 matvec's
      !next values needed in order to run test kreis
      fpar(1) = 1.0E-8		! relative tolerance 1.0E-6
      fpar(2) = 1.0E-8		! absolute tolerance 1.0E-1
      !fpar(1) = 1.0E-6		! relative tolerance 1.0E-6
      !fpar(2) = 1.0E-1		! absolute tolerance 1.0E-1

      nn=0
	do k=1,nkn
         lmax=ilhkv(k)
         do l=1,lmax
	  nn=nn+1
	  guess(nn) = qpnv(l,k)	!deb
	  !guess(nn) = 0.0d0 !deb
	 enddo
	end do
     
       if (itsol .eq. 1) then
         call runrc(matdimmax,rvecnh,rauxnh,ipar,fpar,wksp,guess,csr,
     +	   jcsr,icsr,alu,jlu,ju,dbcg)
       else if (itsol .eq. 2) then
         call runrc(matdimmax,rvecnh,rauxnh,ipar,fpar,wksp,guess,csr,
     +     jcsr,icsr,alu,jlu,ju,bcg)
       else if (itsol .eq. 3) then
         call runrc(matdimmax,rvecnh,rauxnh,ipar,fpar,wksp,guess,csr,
     +     jcsr,icsr,alu,jlu,ju,bcgstab)
         !write(6,*)'runrc'
       else if (itsol .eq. 4) then
         call runrc(matdimmax,rvecnh,rauxnh,ipar,fpar,wksp,guess,csr,
     +     jcsr,icsr,alu,jlu,ju,fom)
       else if (itsol .eq. 5) then
         call runrc(matdimmax,rvecnh,rauxnh,ipar,fpar,wksp,guess,csr,
     +     jcsr,icsr,alu,jlu,ju,tfqmr)
       else if (itsol .eq. 6) then
         call runrc(matdimmax,rvecnh,rauxnh,ipar,fpar,wksp,guess,csr,
     +     jcsr,icsr,alu,jlu,ju,gmres)
       else if (itsol .eq. 7) then
         call runrc(matdimmax,rvecnh,rauxnh,ipar,fpar,wksp,guess,csr,
     +     jcsr,icsr,alu,jlu,ju,fgmres)
       else if (itsol .eq. 8) then
         call runrc(matdimmax,rvecnh,rauxnh,ipar,fpar,wksp,guess,csr,
     +     jcsr,icsr,alu,jlu,ju,dqgmres)
       end if
       
!-----------------------------------------------------------------

      nn=0
      do k=1,nkn
         lmax=ilhkv(k)
         do l=1,lmax
	  nn=nn+1
          rvecnh(nn) = rauxnh(nn)
c	 if(rvecnh(nn).ne.0.)write(6,*)'oaoao',nn,rvecnh(nn),rauxnh(nn)
	 enddo
      end do

	return
   99	continue
	write(6,*) nnzeronh,nkn+1,csrdimnh
	stop 'error stop spk_solve_system: dimension iwork'
	
	end
