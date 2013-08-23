********************************************************************
       subroutine nonhydro_solve_matrix_pard

	implicit none
	include 'common.h'
        include 'nohydlinks.h'
        integer nlvdi,nlv
	 common /level/ nlvdi,nlv
         integer ilhv(1), ilhkv(1)
	 common /ilhv/ilhv, /ilhkv/ilhkv
         real*8 csr(csrdimnh)
        integer icsr(matdimmax+1),jcsr(csrdimnh)
        integer iwork(2*csrdimnh)
        !real*8 ddum(matdimmax)
        real*8 ddum
        integer precision
        logical bdirect
        integer icall
        data icall /0/
        save icall
	integer nn,l,k,lmax

	precision = iprec       !precision - see common.h

	bdirect = precision .eq. 0

!-----------------------------------------------------------------
! coo to csr conversion and sorting     
!-----------------------------------------------------------------
       call coocsr(matdimmax,nnzeronh,conh,ioii1,iojj1,csr,jcsr,icsr)

      if( nnzeronh.gt.csrdimnh.or. matdimmax+1 .gt. 2*csrdimnh ) goto 99

        call csort (matdimmax,csr,jcsr,icsr,iwork,.true.)


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
           !call pardiso_solve(0,nkn,precision,csr,icsr,jcsr,ddum,ddum)
      call pardiso_solve(0,matdimmax,precision,csr,icsr,jcsr,ddum,ddum)
           icall=1
        end if
	
!-----------------------------------------------------------------
! solving	
!-----------------------------------------------------------------

      call pardiso_solve(1,matdimmax,precision,csr,icsr,jcsr,rvecnh,
     +		rauxnh)

	if( bdirect ) then
      call pardiso_solve(3,matdimmax,precision,csr,icsr,jcsr,ddum,ddum)
	end if

      nn=0
      do k=1,nkn
         lmax=ilhkv(k)
         do l=1,lmax
	  nn=nn+1
          rvecnh(nn) = rauxnh(nn)
	 !if(rvecnh(nn).ne.0.)write(6,*)'oaoao',nn,rvecnh(nn),rauxnh(nn)
	 enddo
      end do

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

	return
   99	continue
	write(6,*) nnzeronh,matdimmax+1,csrdimnh
	stop 'error stop pard_solve_system: dimension iwork'
	end
c********************************************************************

