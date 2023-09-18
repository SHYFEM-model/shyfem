!
! $Id: subssv.f,v 1.4 1999/06/22 14:47:52 georg Exp $
!
! utility routines for operations on vectors
!
! contents :
!
! subroutine mima(xx,n,xmin,xmax)	min/max of vector
!
! subroutine mimari(xx,n,xmin,xmax,imin,imax,rnull)
!                                       min/max of vector with null value
!
! subroutine mima2i(text,n,a1,a2)	writes min/max of 2 arrays to terminal
!
! subroutine addvec(a1,a2,n)		adds two vectors
! subroutine mulvec(a1,a2,n)		multiplies two vectors
! subroutine zervec(a1,n)		initializes vector
!
! revision log :
!
! 26.08.1998    ggu	routines mimari transferred from newbcl0
! 31.05.1999    ggu	new comodity routine mima2i
!
!*******************************************
!
!----------------------------------------------------------------------
        module vec_util
!----------------------------------------------------------------------
        contains
!----------------------------------------------------------------------

	subroutine mima(xx,n,xmin,xmax)
!
! computes min/max of vector
!
! xx		vector
! n		dimension of vector
! xmin,xmax	min/max value in vector
!
        implicit none
!
        integer n,i
        double precision xx(n)
        double precision xmin,xmax,x
!
	xmax=xx(1)
	xmin=xmax
!
	do i=1,n
          x=xx(i)
          if(x.gt.xmax) xmax=x
          if(x.lt.xmin) xmin=x
	end do
!
	return
	end
!
!*******************************************
!
        subroutine mimari(xx,n,xmin,xmax,imin,imax,rnull)
!
! computes min/max of vector
!
! xx            vector
! n             dimension of vector
! xmin,xmax     min/max value in vector
! imin,imax     pointer to min/max value in vector
! rnull         invalid value
!
        implicit none
!
        integer n,i,nmin
        integer imin,imax
        double precision xx(n)
        double precision xmin,xmax,x,rnull

        do i=1,n
          if(xx(i).ne.rnull) goto 1
        end do
    1   continue

        if(i.le.n) then
          xmax=xx(i)
          xmin=xx(i)
          imin=i
          imax=i
        else
          xmax=rnull
          xmin=rnull
          imin=0
          imax=0
        end if

        nmin=i+1

        do i=nmin,n
          x=xx(i)
          if(x.ne.rnull) then
            if(x.gt.xmax) then
                xmax=x
                imax=i
            end if
            if(x.lt.xmin) then
                xmin=x
                imin=i
            end if
          end if
        end do
!
        return
        end

!*******************************************

	subroutine mima2i(text,n,a1,a2)

! writes min/max of 2 arrays to terminal

	implicit none

	character*(*) text
	integer n
	double precision a1(1), a2(1)

	double precision r1min,r1max,r2min,r2max

	call mima(a1,n,r1min,r1max)
	call mima(a2,n,r2min,r2max)

	write(6,'(a)') text
	write(6,*) 'min/max: ',r1min,r1max
	write(6,*) 'min/max: ',r2min,r2max

	end

!*******************************************

	subroutine addvec(a1,a2,n)
!
! adds two vectors
!
! results are passed back in a1
!
! a1		first array
! a2		second array
! n		dimension of vectors
!
        implicit none
!
        integer n,i
        double precision a1(n),a2(n)
!
	do i=1,n
          a1(i)=a1(i)+a2(i)
	end do
!
	return
	end
!
!*******************************************
!
	subroutine mulvec(a1,a2,n)
!
! multiplies two vectors
!
! results are passed back in a1
!
! a1		first array
! a2		second array
! n		dimension of vectors
!
        implicit none
!
        integer n,i
        double precision a1(n),a2(n)
!
	do i=1,n
          a1(i)=a1(i)*a2(i)
	end do
!
	return
	end
!
!*******************************************
!
	subroutine zervec(a1,n)
!
! initializes vector
!
! a1		vector
! n		dimension of vector
!
        implicit none
!
        integer n,i
        double precision a1(n)
!
	do i=1,n
          a1(i)=0.
	end do
!
	return
	end
!----------------------------------------------------------------------
        end module vec_util
!----------------------------------------------------------------------
