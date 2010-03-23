c
c $Id: subssv.f,v 1.4 1999/06/22 14:47:52 georg Exp $
c
c utility routines for operations on vectors
c
c contents :
c
c subroutine mima(xx,n,xmin,xmax)	min/max of vector
c
c subroutine mimari(xx,n,xmin,xmax,imin,imax,rnull)
c                                       min/max of vector with null value
c
c subroutine mima2i(text,n,a1,a2)	writes min/max of 2 arrays to terminal
c
c subroutine addvec(a1,a2,n)		adds two vectors
c subroutine mulvec(a1,a2,n)		multiplies two vectors
c subroutine zervec(a1,n)		initializes vector
c
c revision log :
c
c 26.08.1998    ggu	routines mimari transferred from newbcl0
c 31.05.1999    ggu	new comodity routine mima2i
c
c*******************************************
c
	subroutine mima(xx,n,xmin,xmax)
c
c computes min/max of vector
c
c xx		vector
c n		dimension of vector
c xmin,xmax	min/max value in vector
c
        implicit none
c
        integer n,i
        real xx(n)
        real xmin,xmax,x
c
	xmax=xx(1)
	xmin=xmax
c
	do i=1,n
          x=xx(i)
          if(x.gt.xmax) xmax=x
          if(x.lt.xmin) xmin=x
	end do
c
	return
	end
c
c*******************************************
c
        subroutine mimari(xx,n,xmin,xmax,imin,imax,rnull)
c
c computes min/max of vector
c
c xx            vector
c n             dimension of vector
c xmin,xmax     min/max value in vector
c imin,imax     pointer to min/max value in vector
c rnull         invalid value
c
        implicit none
c
        integer n,i,nmin
        integer imin,imax
        real xx(n)
        real xmin,xmax,x,rnull

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
c
        return
        end

c*******************************************

	subroutine mima2i(text,n,a1,a2)

c writes min/max of 2 arrays to terminal

	implicit none

	character*(*) text
	integer n
	real a1(1), a2(1)

	real r1min,r1max,r2min,r2max

	call mima(a1,n,r1min,r1max)
	call mima(a2,n,r2min,r2max)

	write(6,'(a)') text
	write(6,*) 'min/max: ',r1min,r1max
	write(6,*) 'min/max: ',r2min,r2max

	end

c*******************************************

	subroutine addvec(a1,a2,n)
c
c adds two vectors
c
c results are passed back in a1
c
c a1		first array
c a2		second array
c n		dimension of vectors
c
        implicit none
c
        integer n,i
        real a1(n),a2(n)
c
	do i=1,n
          a1(i)=a1(i)+a2(i)
	end do
c
	return
	end
c
c*******************************************
c
	subroutine mulvec(a1,a2,n)
c
c multiplies two vectors
c
c results are passed back in a1
c
c a1		first array
c a2		second array
c n		dimension of vectors
c
        implicit none
c
        integer n,i
        real a1(n),a2(n)
c
	do i=1,n
          a1(i)=a1(i)*a2(i)
	end do
c
	return
	end
c
c*******************************************
c
	subroutine zervec(a1,n)
c
c initializes vector
c
c a1		vector
c n		dimension of vector
c
        implicit none
c
        integer n,i
        real a1(n)
c
	do i=1,n
          a1(i)=0.
	end do
c
	return
	end
