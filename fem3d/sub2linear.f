!
! utility routines to convert between 3D arrays and linear array
!
! contents :
!
! revision log :
!
! 09.11.2018    ggu     written from scartch
!
! usage :
!
!	call count_linear
!	allocate(rlin)
!	read rlin
!	call linear2vals
!
!	call count_linear
!	allocate(rlin)
!	call vals2linear
!	write rlin
!
!************************************************************************

        subroutine count_linear(nlvddi,n,m,il,nlin)

        implicit none

        integer nlvddi,n,m
        integer il(n)
	integer nlin

	integer j,lmax

        nlin = 0

        do j=1,n
          lmax = min(nlvddi,il(j))
	  nlin = nlin + lmax
	end do

	nlin = nlin * m

	end

!************************************************************************

        subroutine vals2linear(nlvddi,n,m,il,vals,rlin,nlin)

        implicit none

        integer nlvddi,n,m
        integer il(n)
        real vals(nlvddi,n,m)
        real rlin(nlvddi*n*m)
        integer nlin

        integer i,j,lmax

        nlin = 0

        do i=1,m
          do j=1,n
            lmax = min(nlvddi,il(j))
            rlin(nlin+1:nlin+lmax) = vals(1:lmax,j,i)
	    nlin = nlin + lmax
          end do
        end do

        end

!************************************************************************

        subroutine linear2vals(nlvddi,n,m,il,vals,rlin,nlin)

        implicit none

        integer nlvddi,n,m
        integer il(n)
        real vals(nlvddi,n,m)
        real rlin(nlvddi*n*m)
        integer nlin

        integer i,j,lmax

        nlin = 0

        do i=1,m
          do j=1,n
            lmax = min(nlvddi,il(j))
            vals(1:lmax,j,i) = rlin(nlin+1:nlin+lmax)
	    nlin = nlin + lmax
          end do
        end do

        end

!************************************************************************

