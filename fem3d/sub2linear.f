!
! utility routines to convert between 3D arrays and linear array
!
! contents :
!
! revision log :
!
! 09.11.2018    ggu     written from scartch
! 11.11.2018    ggu     some bug fixes
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
        real rlin(nlin)
        integer nlin

        integer i,j,lmax,nl,ne

        nl = 0

        do i=1,m
          do j=1,n
            lmax = min(nlvddi,il(j))
	    ne = nl + lmax
	    if( ne > nlin ) goto 99
            rlin(nl+1:ne) = vals(1:lmax,j,i)
	    nl = ne
          end do
        end do

	nlin = nl

	return
   99	continue
	write(6,*) nl,ne,nlin
	stop 'error stop vals2linear: nl>nlin'
        end

!************************************************************************

        subroutine linear2vals(nlvddi,n,m,il,vals,rlin,nlin)

        implicit none

        integer nlvddi,n,m
        integer il(n)
        real vals(nlvddi,n,m)
        real rlin(nlvddi*n*m)
        integer nlin

        integer i,j,lmax,nl,ne

        nl = 0

        do i=1,m
          do j=1,n
            lmax = min(nlvddi,il(j))
	    ne = nl + lmax
	    if( ne > nlin ) goto 99
            vals(1:lmax,j,i) = rlin(nl+1:ne)
	    nl = ne
          end do
        end do

	if( nl /= nlin ) stop 'error stop linear2vals: nl/=nlin'

	return
   99	continue
	write(6,*) nl,ne,nlin
	stop 'error stop linear2vals: nl>nlin'
        end

!************************************************************************

