
	program array_assign

! tests array assignments

	implicit none

	interface
	subroutine write_array(a1,a2)
        integer a1(:),a2(:)
	end
	end interface


	integer i
	integer, allocatable :: a(:)
	integer, allocatable :: b(:)
	integer, allocatable :: c(:)

	allocate(a(5),b(7),c(9))

	do i=1,7
	  b(i) = i
	end do

	c = 0
	call write_array(a,b)
	call write_array(c,b)

	end

	subroutine write_array(a1,a2)

	implicit none

	integer a1(:),a2(:)

	a1 = a2
	write(6,*) a1

	end

