
c********************************************************************
c version for reals
c********************************************************************

	subroutine set_3d_r_array(n1,n2,n3,ra,rval)

	implicit none

	integer n1,n2,n3
	real ra(n1,n2,n3)
	real rval

	integer i1,i2,i3

	do i3=1,n3
	  do i2=1,n2
	    do i1=1,n1
	      ra(i1,i2,i3) = rval
	    end do
	  end do
	end do

	end

c********************************************************************

	subroutine set_2d_r_array(n1,n2,ra,rval)

	implicit none

	integer n1,n2
	real ra(n1,n2)
	real rval

	integer i1,i2

	do i2=1,n2
	  do i1=1,n1
	    ra(i1,i2) = rval
	  end do
	end do

	end

c********************************************************************

	subroutine set_1d_r_array(n,ra,rval)

	implicit none

	integer n
	real ra(n)
	real rval

	integer i

	do i=1,n
	  ra(i) = rval
	end do

	end

c********************************************************************
c version for double precision
c********************************************************************

	subroutine set_3d_d_array(n1,n2,n3,da,dval)

	implicit none

	integer n1,n2,n3
	double precision da(n1,n2,n3)
	double precision dval

	integer i1,i2,i3

	do i3=1,n3
	  do i2=1,n2
	    do i1=1,n1
	      da(i1,i2,i3) = dval
	    end do
	  end do
	end do

	end

c********************************************************************

	subroutine set_2d_d_array(n1,n2,da,dval)

	implicit none

	integer n1,n2
	double precision da(n1,n2)
	double precision dval

	integer i1,i2

	do i2=1,n2
	  do i1=1,n1
	    da(i1,i2) = dval
	  end do
	end do

	end

c********************************************************************

	subroutine set_1d_d_array(n,da,dval)

	implicit none

	integer n
	double precision da(n)
	double precision dval

	integer i

	do i=1,n
	  da(i) = dval
	end do

	end

c********************************************************************

