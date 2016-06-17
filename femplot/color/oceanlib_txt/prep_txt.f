
	program prep

	implicit none

	real col(3)

	do
	  read(5,*,end=1) col
	  write(6,1000) '[',col(1),',',col(2),',',col(3),'],'
	end do
    1	continue

	stop
 1000	format(a,3(f9.6,a))
	end

