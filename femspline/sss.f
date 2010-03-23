
	integer ndim
	parameter(ndim=10000)

	logical bperiod
	real t(ndim)
	real v(ndim)

	bperiod = .false.
	sigma = 1.
	sigma = 0.25
	i = 0

    1	continue
	  read(5,*,end=2) taux,vaux
	  i = i + 1
	  if( i .gt. ndim ) stop 'dimension...'
	  t(i) = taux
	  v(i) = vaux
	goto 1
    2	continue
	n = i

	call gsmooth(n,t,v,sigma,bperiod)

	do i=1,n
	  write(6,*) t(i),v(i)
	end do

	end

c********************************************************

