
	real grand

	t = 0
	dt = 1./24.

    1	continue
	  read(5,'(a4,i4)',end=2) word,ival
	  i = i + 1
	  t = t + dt
	  val = ival
	  r = dt * ( grand(iseed) - 0.5 )
	  taux = t + r
	  write(6,*) taux,val
	goto 1
    2	continue
	end
