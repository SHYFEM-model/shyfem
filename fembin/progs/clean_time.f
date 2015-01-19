
	program clean_time

c cleans time series from syncronization time records

	implicit none

	integer ndim	!number of values to read
	parameter(ndim=6)

	integer idtsyn,time
	integer ival(2)
	real rval(3)

	read(5,*) idtsyn

    1	continue
	  read(5,*,end=2) time,rval,ival
	  if( mod(time,idtsyn) .eq. 0 ) goto 1
	  write(6,*) time,rval,ival
	  goto 1
    2	continue

	end

