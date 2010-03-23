
c***********************************************************************

	subroutine read_ts(file,ndim,nvar,n,itime,array)

c reads time series file

	implicit none

	character*(*) file		!file name
	integer ndim			!dimension of arrays
	integer nvar			!number of expected variables in file
	integer n			!number of records read (return)
	integer itime(ndim)		!time column
	real array(ndim,nvar)		!value column(s)

	integer naux
	parameter(naux=100)

	integer iunit,i,it
	real aux(naux)
	integer ifileo

	if( nvar .gt. naux ) stop 'error stop read_ts: naux'

	iunit = 55
        iunit = ifileo(iunit,file,'form','old')
        if( iunit .le. 0 ) stop 'error stop read_ts: no such file'

	n = 0
    1	continue
	  read(iunit,*,end=2) it,(aux(i),i=1,nvar)
	  n = n + 1
	  if( n .gt. ndim ) stop 'error stop read_ts: ndim'
	  itime(n) = it
	  do i=1,nvar
	    array(n,i) = aux(i)
	  end do
	goto 1
    2	continue

	close(iunit)

	end

c***********************************************************************

