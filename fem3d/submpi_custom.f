
!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine make_custom_domain_area

	use basin
	use shympi

	implicit none

	integer ie,ii
	real r,h

	node_area = 0

	if( nkn /= 225 .or. nel /= 384 ) then
	  if( n_threads > 1 ) then
	    write(6,*) 'nkn,nel: ',nkn,nel
	    write(6,*) 'expecting: ',225,384
	    stop 'error stop make_domain_area: wrong basin'
	  end if
	end if

	if( n_threads == 1 ) then
	  return
	else if( n_threads == 2 ) then
	  call make_domain_area_2
	else if( n_threads == 3 ) then
	  call make_domain_area_3
	else if( n_threads == 4 ) then
	  call make_domain_area_4
	else
	  write(6,*) 'n_threads = ',n_threads
	  stop 'error stop make_domain_area: cannot handle'
	end if

	do ie=1,nel
	  do ii=1,3
	    call random_number(r)
	    h = 10.* r
	    h = 2.* r
	    !hm3v(ii,ie) = h
	    !write(6,*) ie,ii,r,h
	  end do
	end do

	end

!*****************************************************************

	subroutine make_domain_area_2

	use basin
	use shympi

	implicit none

	integer k

	do k=1,nkn
	  if( ygv(k) > 3100.  ) then
	    node_area(k) = 1
	  end if
	end do

	end

!*****************************************************************

	subroutine make_domain_area_3

	use basin
	use shympi

	implicit none

	integer k

	do k=1,nkn
	  if( ygv(k) > 4100.  ) then
	    node_area(k) = 2
	  else if( ygv(k) > 2100.  ) then
	    node_area(k) = 1
	  end if
	end do

	end

!*****************************************************************

	subroutine make_domain_area_4

	use basin
	use shympi

	implicit none

	integer k

	do k=1,nkn
	  if( ygv(k) > 3100.  ) then
	    node_area(k) = 2
	  end if
	  if( xgv(k) > 100.  ) then
	    node_area(k) = node_area(k) + 1
	  end if
	end do

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

