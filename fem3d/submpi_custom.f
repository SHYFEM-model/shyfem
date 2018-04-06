
!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine make_custom_domain_area(area_node)

	use basin
	use shympi

	implicit none

	integer area_node(nkn)

	integer ie,ii
	real r,h

	if( nkn /= 225 .or. nel /= 384 ) then
	  if( n_threads > 1 ) then
	    write(6,*) 'nkn,nel: ',nkn,nel
	    write(6,*) 'expecting: ',225,384
	    write(6,*) 'cannot make custom domain'
	    !stop 'error stop make_domain_area: wrong basin'
	    return
	  end if
	end if

	if( n_threads == 1 ) then
	  return
	else if( n_threads == 2 ) then
	  call make_domain_area_2(area_node)
	else if( n_threads == 3 ) then
	  call make_domain_area_3(area_node)
	else if( n_threads == 4 ) then
	  call make_domain_area_4(area_node)
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

	subroutine make_domain_area_2(area_node)

	use basin
	use shympi

	implicit none

	integer k
	integer area_node(nkn)

	do k=1,nkn
	  if( ygv(k) > 3100.  ) then
	    area_node(k) = 1
	  end if
	end do

	end

!*****************************************************************

	subroutine make_domain_area_3(area_node)

	use basin
	use shympi

	implicit none

	integer k
	integer area_node(nkn)

	do k=1,nkn
	  if( ygv(k) > 4100.  ) then
	    area_node(k) = 2
	  else if( ygv(k) > 2100.  ) then
	    area_node(k) = 1
	  end if
	end do

	end

!*****************************************************************

	subroutine make_domain_area_4(area_node)

	use basin
	use shympi

	implicit none

	integer k
	integer area_node(nkn)

	do k=1,nkn
	  if( ygv(k) > 3100.  ) then
	    area_node(k) = 2
	  end if
	  if( xgv(k) > 100.  ) then
	    area_node(k) = area_node(k) + 1
	  end if
	end do

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

