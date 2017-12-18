!
! implements flood fill for finite element grid
!
!*********************************************************************

! fills color starting from some colored points
!
! ic = color(k), k=1,nkn
!
! ic == 0 indicates a point that still has to be colored
! ic > 0  indicates a node colored with color ic
! ic < 0  indicates a point that not has to be touched
!
! at the end of the routine all points are different from 0

!*********************************************************************

	subroutine flood_fill(color)

! classic flood fill

	use basin

	implicit none

	integer color(nkn)

	integer iter

	call flood_fill_internal(color,.false.,iter)

	end

!*********************************************************************

	subroutine flood_fill_progessive(color)

! increase progressively color by one

	use basin

	implicit none

	integer color(nkn)

	integer iter

	call flood_fill_internal(color,.true.,iter)

	end

!*********************************************************************

	subroutine flood_fill_iter(color,iter)

! increase progressively color by one

	use basin

	implicit none

	integer color(nkn)

	integer iter

	call flood_fill_internal(color,.false.,iter)

	end

!*********************************************************************

	subroutine flood_fill_internal(color,bprogress,iter)

! fills color starting from some colored points
!
! ic = color(k), k=1,nkn
!
! ic == 0 indicates a point that still has to be colored
! ic > 0  indicates a node colored with color ic
! ic < 0  indicates a point that not has to be touched
!
! at the end of the routine all points are different from 0

	use basin

	implicit none

	integer color(nkn)
	logical bprogress
	integer iter

	integer ngood,nbad
	integer ic,inc
	integer it,is,isc
	integer ie,ii,k
	integer coloraux(nkn)

	coloraux = 0

	iter = 0
	inc = 0
	if( bprogress ) inc = 1

	nbad = count( color < 0 )

	do

	  ngood = count( color > 0 )

          !write(6,*) 'influence: ',iter,ngood,nbad,nkn,nkn-ngood-nbad

	  if( nkn-ngood-nbad == 0 ) exit

	  iter = iter + 1

          do ie=1,nel

            it = 0
            is = 0
            do ii=1,3
              k = nen3v(ii,ie)
              if( color(k) > 0 ) then
                it = it + 1
                is = is + ii
              end if
            end do

            if( it == 1 ) then        !one node is colored
              k = nen3v(is,ie)
              ic = color(k)
              do ii=1,3
                k = nen3v(ii,ie)
                coloraux(k) = ic + inc
              end do
            else if( it == 2 ) then	!two nodes are colored
              is = 6 - is		!this is the uncolored node
              isc = mod(is,3) + 1	!take color from this node
              k = nen3v(isc,ie)
              ic = color(k)
              k = nen3v(is,ie)
              coloraux(k) = ic + inc
            end if

            if( it > 0 ) then        !check
              do ii=1,3
                k = nen3v(ii,ie)
                if( color(k) == 0 .and. coloraux(k) == 0 ) then
                  write(6,*) 'internal error...... '
                  write(6,*) ie,it,is
                  stop 'error stop: internal error'
                end if
              end do
            end if

          end do

	  where( coloraux > 0 .and. color == 0 ) color = coloraux

        end do

	end

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine fill_areas(color,ncol)

! fills progressively one area (ic) starting from 1

	use basin
	use mod_geom

	implicit none

	integer color(nkn)
	integer ncol

	integer ic
	integer color_aux(nkn)
	integer color_out(nkn)
	integer kant(2,nkn)

	color_out = 0
	if( .not. allocated(kantv) ) goto 99
	kant = kantv
	call adjust_kant(color,kant)

	do ic=1,ncol
	  color_aux = color
	  call fill_area(color_aux,kant,ic)
	  where( color_aux > 0 ) color_out = color_aux
	end do

	color = color_out

	return
   99	continue
	stop 'error stop fill_border: kantv not allocated'
	end

!*********************************************************************

	subroutine fill_area(color,kant,ic)

! fills progressively one area (ic) starting from 1

	use basin

	implicit none

	integer color(nkn)
	integer kant(2,nkn)
	integer ic

	integer color_aux(nkn)

	color_aux = -1
	where( color == ic ) color_aux = 0
	call fill_border(color_aux,kant,1)
	call flood_fill_progessive(color_aux)

	color = color_aux

	end

!*********************************************************************

	subroutine fill_border(color,kant,ic)

! fills border of 0 area with ic

	use basin

	implicit none

	integer color(nkn)
	integer kant(2,nkn)
	integer ic

	integer ie,ii,k,kb
	logical bhaszero,bhasneg
	integer kanta(2,nkn)

!----------------------------------------------
! fill border to negative values with ic
!----------------------------------------------

	do ie=1,nel
	  bhaszero = .false.
	  bhasneg = .false.
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( color(k) == 0 ) bhaszero = .true.
	    if( color(k) < 0 ) bhasneg = .true.
	  end do
	  if( bhaszero .and. bhasneg ) then
	    do ii=1,3
	      k = nen3v(ii,ie)
	      if( color(k) == 0 ) color(k) = ic
	    end do
	  end if
	end do
	
!----------------------------------------------
! fill external border with ic - only those that are not compl. in area
!----------------------------------------------

	do k=1,nkn
	  kb = kant(1,k)
	  if( kb > 0 ) then
	    if( color(k) == 0 ) color(k) = ic
	  end if
	end do

	end

!*********************************************************************

	subroutine adjust_kant(color,kant)

	use basin

	implicit none

	integer color(nkn)
	integer kant(2,nkn)

	logical bunique
	integer k,kk,i,n,col
	integer list(nkn)

	do k=1,nkn
	  if( kant(1,k) > 0 ) then
	    call get_kant_node_list(kant,k,list,n)
	    bunique = .true.
	    col = color(list(1))
	    do i=2,n
	      if( color(list(i)) /= col ) bunique = .false.
	    end do
	    if( bunique ) then
	      do i=1,n
		kk = list(i)
	        kant(:,kk) = 0
	      end do
	    else
	      do i=1,n
		kk = list(i)
	        kant(:,kk) = -kant(:,kk)
	      end do
	    end if
	  end if
	end do

	where( kant < 0 ) kant = -kant

	end

!*********************************************************************

	subroutine get_kant_node_list(kant,k,list,n)

	use basin

	implicit none

	integer k,n
	integer kant(2,nkn)
	integer list(nkn)

	integer kk

	n = 0
	kk = k

	do
	  n = n + 1
	  list(n) = kk
	  kk = kant(1,kk)
	  if( kk .eq. k ) exit
	end do

	end

!*********************************************************************

