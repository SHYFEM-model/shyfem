!
! implements flood fill for finite element grid
!
!*********************************************************************

	subroutine flood_fill(color)

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

	integer ngood,nbad
	integer ic
	integer it,is,isc
	integer ie,ii,k
	integer coloraux(nkn)

	coloraux = 0

	nbad = count( color < 0 )

	do

	  ngood = count( color > 0 )

          write(6,*) 'influence: ',ngood,nbad,nkn,nkn-ngood-nbad

	  if( nkn-ngood-nbad == 0 ) exit

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
                coloraux(k) = ic
              end do
            else if( it == 2 ) then	!two nodes are colored
              is = 6 - is		!this is the uncolored node
              isc = mod(is,3) + 1	!take color from this node
              k = nen3v(isc,ie)
              ic = color(k)
              k = nen3v(is,ie)
              coloraux(k) = ic
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

