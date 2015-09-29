
!==================================================================
	module mod_nudging
!==================================================================

	implicit none

	integer, private, save :: nkn_nudging = 0
	integer, private, save :: nel_nudging = 0
	integer, private, save :: nlv_nudging = 0
        
	real, save :: anpar = 0.		!implicit parameter for nudging
	real, save :: taudefvel = 0.		!default tau for velocities

        real, allocatable, save :: andgzv(:)	!contribution to zeta
        real, allocatable, save :: tauvel(:,:)	!weighting for vel
        real, allocatable, save :: uobs(:,:)	!observations for x vel
        real, allocatable, save :: vobs(:,:)	!observations for y vel

!==================================================================
	contains
!==================================================================

	subroutine mod_nudging_init(nkn,nel,nlv)

	integer nkn,nel,nlv

        if( nkn == nkn_nudging .and. nel == nel_nudging .and.
     +      nlv == nlv_nudging ) return

        if( nel > 0 .or. nkn > 0 .or. nlv > 0 ) then
          if( nel == 0 .or. nkn == 0 .or. nlv == 0 ) then
            write(6,*) 'nel,nkn,nlv: ',nel,nkn,nlv
            stop 'error stop mod_nudging_init: incompatible parameters'
          end if
        end if

        if( nkn == nkn_nudging ) return

	if( nkn_nudging > 0 ) then
          deallocate(andgzv)
          deallocate(tauvel)
          deallocate(uobs)
          deallocate(vobs)
        end if

        nkn_nudging = nkn
        nel_nudging = nel
        nlv_nudging = nlv
        
        if( nkn == 0 ) return
        
        allocate(andgzv(nkn))
        allocate(tauvel(nlv,nel))
        allocate(uobs(nlv,nel))
        allocate(vobs(nlv,nel))

	andgzv = 0.
	tauvel = 0.
	uobs = 0.
	vobs = 0.
        
        end subroutine mod_nudging_init 

!==================================================================
        end module mod_nudging
!==================================================================

