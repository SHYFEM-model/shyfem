        module mod_geom_dynamic

        implicit none

	integer, private, save :: nkn_geom_dynamic = 0
        integer, private, save :: nel_geom_dynamic = 0
	
	integer, allocatable, save :: iwegv(:)
	integer, allocatable, save :: iwetv(:)
	integer, allocatable, save :: inodv(:)

	contains

!******************************************************

	subroutine mod_geom_dynamic_init(nkn,nel)

	integer nkn
	integer nel

	if( nkn == nkn_geom_dynamic .and. 
     +  nel == nel_geom_dynamic ) return 

	if( nkn > 0 .or. nel > 0 ) then
          if( nkn == 0 .or. nel == 0 ) then
            write(6,*) 'nkn,nel: ',nkn,nel
            stop 'error stop mod_geom_dynamic_init: 
     +            incompatible parameters'
          end if
        end if

	if( nel_geom_dynamic > 0 ) then
	  deallocate(iwegv)
	  deallocate(iwetv)
          deallocate(inodv)
        end if
	
	nel_geom_dynamic = nel
	nkn_geom_dynamic = nkn	 
	 
	if( nkn == 0 ) return
    	
	allocate(iwegv(nel))
	allocate(iwetv(nel))
	allocate(inodv(nkn))
	
       	end subroutine mod_geom_dynamic_init	

!*****************************************************

	end module mod_geom_dynamic


	

	


