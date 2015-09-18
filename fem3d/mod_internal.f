
	module mod_internal

	implicit none

        integer, private, save :: nkn_internal = 0
        integer, private, save :: nel_internal = 0
        integer, private, save :: nlv_internal = 0
        
        real, allocatable, save :: rdistv(:)
        real, allocatable, save :: fcorv(:)
        real, allocatable, save :: fxv(:,:)
        real, allocatable, save :: fyv(:,:)
        real, allocatable, save :: momentxv(:,:)
        real, allocatable, save :: momentyv(:,:)
        real, allocatable, save :: iuvfix(:)
        double precision, allocatable, save :: ddxv(:,:)
        double precision, allocatable, save :: ddyv(:,:)

        contains

*******************************************************************

	subroutine mod_internal_init(nkn,nel,nlv)

	integer nkn
        integer nel
        integer nlv

        if( nkn == nkn_internal .and. nel == nel_internal
     +  .and. nlv == nlv_internal) return

        if( nkn > 0 .or. nel > 0 .or. nlv > 0) then
          if( nkn == 0 .or. nel == 0 .or. nlv == 0) then
            write(6,*) 'nkn,nel,nlv: ',nkn,nel,nlv
            stop 'error stop mod_internal_init: incompatible params'
          end if
        end if

	if( nkn_internal > 0 ) then
         deallocate(rdistv)
         deallocate(fcorv)
         deallocate(fxv)
         deallocate(fyv)
         deallocate(momentxv)
         deallocate(momentyv)
         deallocate(iuvfix)
         deallocate(ddxv)
         deallocate(ddyv)
        end if

        nel_internal = nel
        nkn_internal = nkn 
        nlv_internal = nlv       
        
        if( nkn == 0 ) return
        
         allocate (rdistv(nkn))
         allocate (fcorv(nel))
         allocate (fxv(nlv,nel))
         allocate (fyv(nlv,nel))
         allocate (momentxv(nlv,nkn))
         allocate (momentyv(nlv,nkn))
         allocate (iuvfix(nel))
         allocate (ddxv(2*nlv,nel))
         allocate (ddyv(2*nlv,nel))
        
        end subroutine mod_internal_init 

!*****************************************************

        end module mod_internal

