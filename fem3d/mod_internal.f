
	module mod_internal

	implicit none

!	real rdistv(nkndim)
!       common /rdistv/rdistv

!       real fcorv(neldim)
!       common /fcorv/fcorv
!       real fxv(nlvdim,neldim)          !new HYDRO debora
!       common /fxv/fxv
!       real fyv(nlvdim,neldim)
!       common /fyv/fyv

!       save /rdistv/,/fcorv/,/fxv/,/fyv/

!       integer iuvfix(neldim)
!       common /iuvfix/iuvfix
!       save /iuvfix/

!       double precision ddxv(2*nlvdim,neldim)  !ASYM
!       double precision ddyv(2*nlvdim,neldim)  !ASYM
!       common /ddxv/ddxv, /ddyv/ddyv
!       save /ddxv/, /ddyv/


        integer, private, save :: nkn_internal = 0
        integer, private, save :: nel_internal = 0
        integer, private, save :: nlv_internal = 0
        
        real, allocatable, save :: rdistv(:)
        real, allocatable, save :: fcorv(:)
        real, allocatable, save :: fxv(:,:)
        real, allocatable, save :: fyv(:,:)
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
     +  .and. nel == nel_internal) return

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
         allocate (iuvfix(nel))
         allocate (ddxv(2*nlv,nel))
         allocate (ddyv(2*nlv,nel))
        
        end subroutine mod_internal_init 

!*****************************************************

        end module mod_internal

