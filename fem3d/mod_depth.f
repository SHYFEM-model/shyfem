
	module mod_depth

	implicit none

! 	 real hkv(nkndim)
!        common /hkv/hkv
!        real hev(neldim)
!        common /hev/hev

!        real hdknv(nlvdim,nkndim)
!        common /hdknv/hdknv
!        real hdkov(nlvdim,nkndim)
!        common /hdkov/hdkov

!        real hdenv(nlvdim,neldim)
!        common /hdenv/hdenv
!        real hdeov(nlvdim,neldim)
!        common /hdeov/hdeov

!        save /hkv/,/hev/,/hdknv/,/hdkov/,/hdenv/,/hdeov/

!        real hkv_min(nkndim), hkv_max(nkndim)
!        common /hkv_min/hkv_min, /hkv_max/hkv_max
!        save /hkv_min/,/hkv_max/

        integer, private, save :: nkn_depth = 0
        integer, private, save :: nel_depth = 0
        integer, private, save :: nlv_depth = 0
        
        real, allocatable, save :: hkv(:)
        real, allocatable, save :: hev(:)
        real, allocatable, save :: hdknv(:,:)
        real, allocatable, save :: hdkov(:,:)
        real, allocatable, save :: hdenv(:,:)
        real, allocatable, save :: hdeov(:,:)

        real, allocatable, save :: hkv_min(:)
        real, allocatable, save :: hkv_max(:)

        contains

*******************************************************************

	subroutine mod_depth_init(nkn,nel,nlv)

	integer nkn
        integer nel
        integer nlv

        if( nkn == nkn_depth .and. nel == nel_depth
     +  .and. nlv == nlv_depth ) return

        if( nkn > 0 .or. nel > 0 .or. nlv > 0 ) then
          if( nkn == 0 .or. nel == 0 .or. nlv == 0 ) then
            write(6,*) 'nkn,nel,nlv: ',nkn,nel,nlv
            stop 'error stop mod_depth_init: incompatible params'
          end if
        end if

	if( nkn_depth > 0 ) then
         deallocate(hkv)
         deallocate(hdknv)
         deallocate(hdkov)
         deallocate(hkv_min)
         deallocate(hkv_max)
        end if

	if( nel_depth > 0 ) then
         deallocate(hev)
         deallocate(hdenv)
         deallocate(hdeov)
        end if

         nkn_depth = nkn 
         nel_depth = nel
         nlv_depth = nlv       
        
        if( nkn == 0 ) return
        
         allocate(hkv(nkn))
         allocate(hev(nel))
         allocate(hdknv(nlv,nkn))
         allocate(hdkov(nlv,nkn))
         allocate(hdenv(nlv,nel))
         allocate(hdeov(nlv,nel))
         allocate(hkv_min(nkn))
         allocate(hkv_max(nkn))

        end subroutine mod_depth_init 

!*****************************************************

        end module mod_depth

