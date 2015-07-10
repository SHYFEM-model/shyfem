
        module mod_bclfix

        implicit none

        !integer ielfix(0:3,neldim)
        !common /ielfix/ielfix
        !save /ielfix/

        !integer tnudgev(neldim)
        !common /tnudgev/tnudgev
        !save /tnudgev/

        !real ubound(nlvdim,nkndim)
        !real vbound(nlvdim,nkndim)
        !common /ubound/ubound
        !common /vbound/vbound
	!save /ubound/,/vbound/

        integer, private, save :: nel_bclfix = 0
        integer, private, save :: nkn_bclfix = 0
        integer, private, save :: nlv_bclfix = 0

        real, allocatable, save :: ielfix(:,:)	! total number of nodes for ele
        real, allocatable, save :: tnudgev(:)	! nudging coefficient array
        real, allocatable, save :: ubound(:,:)	! current velocity in x for boundary [m/s]
        real, allocatable, save :: vbound(:,:)	! current velocity in y for boundary [m/s]

        contains

!************************************************************

        subroutine mod_bclfix_init(nkn,nel,nlv)

        integer  :: nkn
        integer  :: nel
        integer  :: nlv

        if( nlv == nlv_bclfix .and. nel == nel_bclfix .and. 
     +      nkn == nkn_bclfix ) return

        if( nlv > 0 .or. nel > 0 .or. nkn > 0 ) then
          if( nlv == 0 .or. nel == 0 .or. nkn == 0 ) then
            write(6,*) 'nlv,nel,nkn: ',nlv,nel,nkn
            stop 'error stop mod_bclfix_init: incompatible parameters'
          end if
        end if

        if( nkn_bclfix > 0 ) then 
          deallocate(ielfix)
          deallocate(tnudgev)
          deallocate(ubound)
          deallocate(vbound)
        end if

        nlv_bclfix = nlv
        nel_bclfix = nel
        nkn_bclfix = nkn

        if( nkn == 0 ) return

        allocate(ielfix(0:3,nel))
        allocate(tnudgev(nel))
        allocate(ubound(nlv,nkn))
        allocate(vbound(nlv,nkn))

        end subroutine mod_bclfix_init

!************************************************************

        end module mod_bclfix

