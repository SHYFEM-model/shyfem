
        module bclfix_util

        implicit none

        integer, private, save :: nel_bclfix = 0
        integer, private, save :: nkn_bclfix = 0
        integer, private, save :: nlv_bclfix = 0

        double precision, allocatable, save :: ielfix(:,:)	! total number of nodes for ele
        double precision, allocatable, save :: tnudgev(:)	! nudging coefficient array
        double precision, allocatable, save :: ubound(:,:)	! current velocity in x for boundary [m/s]
        double precision, allocatable, save :: vbound(:,:)	! current velocity in y for boundary [m/s]

        contains

!************************************************************

        subroutine mod_bclfix_util_init(nkn,nel,nlv)

        integer  :: nkn
        integer  :: nel
        integer  :: nlv

        if( nlv == nlv_bclfix .and. nel == nel_bclfix .and. nkn == nkn_bclfix ) return

        if( nlv > 0 .or. nel > 0 .or. nkn > 0 ) then
          if( nlv == 0 .or. nel == 0 .or. nkn == 0 ) then
            write(6,*) 'nlv,nel,nkn: ',nlv,nel,nkn
            stop 'error stop mod_bclfix_util_init: incompatible parameters'
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

        end subroutine mod_bclfix_util_init

!************************************************************

        end module bclfix_util

