!----------------------------------------------------------------------
        module befor_after
!----------------------------------------------------------------------
        contains
!----------------------------------------------------------------------
!********************************************************************

	subroutine do_befor

! to do in time loop before time step

        use tidef
        use chezy
        use modls

	implicit none

	include 'modules.h'

	include 'femtime.h'

	call modules(M_BEFOR)

        call tideforce(it)       !tidal potential !ccf

	call adjust_chezy

	end

!********************************************************************

	subroutine do_after

! to do in time loop after time step

        use restart
        use flux
        use volume
        use modls
        use ous_admin
        use flux_box
        use residual
        use custom_admin

        implicit none

        include 'modules.h'

        include 'femtime.h'

        call modules(M_AFTER)

!	call wrouta
        call wrousa
!	call wrexta(it)
        call wrflxa(it)
        call wrvola(it)
        call wrboxa(it)

        call resid
        call rmsvel

        call admrst             !restart

!        call tsmed
        call ts_shell

!	call wrnetcdf		!output in netcdf format - not supported

        call custom(it)

        end

!*******************************************************************

!----------------------------------------------------------------------
        end module befor_after
!----------------------------------------------------------------------
