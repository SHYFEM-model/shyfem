!
! $Id: ecological_dummy.f,v 1.1 2008-04-29 15:31:32 georg Exp $
!
! ecological dummy module
!
! revision log :
!
! 29.04.2008    ggu     bfm model integrated in main branch
! 18.02.2011    ggu     general framework for ecological model
!
!**************************************************************
!--------------------------------------------------------------
        module ecological
!--------------------------------------------------------------
        contains
!--------------------------------------------------------------

        subroutine ecological_module(it,dt)

! general interface to ecological module

        use para

        implicit none

        integer it
        double precision dt

	integer ibfm,ibio

	integer icall
	save icall
	data icall / 0 /

	if( icall .eq. -1 ) return

        ibfm = nint(getpar('ibfm'))
	ibio = nint(getpar('ibio'))

	if( ibfm .gt. 0 ) then
	  write(6,*) 'BFM module has not been linked'
	  write(6,*) 'ibfm = ',ibfm
	  write(6,*) 'You must enable this feature in Rules.make'
	  stop 'error stop ecological_module: ibfm'
	end if

	if( ibio .gt. 0 ) then
	  write(6,*) 'No ecological module has been linked'
	  write(6,*) 'ibio = ',ibio
	  write(6,*) 'You must enable this feature in Rules.make'
	  stop 'error stop ecological_module: ibio'
	end if

	icall = -1

        end

!**************************************************************

        subroutine write_restart_eco(iunit)
        implicit none
	integer iunit
	integer nstate,nkn,i
	nstate = 0
	nkn = 0
        write(iunit) nstate,nkn
	end
        subroutine skip_restart_eco(iunit)
        implicit none
	integer iunit
	integer nstate,nkn,i
        read(iunit) nstate,nkn
        do i=1,nstate
          read(iunit)
        end do
	end
        subroutine read_restart_eco(iunit)
        implicit none
	integer iunit
	call skip_restart_eco(iunit)
	end

!**************************************************************

!--------------------------------------------------------------
        end module ecological
!--------------------------------------------------------------
