c
c $Id: ecological_dummy.f,v 1.1 2008-04-29 15:31:32 georg Exp $
c
c ecological dummy module
c
c revision log :
c
c 29.04.2008    ggu     bfm model integrated in main branch
c 18.02.2011    ggu     general framework for ecological model
c
c**************************************************************

        subroutine ecological_module(it,dt)

c general interface to ecological module

        implicit none

        integer it
        real dt

	integer ibfm,ibio
	real getpar

	integer icall
	save icall
	data icall / 0 /

	if( icall .eq. -1 ) return

        ibfm = nint(getpar('ibfm'))
	ibio = nint(getpar('ibio'))

	ibfm = 0
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

c**************************************************************

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

c**************************************************************

