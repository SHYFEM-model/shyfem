c
c $Id: newdum.f,v 1.5 2010-02-22 15:38:36 georg Exp $
c
c dummy routines for compatibility with 2D version
c
c contents :
c
c subroutine inclos             initializes closing sections
c subroutine rdclos(isc)        reads closing sections
c subroutine ckclos             post-processes closing sections
c subroutine prclos             prints info on closing sections
c subroutine tsclos             tests closing sections
c
c subroutine inoxy              initializes oxygen module
c subroutine rdoxy              reads oxygen module
c subroutine ckoxy              post-processes oxygen module
c subroutine proxy              prints info on oxygen module
c subroutine tsoxy              tests oxygen module
c
c revision log :
c
c 22.05.1998	ggu	created for closing sections
c 22.01.1999	ggu	oxygen module added
c 19.02.2010	ggu	massconc eliminated
c
c*****************************************************************

c dummies for closing section

c        subroutine inclos
c	end
c
c        subroutine rdclos(isc)
c	integer isc
c	write(6,*) 'Closing sections not yet supported in 3d version'
c	write(6,*) 'closing section : ',isc
c	stop 'error stop rdclos'
c	end
c
c        subroutine ckclos
c	end
c
c        subroutine prclos
c	end
c
c        subroutine tsclos
c	end

c*****************************************************************

c dummies for oxygen module

        subroutine inoxy
	implicit none
	end

        subroutine rdoxy
	implicit none
	write(6,*) 'Oxygen module not yet supported in 3d version'
	stop 'error stop rdoxy'
	end

        subroutine ckoxy
	implicit none
	end

        subroutine proxy
	implicit none
	end

        subroutine tsoxy
	implicit none
	end

c*****************************************************************

	subroutine basin_read(nin)

	implicit none

	integer nin

	include 'param.h'

	call sp13rr(nin,nkndim,neldim)

	end

c*****************************************************************

