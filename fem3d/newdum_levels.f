c
c $Id: newdum.f,v 1.5 2010-02-22 15:38:36 georg Exp $
c
c*****************************************************************


	subroutine levels_init(nkn,nel,nl)
	use levels, only : nlvdi,nlv

	implicit none
	integer nkn,nel,nl
	include 'param.h'
	if( nl > nlvdim ) stop 'error stop levels_init: nlvdim'
	nlvdi = nlvdim
	end

	subroutine levels_reinit(nl)
	use levels, only : nlvdi,nlv

	implicit none
	integer nl
	include 'param.h'
	if( nl > nlvdim ) stop 'error stop levels_reinit: nlvdim'
	nlvdi = nlvdim
	end

	subroutine levels_hlv_init(nl)
	use levels, only : nlvdi,nlv

	implicit none
	integer nl
	include 'param.h'
	if( nl > nlvdim ) stop 'error stop levels_reinit: nlvdim'
	nlvdi = nlvdim
	end

c*****************************************************************

        subroutine levels_get_dimension(nl)

! returns vertical dimension (static)

        integer nl

        include 'param.h'

        nl = nlvdim

        end subroutine levels_get_dimension

c*****************************************************************

