c
c $Id: newdum.f,v 1.5 2010-02-22 15:38:36 georg Exp $
c
c*****************************************************************


	subroutine levels_init(nkn,nel,nl)
	implicit none
	integer nkn,nel,nl
	include 'param.h'
	include 'nlevel.h'
	if( nl > nlvdim ) stop 'error stop levels_init: nlvdim'
	nlvdi = nlvdim
	end

	subroutine levels_reinit(nl)
	implicit none
	integer nl
	include 'param.h'
	include 'nlevel.h'
	if( nl > nlvdim ) stop 'error stop levels_reinit: nlvdim'
	nlvdi = nlvdim
	end

	subroutine levels_hlv_init(nl)
	implicit none
	integer nl
	include 'param.h'
	include 'nlevel.h'
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

