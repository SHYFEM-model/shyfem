c
c $Id: subetsa.f,v 1.6 2001/11/16 07:35:43 georg Exp $
c
c ETS file administration routines (deals with section $EXTTS)
c
c revision log :
c
c 24.01.2014    ggu     copied from subexta.f
c
c******************************************************************
c******************************************************************
c******************************************************************

!==================================================================
        module ets
!==================================================================

        implicit none

        integer, save :: nets = 0
        integer, save, allocatable :: nkets(:)		!node numbers
        character*80, save, allocatable :: chets(:)	!description

        real, save, allocatable :: xets(:)		!x-coordinates
        real, save, allocatable :: yets(:)		!y-coordinates

        integer, save, allocatable :: ilets(:)		!layers
        real, save, allocatable :: hets(:)		!depth

        real, save, allocatable :: outets(:,:)		!aux array

!==================================================================
        contains
!==================================================================

	subroutine ets_init_module(n)

	integer n

	nets = n

        allocate(nkets(n))
        allocate(chets(n))

        allocate(xets(n))
        allocate(yets(n))
        allocate(ilets(n))
        allocate(hets(n))

	end subroutine ets_init_module

!==================================================================
        end module ets
!==================================================================

