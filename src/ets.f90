!
! $Id: subetsa.f,v 1.6 2001/11/16 07:35:43 georg Exp $
!
! ETS file administration routines (deals with section $EXTTS)
!
! revision log :
!
! 24.01.2014    ggu     copied from subexta.f
!
!******************************************************************
!******************************************************************
!******************************************************************

!==================================================================
        module ets
!==================================================================

        implicit none

        integer, save :: nets = 0
        integer, save, allocatable :: nkets(:)		!node numbers
        character*80, save, allocatable :: chets(:)	!description

        double precision, save, allocatable :: xets(:)		!x-coordinates
        double precision, save, allocatable :: yets(:)		!y-coordinates

        integer, save, allocatable :: ilets(:)		!layers
        double precision, save, allocatable :: hets(:)		!depth

        double precision, save, allocatable :: outets(:,:)		!aux array

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

