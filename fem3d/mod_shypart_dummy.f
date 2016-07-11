!
! mpi routines
!
! contents :
!
! revision log :
!
! 16.02.2016    CMCC     project started
!
!******************************************************************

!==================================================================
        module shypart
!==================================================================

	implicit none

	public

        integer, allocatable, save, dimension(:) :: allPartAssign

        integer, allocatable, save :: temp_ilinkv(:)
        integer, allocatable, save :: temp_lenkv(:)
        integer, allocatable, save :: temp_lenkiiv(:)
        integer, allocatable, save :: temp_linkv(:)

        integer, allocatable, save :: temp_ieltv(:,:)
        integer, allocatable, save :: temp_kantv(:,:)
 
        integer, allocatable, save :: temp_ilhv(:)


        integer, allocatable, save, dimension(:) :: tempIpv,tempIpev
        real, allocatable, save, dimension(:) :: tempXgv,tempYgv
        integer, allocatable, save, dimension(:) :: tempIarv,tempIarnv
        integer, allocatable, save, dimension(:,:) :: tempNen3v
        real, allocatable, save, dimension(:,:) :: tempHm3v

        integer,save :: n_threads = 1
        integer,save :: my_id = 0
        logical,save :: repart = .false.
        character*(3) :: phg_method='AGG'

!==================================================================
        contains
!==================================================================

	subroutine shypart_init(b_use_mpi)

	logical b_use_mpi

        return

	end subroutine shypart_init

!*****************************************************************

	subroutine shypart_setup

	implicit none

        return
 
	end subroutine shypart_setup

!*****************************************************************

	subroutine deallocate_array

	implicit none

        return

	end subroutine deallocate_array

!******************************************************************

        subroutine shypart_finalize

        implicit none

        return

        end subroutine shypart_finalize
        

!******************************************************************

        subroutine partPHG

        implicit none

        return

        end subroutine partPHG
        
!==================================================================
        end module shypart
!==================================================================
