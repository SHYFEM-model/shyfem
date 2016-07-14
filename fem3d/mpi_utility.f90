
module MPI_Utility

   use mpi

   contains

!######################################################################!
!*******************  start next_mpi_tag function  ********************!
!######################################################################!
!######################################################################!
!* This function computes the next tag for mpi communications *********!
!######################################################################!

   integer function next_mpi_tag()

      use communicationStruct

      implicit none

      integer, save::last_tag=0, tag_ub=0


!      if(tag_ub==0) then
!        call MPI_Attr_get(MPI_COMM_WORLD, MPI_TAG_UB, tag_ub, flag, ierr)
!      end if

      last_tag = last_tag + numprocs + 1

!      last_tag = mod(last_tag+1, tag_ub)
!      if(last_tag==0) then
!        last_tag = last_tag+1
!      end if

      next_mpi_tag = last_tag

   end function next_mpi_tag


end module MPI_Utility


!######################################################################!
!************* start pending_communication_communicator  **************!
!######################################################################!
!######################################################################!
!* Return whether there is a pending communication for the supplied   *!
!* communicator.                                                      *!
!######################################################################!

  function pending_communication_communicator(communicator) result(pending)

    integer, optional, intent(in) :: communicator

    logical :: pending

    integer :: lcommunicator

    integer :: ierr, ipending

    if(present(communicator)) then
      lcommunicator = communicator
    else
      lcommunicator = MPI_COMM_WORLD
    end if

    call mpi_iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, lcommunicator, ipending, MPI_STATUS_IGNORE, ierr)
    if(ierr .ne. MPI_SUCCESS) then
       call abort
    end if

    pending = (ipending /= 0)

    ! Note - removing this mpi_barrier could result in a false
    ! positive on another process.
    call mpi_barrier(lcommunicator, ierr)
    if(ierr .ne. MPI_SUCCESS) then
       call abort
    end if

  end function pending_communication_communicator

