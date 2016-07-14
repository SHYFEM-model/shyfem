
module mpi_utility

   use mpi

   contains

!######################################################################!
!*******************  start next_mpi_tag function  ********************!
!######################################################################!
!######################################################################!
!* This function computes the next tag for mpi communications *********!
!######################################################################!

   integer function next_mpi_tag()

      use mpi_communication_struct

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


end module mpi_utility

