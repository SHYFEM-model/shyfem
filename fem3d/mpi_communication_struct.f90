
!#######################################################################!
!*******************  start module communicationStruc  ******************!
!***********************************************************************!
!#######################################################################!
!------------------------------------------------------------------------
! communicationStruc - common data structures useful for communication
! and the creation of temporary data structures
!------------------------------------------------------------------------

module communicationStruct

   use mpi

   ! MPI variables rankID and processes number
   integer, save, public :: my_rank, numprocs

   type communication_info
      integer, public :: numberID
      integer, public :: totalID
      ! contains the global ID of the elements belonging to the process
      integer, public, dimension(:), allocatable :: globalID, localID ! contains the global ID of the elements
      ! contains the rank ID of the process to which belong the global ID
      ! elements with which to communicate
      integer, public, dimension(:,:), allocatable :: rankID
      integer, public, dimension(:,:), allocatable :: neighborID
   end type communication_info

   type sendsID
      integer numItems   ! number of items to send
      integer, dimension(:), allocatable :: items   ! my local ID of the data to send
   end type sendsID
   
   type sendInfo
      integer sends   ! number of processes to which to send
      integer maxItems   ! max number of items to send 
      integer, dimension(:), allocatable :: process
      type (sendsID), dimension(:), allocatable :: part_send, all_send, node_send, node_temp  ! process to send
   end type sendInfo

   type receivesID
      integer numItems   ! number of items to receive
      integer, dimension(:), allocatable :: items   ! my local ID of the data to receive
   end type receivesID

   type receiveInfo
      integer receives   ! number of processes from which to receive
      integer maxItems   ! max number of items to receive
      integer, dimension(:), allocatable :: process
      type (receivesID), dimension(:), allocatable :: part_receive, all_receive, node_receive   ! process to receive 
   end type receiveInfo
   
        
   type part
!      type (communication_info) myele, mynodes
      type (sendInfo) mysend
      type (receiveInfo) myreceive
   end type part

   type (communication_info), SAVE :: myele, mynodes, univocal_nodes

   type (part), SAVE :: mypart

   integer, public, allocatable, save, dimension(:) :: sreqUtlnv, rreqUtlnv
   integer, public, allocatable, save, dimension(:) :: sreqVtlnv, rreqVtlnv
   integer, public, allocatable, save, dimension(:) :: sreqWdifhv, rreqWdifhv
   integer, public, allocatable, save, dimension(:) :: statuses
   real, dimension(:,:,:),allocatable :: data_send_utlnv, data_recv_utlnv
   real, dimension(:,:,:),allocatable :: data_send_vtlnv, data_recv_vtlnv
   real, dimension(:,:,:,:),allocatable :: data_send_wdifhv, data_recv_wdifhv

!   type(MPI_Status), public, allocatable, save, dimension(:) :: statuses


end module communicationStruct 

!------------------------------------------------------------------------
! end of commonStructures module
!------------------------------------------------------------------------
