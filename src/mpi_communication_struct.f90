
!#######################################################################!
!*******************  start module communicationStruc  ******************!
!***********************************************************************!
!#######################################################################!
!------------------------------------------------------------------------
! communicationStruc - common data structures useful for communication
! and the creation of temporary data structures
!------------------------------------------------------------------------

module mpi_communication_struct

   use mpi

   type commItem
      integer prank !contains the process ID to receive from
      integer numItems   ! number of items to recv
      integer, dimension(:), allocatable :: items   ! my local ID of the data to send
      double precision, dimension(:), allocatable :: x  ! my local ID of the data to send
      double precision, dimension(:), allocatable :: y  ! my local ID of the data to send
   end type commItem

   type mapping_info
      integer, public :: univocalID
      integer, public :: domainID
      integer, public :: totalID
      double precision, public :: xmin,xmax
      double precision, public :: ymin,ymax
      ! contains the global ID of the elements belonging to the process
      integer, public, dimension(:), allocatable :: globalID ! contains the global ID of the elements
      integer, public, dimension(:), allocatable :: mapID ! contains the global ID of the elements
      ! contains the rank ID of the process to which belong the global ID
      ! elements with which to communicate
      !integer, public, dimension(:,:), allocatable :: rankID
      !integer, public, dimension(:,:), allocatable :: neighborID
   end type mapping_info

   type domainInfo
      integer exchanges   ! number of processes to which to send
      integer maxItems   ! max number of items to send 
      integer pos
      integer, dimension(:), allocatable :: process
      type (commItem), dimension(:), allocatable :: ele_temp,ele_receive, ele_send, node_send, node_temp  ! process to send
      type (commItem), dimension(:), allocatable :: halo_temp,halo_temp2,halo_send_e,halo_recv_e
      type (commItem), dimension(:), allocatable :: halo_temp3,halo_temp4,halo_send_k,halo_recv_k
      type (mapping_info) :: elems, nodes
   end type domainInfo

   type (domainInfo), SAVE :: domain

   type tvdInfo
     integer send   ! number of processes to which to send
     integer recv   ! number of processes to receive from
     integer maxsend   ! max number of items to send 
     integer maxrecv   ! max number of items to receive
     ! contains the process ID with which to send information
     !integer, dimension(:), allocatable :: psend
     ! contains the process ID to receive from
     !integer, dimension(:), allocatable :: precv
     ! contains contains local identifier of its element upwind or the global
     ! identifier of an element belonging to another domain
     integer,dimension(:,:,:),allocatable :: eleID !number elements dimension
     ! true if the upwind element is mine false otherwise
     integer,dimension(:,:,:),allocatable :: myele !number elements dimension
     type (commItem), dimension(:), allocatable :: relem, selem  ! process to send,receive
   end type tvdInfo

   type (tvdInfo), save :: tvdmpi

   type bounds_info
      integer rankfirst         ! rank of the process that got the first node of
                                !first boundary
      integer iobound           !id open boundary
      integer nob               !number of open boundaries
      integer nnob              !total number nodes open boundaries
      integer nneighbors         !number boundary neighbor nodes in other domains
      integer send,recv         !number of process to communicate with
      !node id for all open boundaries
      integer, dimension(:), allocatable :: niob   !nnob dimension
      integer, dimension(:), allocatable :: bneighbors   !number of neighbors for boundary (nob dimension)
      integer, dimension(:), allocatable :: first  !nnob dimension
      integer, dimension(:), allocatable :: last  !nnob dimension
      integer, dimension(:), allocatable :: neighbors   !nneighbors dimension
      !dimension of bstart and bend equal to nob(number of open boundarier)
      !they represent the start and end point of the nodes in niob for each boundary 
      integer, dimension(:), allocatable :: bstart !nob dimension
      integer, dimension(:), allocatable :: bend   !nob dimension
      !total number nodes for each boundary for all processes
      integer, dimension(:), allocatable :: tnob  ! nob dimension
      type (commItem), dimension(:), allocatable :: rnode, snode  ! process to send,receive
   end type bounds_info

   type (bounds_info), SAVE :: bounds

        
!   type part
!      type (communication_info) myele, mynodes
!      type (sendInfo) mysend
!      type (receiveInfo) myreceive
!   end type part

!   type (mapping_info), SAVE :: myele, mynodes

!   type (part), SAVE :: mypart

   integer, public, allocatable, save, dimension(:,:) :: statuses

!   type(MPI_Status), public, allocatable, save, dimension(:) :: statuses

end module mpi_communication_struct

!------------------------------------------------------------------------
!------------------------------------------------------------------------
