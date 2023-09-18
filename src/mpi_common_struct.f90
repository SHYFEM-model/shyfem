
!#######################################################################!
!*******************  start module mpi_common_struct  ******************!
!***********************************************************************!
!#######################################################################!
!------------------------------------------------------------------------
! mpi_common_struct - common structures are declared in this module
! allowing to eliminate the deprecate common block in order to use
! the more flexible allocatable structures 
!------------------------------------------------------------------------

module mpi_common_struct

   use mpi
   integer,save :: n_threads = 1
   integer,save :: my_id = 0
   integer,save :: my_unit = 0

   logical, save :: bmpi = .false.         ! bmpi is true for processes number > 1
   logical, save :: bmpi_debug = .false.
   logical, save :: b_use_mpi = .false.    ! b_use_mpi is true if use mpi

   integer,save :: nkn_global = 0		!total basin
   integer,save :: nel_global = 0          
   integer,save :: nkn_local = 0		!local domain + halo for mpi domain
   integer,save :: nel_local = 0
   integer,save :: nkn_inner = 0		!only proper, no ghost
   integer,save :: nel_inner = 0           ! this domain

   integer, allocatable, save :: sreq(:),rreq(:)
   integer, allocatable, save :: sreq_ut(:),rreq_ut(:)
   integer, allocatable, save :: sreq_vt(:),rreq_vt(:)
   integer, allocatable, save :: tilinkv(:) !total (whole basin) struct of ilinkv
   integer, allocatable, save :: tlinkv(:)  !total (whole basin) struct of linkv
   double precision, allocatable, save :: data_send_d(:,:,:)
   double precision, allocatable, save :: data_recv_d(:,:,:)
   double precision, allocatable, save :: data_send_ut(:,:,:)
   double precision, allocatable, save :: data_recv_ut(:,:,:)
   double precision, allocatable, save :: data_send_vt(:,:,:)
   double precision, allocatable, save :: data_recv_vt(:,:,:)

   integer, allocatable, save, dimension(:) :: allPartAssign
   integer,save,allocatable :: ival(:)
   integer, allocatable, save :: total_ieltv(:,:)

   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: myieltv
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: numberElements
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: numberNodes
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: procNodes
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: numberLevels
   INTEGER(1), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: allNodesAssign
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: univocalNodesAssign
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: scatterNodes

   integer,public,save :: maxelements,maxlevelsproc,totalnodes,maxNodes

end module mpi_common_struct

!------------------------------------------------------------------------
! end of mpi_common_struct module
!------------------------------------------------------------------------
