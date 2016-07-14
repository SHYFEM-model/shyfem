
!#######################################################################!
!*******************  start module Global_Graph_Data  ******************!
!#######################################################################!

module mpi_global_graph_node

  use mpi
  use zoltan
  use basin
  use shypart
   ! Structure to hold graph

   ! MPI variables rankID and processes number
   !integer, save, public :: my_id, n_threads
   !INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: allPartAssign

   type GRAPH_DATA
     integer :: numMyVertices  ! total vertices in in my partition
     integer :: numAllNbors   ! total number of neighbors of my vertices
     integer(ZOLTAN_INT), dimension(:), allocatable :: vertexGID    ! global ID of each of my vertices
     integer, dimension(:),allocatable :: nborIndex    ! nborIndex[i] is location of start of neighbors for vertex i
     integer(ZOLTAN_INT), dimension(:), allocatable :: nborGID      ! nborGIDs[nborIndex[i]] is first neighbor of vertex i
     integer, dimension(:),allocatable :: nborProc     ! process owning each nbor in nborGID
   end type GRAPH_DATA

   type IELTV_GRID
     integer :: totalNeighbor   ! total number of neighbors
     integer, dimension(:), allocatable :: nNeighbor  ! array containing the number of neighbor for each element
   end type IELTV_GRID     

   !! Zoltan data to store in module
   LOGICAL :: changes
   INTEGER(Zoltan_INT) :: numGidEntries, numLidEntries
   INTEGER(Zoltan_INT) :: numImport, numExport
   INTEGER(Zoltan_INT), POINTER, DIMENSION(:) :: importGlobalGids, exportGlobalGids
   INTEGER(Zoltan_INT), POINTER, DIMENSION(:) :: importLocalGids, exportLocalGids
   INTEGER(Zoltan_INT), POINTER, DIMENSION(:) :: importProcs, exportProcs
   INTEGER(Zoltan_INT), POINTER, DIMENSION(:) :: importToPart, exportToPart

   ! Declare a global Mesh data structure. 
   type(GRAPH_DATA) :: graph

contains

!######################################################################!
!***************  start get_number_of_vertices function  **************!
!######################################################################!

  ! Application defined query functions
  function get_number_of_vertices(mygraph, ierr) result (nVert)

    implicit none

    ! in
    integer(Zoltan_INT) :: mygraph(1)

    ! out
    integer, intent(out):: ierr
    integer :: nVert

    ierr = ZOLTAN_OK
    nVert = graph%numMyVertices
  end function get_number_of_vertices

!######################################################################!
!******************  start get_vertex_list subroutine  ****************!
!######################################################################!

  subroutine get_vertex_list(mygraph, sizeGID, sizeLID, globalID, localID, wgt_dim, obj_wgts, ierr)

    use levels

    implicit none

    ! in
    integer(Zoltan_INT), intent(in) :: mygraph(1)
    integer(Zoltan_INT), intent(in) :: sizeGID, sizeLID, wgt_dim
    
    ! out
    integer(Zoltan_INT), intent(out), dimension(*) :: globalID, localID  
    real(Zoltan_FLOAT), dimension(*) :: obj_wgts
    !real(Zoltan_FLOAT), intent(out), dimension(*) :: obj_wgts
    integer(Zoltan_INT), intent(out) :: ierr
    
    ! local
    integer :: i
   
    ierr = ZOLTAN_OK

    ! In this example, return the IDs of our vertices, but no weights.
    ! Zoltan will assume equally weighted vertices.

    do i=1,graph%numMyVertices
      globalID(i) = graph%vertexGID(i)
      localID(i) = i
      obj_wgts(i) = 0.650528 + 0.349472 * ilhkv(globalID(i))
      !obj_wgts(i) = 0.70 + 0.30 * ilhkv(globalID(i))
    end do
  end subroutine get_vertex_list

!######################################################################!
!****************  start get_num_edges_list subroutine  ***************!
!######################################################################!

  subroutine get_num_edges_list(mygraph, sizeGID, sizeLID, num_obj, globalID, localID, numEdges, ierr)

    implicit none
    
    ! in
    integer(Zoltan_INT) :: mygraph(1)
    integer(Zoltan_INT), intent(in) :: sizeGID, sizeLID, num_obj
    integer(Zoltan_INT), dimension(*), intent(in) :: globalID,localID
    
    ! out
    integer(Zoltan_INT), dimension(*), intent(out) :: numEdges
    integer(Zoltan_INT), intent(out) :: ierr

    ! local
    integer :: i, idx

    if ((sizeGID /= 1) .or. (sizeLID /= 1) .or. (num_obj /= graph%numMyVertices)) then
      ierr = ZOLTAN_FATAL
      return
    end if
    
    do i=1,num_obj
      idx = localID(i)
      numEdges(i) = graph%nborIndex(idx+1) - graph%nborIndex(idx)
    end do

    ierr = ZOLTAN_OK
  end subroutine get_num_edges_list

!######################################################################!
!*******************  start get_edge_list subroutine  *****************!
!######################################################################!

  subroutine get_edge_list(mygraph,sizeGID,sizeLID,num_obj,globalID,localID,num_edges,nborGID,nborProc,wgt_dim,ewgts, ierr)

    implicit none

    !in
    integer(Zoltan_INT), intent(in) :: mygraph(1)
    integer(Zoltan_INT), intent (in) :: sizeGID, sizeLID, num_obj, wgt_dim
    integer(Zoltan_INT), dimension(*), intent(in) :: globalID,localID, num_edges
    
    !out
    integer(Zoltan_INT), dimension(*), intent(out) :: nborGID, nborProc
    real(Zoltan_FLOAT), dimension(*) :: ewgts
    integer(Zoltan_INT), intent (out) :: ierr
    
    !local
    integer i, j, from, to

    ierr = ZOLTAN_OK

	if ((sizeGID /= 1).or.(sizeLID /= 1).or.(num_obj /= graph%numMyVertices).or.(wgt_dim /= 0)) then
	    ierr = ZOLTAN_FATAL
	    return
	end if

    do i=1, num_obj

      ! In this example, we are not setting edge weights.
      ! Zoltan will set each edge to weight 1.0.

      to = graph%nborIndex(localID(i)+1)
      from = graph%nborIndex(localID(i))
      if ((to - from) /= num_edges(i)) then
        ierr = ZOLTAN_FATAL
        exit
      end if
      do j=from+1, to
        nborGID(j) = graph%nborGID(j)
        nborProc(j) = graph%nborProc(j)
      end do
    end do
    
  end subroutine get_edge_list

  subroutine mpi_get_param

    use basin
    use mod_geom

    implicit none

    integer param(10)
    integer ierr

    if(my_id .eq. 0) then
      param(1) = nkndi
      param(2) = neldi
      param(3) = ngr
      param(4) = mbw
      param(5) = nkn
    end if

    call MPI_BCAST(param, 5, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)  

    if(my_id .ne. 0) then
      nkndi = param(1)
      neldi = param(2)
      ngr = param(3)
      mbw = param(4)
      nkn = param(5)
    end if 

 

  end subroutine mpi_get_param

!######################################################################!
!********************  end Module Global_Graph_Data  ******************!
!######################################################################!

  end module mpi_global_graph_node
