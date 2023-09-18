!------------------------------------------------------------------------
module partitioning
!------------------------------------------------------------------------

  use basin
  use util_part
  use para

  contains

!------------------------------------------------------------------------

  subroutine partMetis()  ! MAIN subroutine for partitioning

    implicit none

    ! global variables read by namelist
    integer                     :: nn_wvert     !: weight of the vertices read by namelist
    integer                     :: nn_nparts    !: number of sub-domain to partition the domain read by namelist

    !loop control
    integer je,jke

    !local variables
    integer                     :: iidx

!------------------  Metis struct  -----------------
    integer                     :: iopts(40)            !- metis iopts
    integer                     :: iobjval              !- total communication volume on return
    integer, allocatable        :: ieptr(:)             !- pointer for ieind
    integer, allocatable        :: ieind(:)             !- list of vertices in elements
    integer, allocatable        :: iepart(:)            !- partition of the elements on return
    integer, allocatable        :: inpart(:)            !- partition of the nodes on return
    integer, allocatable        :: ivwgt(:)             !- weights of the vertices
    integer, allocatable        :: ivsize(:)            !- size of communication volume of the nodes
    double precision, pointer   :: ztpwgts(:)=>null()   !- weight for each partition

    character*(20)              :: filename
    character*(20)              :: format_string
    character*(10)              :: what

!C-------------------  allocate Metis struct  ------------------

    allocate(ieptr(nel+1))
    allocate(ieind(3*nel))
    allocate(iepart(nel))
    allocate(inpart(nkn))
    allocate(ilhv(nel))
    allocate(ilhkv(nkn))

!C---------------------  read by namelist  ---------------------
    nn_nparts = nint(getpar('nn_nparts'))
    nn_wvert = nint(getpar('nn_wvert'))

!---------------------------------------------------------------

    call domain_setup

!C--------------------  set Metis options  ---------------------

    write(6,*) '---------------- Calling Metis Partitioning -----------------'

    call METIS_SetDefaultOptions(iopts)

    iopts(1) = 1                !METIS_OPTION_PTYPE (1=kway)
    iopts(2) = 1                !METIS_OPTION_OBJTYPE (1=vol)
    iopts(3) = 1                !METIS_OPTION_CTYPE (1=Sorted heavy-edge matching)
    iopts(6) = 1                !METIS_OPTION_DBGLVL (1=Prints various diagnostic messages)
    iopts(7) = 100              !METIS OPTION NITER (100=number of iterations for the refinement algorithms)
    iopts(8) = 20               !METIS OPTION NCUTS
    iopts(9) = 1               !METIS_OPTION_SEED (selects the seed of the random number generator.)
    iopts(11) = 1              !METIS OPTION MINCONN (Explicitly minimize the maximum connectivity)
    iopts(17) = 1               !METIS_OPTION_UFACTOR (1=1/1000 Specifies the maximum allowed load imbalance among the partitions.)
    iopts(18) = 1               !METIS_OPTION_NUMBERING (1 Fortran-style)

!C-----------------  set Metis input struct  -------------------
!
    ieptr(1)=0
    do je = 1,nel
      do jke = 1,3
        iidx = ieptr(je) + jke 
        ieind(iidx) = nen3v(jke,je)
      end do
      ieptr(je+1) = ieptr(je) + 3
    end do

    allocate(ivwgt(nkn))
    allocate(ivsize(nkn))

    ivwgt=nn_wvert+ilhkv
    ivsize=10+ilhkv

!C-----------------------  call Metis  -------------------------

    call METIS_PartMeshNodal(nel, nkn, ieptr, ieind, ivwgt, ivsize, nn_nparts, ztpwgts, iopts, iobjval, iepart, inpart)

     write(6,*) '----------------- End Metis Partitioning -------------------'

!C-------------   print metis partitioning file   --------------

    if(nn_nparts .gt. 999) then
      format_string = "(A10,I4)"
      write(filename,format_string)'first_part',nn_nparts
    else if(nn_nparts .gt. 99) then
      format_string = "(A10,I3)"
      write(filename,format_string)'first_part',nn_nparts
    else if(nn_nparts .gt. 9) then
      format_string = "(A10,I2)"
      write(filename,format_string)'first_part',nn_nparts
    else
      format_string = "(A10,I1)"
      write(filename,format_string)'first_part',nn_nparts
    end if


    what='elems'

    open(unit=1500, file=filename, action='write')
    write(1500,fmt='(i12,i12,i12,A10)'),nkndi,neldi,nn_nparts,what
    write(1500,fmt="(i12,i12,i12,i12,i12,i12)"),iepart-1
    close(1500)

!C--------------   call refinement routine   -------------------

    call refinement_part(iepart,nn_nparts)

!C------------   print finished partitioning file   ------------

    if(nn_nparts .gt. 999) then
      format_string = "(A10,I4)"
      write(filename,format_string)'part_elems',nn_nparts
    else if(nn_nparts .gt. 99) then
      format_string = "(A10,I3)"
      write(filename,format_string)'part_elems',nn_nparts
    else if(nn_nparts .gt. 9) then
      format_string = "(A10,I2)"
      write(filename,format_string)'part_elems',nn_nparts
    else
      format_string = "(A10,I1)"
      write(filename,format_string)'part_elems',nn_nparts
    end if

    what='elems'

    open(unit=1500, file=filename, action='write')
    write(1500,fmt='(i12,i12,i12,A10)'),nkndi,neldi,nn_nparts,what
    write(1500,fmt="(i12,i12,i12,i12,i12,i12)"),iepart
    close(1500)

    return 
    
  end subroutine

end module

