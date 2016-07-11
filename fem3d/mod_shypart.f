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

        use mpi

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

	use basin

	logical b_use_mpi

	integer size

        if(b_use_mpi) then

	  call shympi_init_internal(my_id,n_threads)

        end if

	call shympi_get_status_size_internal(size)

	write(6,*) 'shypart initialized: ',my_id,n_threads
	flush(6)

	end subroutine shypart_init

!*****************************************************************

	subroutine shypart_setup

        use basin
        use mod_geom
        use levels

	implicit none

        integer i,nlkdi
        integer, dimension(:,:), allocatable :: mystruct


        integer nlv_est,nlv_read,nlv_final
        integer nlv_e,nlv_k
        real hmax
        real, allocatable :: hlv_aux(:)
        integer ierr
        integer ounit,error
        character*(20) filename
        character*(20) format_string

        nlkdi = 3*neldi+2*nkndi

        nkn = nkndi
        nel = neldi

c------------------------------------------------------------------
c sanity check
c------------------------------------------------------------------
        call mod_geom_init(nkndi,neldi,ngr)
        call get_hmax_global(hmax)
        if(my_id .eq. 0) then
          write(6,*) 'maximum depth: ',hmax
        end if

        nlv_est = nlv
        call estimate_nlv(nlv_est,hmax)

        if(my_id .eq. 0) then
           write(6,*) 'nlv,nlv_est,nlvdi: ',nlv,nlv_est,nlvdi
        end if
        call check_nlv

        if( nlv > 0 ) then
          allocate(hlv_aux(nlv))
          hlv_aux(1:nlv) = hlv(1:nlv)
          call levels_hlv_init(0)
        end if

        call levels_init(nkn,nel,nlv_est)

        if( nlv > 0 ) then
          hlv(1:nlv) = hlv_aux(1:nlv)
          deallocate(hlv_aux)
        end if

c------------------------------------------------------------------
c levels read in from $levels section
c------------------------------------------------------------------

        call adjust_levels(hmax)        !sets hlv, hldv, nlv, sigma_info, etc.

c------------------------------------------------------------------
c set up layer vectors
c------------------------------------------------------------------

        call mklenk(nlkdi,nkndi,neldi,nen3v,ilinkv,lenkv)

        call mklink(nkndi,ilinkv,lenkv,linkv)

        call mkielt(nkndi,neldi,ilinkv,lenkv,linkv,ieltv)

        call set_ilhv

        call set_last_layer

        call set_ilhkv

        allocate(tempNen3v(3,neldi))
        tempNen3v=nen3v

        allocate(temp_linkv(nlkdi))
        allocate(temp_lenkv(nlkdi))
        allocate(temp_ilinkv(nkndi+1))
        allocate(temp_ieltv(3,neldi))
        allocate(temp_ilhv(neldi))
        temp_linkv=linkv
        temp_lenkv=lenkv
        temp_ilinkv=ilinkv
        temp_ieltv=ieltv
        temp_ilhv=ilhv

        return
 
	end subroutine shypart_setup

!*****************************************************************

	subroutine deallocate_array

        use mod_geom

	implicit none

        deallocate(ilinkv)
        deallocate(linkv)
        deallocate(ieltv)

	end subroutine deallocate_array

!******************************************************************

        subroutine shypart_finalize

        implicit none

        integer ierr

	call MPI_FINALIZE(ierr)

        end subroutine shypart_finalize
        
!==================================================================
        end module shypart
!==================================================================
