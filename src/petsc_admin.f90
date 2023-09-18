module petsc_admin

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscpc.h>
#include <petsc/finclude/petscviewer.h>
#include <petsc/finclude/petscvec.h90>

   Mat                :: mat_solver
   Vec                :: rhs, sol
   KSP                :: ksp
   PetscViewer        :: viewerr, matviewer     
   PetscViewer        :: vecviewer2     
   integer,allocatable,dimension(:),save :: rhs_glob

   logical debug     

contains

   subroutine map_rhs(nkn,nkn_glob,r)

     use shympi


#ifdef DEBUGON
     integer r(nkn_glob)
#else
     integer r(nkn)
#endif
     integer i,nkn,nkn_glob

#ifdef DEBUGON
     do i=1, nkn_glob
        r(i) = i-1
     end do
#else
     if(bmpi) then        
       do i=1, nkn
          r(i) = domain%nodes%globalID(i)-1
       end do
     else           
       do i=1, nkn
          r(i) = i-1
       end do
     end if
#endif
    
     return             

   end subroutine map_rhs          

   subroutine petsc_init(nkn_glob,nkn_in)

      use shympi
      use basin

      implicit none

      integer,intent(in) :: nkn_glob, nkn_in
      PetscInt           :: ierr

      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

#ifdef DEBUGON
      if(my_id.eq.0) then
      call MatCreate( PETSC_COMM_SELF, mat_solver, ierr) 
      call VecCreate( PETSC_COMM_SELF, rhs, ierr)
      call MatSetSizes(mat_solver,nkn_glob,nkn_glob,nkn_glob,nkn_glob,ierr )  
      call MatSetType(mat_solver, MATSEQAIJ, ierr)  
      call MatSEQAIJSetPreallocation(mat_solver,ngr+1,PETSC_NULL_INTEGER,ierr)  
      call VecSetSizes(rhs, nkn_glob, nkn_glob,ierr)  
      call VecSetFromOptions(rhs,ierr)
      call VecDuplicate(rhs, sol, ierr)  

      call KSPCreate( PETSC_COMM_SELF, ksp, ierr)
      call KSPSetOperators(ksp, mat_solver, mat_solver, ierr)
      call KSPSetFromOptions(ksp, ierr)

      allocate(rhs_glob(nkn_glob))
      call map_rhs(nkn_in,nkn_glob,rhs_glob)
      end if
#else
      call MatCreate( PETSC_COMM_WORLD, mat_solver, ierr) 
      call VecCreate( PETSC_COMM_WORLD, rhs, ierr)
      call MatSetSizes(mat_solver,nkn_in,nkn_in,nkn_glob,nkn_glob,ierr )  
      call MatSetType(mat_solver, MATMPIAIJ, ierr)  
      call MatMPIAIJSetPreallocation(mat_solver,ngr+1,PETSC_NULL_INTEGER,ngr,PETSC_NULL_INTEGER,ierr)  
      call VecSetSizes(rhs, nkn_in, nkn_glob,ierr)  
      call VecSetFromOptions(rhs,ierr)
      call VecDuplicate(rhs, sol, ierr)  

      call KSPCreate( PETSC_COMM_WORLD, ksp, ierr)
      call KSPSetOperators(ksp, mat_solver, mat_solver, ierr)
      call KSPSetFromOptions(ksp, ierr)

      allocate(rhs_glob(nkn_in))
      call map_rhs(nkn_in,nkn_glob,rhs_glob)  
#endif

  

   end subroutine

   subroutine sol_petsc_debug

      use shympi
      use system_matrix
      use mpi_io_admin
      use timing

      implicit none

      include 'femtime.h'

      integer                   :: i,j
      double precision, pointer :: xx_v(:)
      double precision, allocatable :: rvec2d_glob(:)
      PetscInt          :: ierr, its
      Vec               :: residual
  
      character(len=40) :: fname, fname2,fname3
      character*20 formato
      double precision time1, time2

      if(ln_timing) time1 = shympi_wtime()

      debug = .false.  
      debug = .true.  


      call gather_matrix

      if(ln_timing) then
         time2 = shympi_wtime()
         comm_solver_time = comm_solver_time + time2 - time1
      endif

      if(my_id.eq.0) then
         call MatZeroEntries(mat_solver,ierr)

         do j=1,n_threads
            do i=displs_coo(j)+1,displs_coo(j)+n2_local(j)
               call MatSetValue(mat_solver,i2coo_glob(i)-1,j2coo_glob(i)-1,c2coo_glob(i),ADD_VALUES,ierr)
            end do  
         end do

         call MatAssemblyBegin(mat_solver,MAT_FINAL_ASSEMBLY,ierr)
         call MatAssemblyEnd(mat_solver,MAT_FINAL_ASSEMBLY,ierr)

         allocate(rvec2d_glob(nkn_global))
         call rebuild_2d_nodes(rvec2d,rvec2d_glob)

         do i=1,nkn_global
            call VecSetValue(rhs,rhs_glob(i),rvec2d_glob(i),INSERT_VALUES,ierr)
         end do

         call VecAssemblyBegin(rhs, ierr)  
         call VecAssemblyEnd(rhs, ierr)  
         call KSPSolve(ksp, rhs, sol, ierr)
         call KSPGetIterationNumber(ksp, its, ierr)

         if (debug) then 
            if (mod(niter,120) == 0 .or. niter == 1) then
               call PetscViewerCreate(PETSC_COMM_SELF,viewerr,ierr);
               write(fname,'(a,i0,a,i0,a)')'rhs_petsc_mpi',n_threads,'iter',niter,'.txt'
               call PetscViewerASCIIOpen(PETSC_COMM_SELF,fname,viewerr,ierr)
               call PetscViewerPushFormat(viewerr,PETSC_VIEWER_ASCII_MATLAB,ierr)
               !call PetscViewerPushFormat(viewerr,PETSC_VIEWER_BINARY_MATLAB,ierr)
               call VecView(rhs,viewerr,ierr)  
               call PetscViewerDestroy(viewerr,ierr) 
           
               call PetscViewerCreate(PETSC_COMM_SELF,viewerr,ierr);
               write(fname,'(a,i0,a,i0,a)')'solution_petsc_mpi',n_threads,'iter',niter,'.txt'
               call PetscViewerASCIIOpen(PETSC_COMM_SELF,fname,viewerr,ierr)
               call PetscViewerPushFormat(viewerr,PETSC_VIEWER_ASCII_MATLAB,ierr)
               !call PetscViewerPushFormat(viewerr,PETSC_VIEWER_BINARY_MATLAB,ierr)
               call VecView(sol,viewerr,ierr)  
               call PetscViewerDestroy(viewerr,ierr)  
           
               call PetscViewerCreate(PETSC_COMM_SELF,matviewer,ierr);
               write(fname2,'(a,i0,a,i0,a)')'mat_mpi',n_threads,'_iter',niter,'.txt'
               call PetscViewerASCIIOpen(PETSC_COMM_SELF,fname2,matviewer,ierr)
               call PetscViewerPushFormat(matviewer,PETSC_VIEWER_ASCII_MATLAB,ierr)
               call MatView(mat_solver,matviewer,ierr)  
               call PetscViewerDestroy(matviewer,ierr)   
            end if
         end if

         call VecGetArrayReadF90(sol,xx_v,ierr)
         do i=1,nkn_global
            rvec2d_glob(i) = xx_v(i)     
         end do
         call VecRestoreArrayF90(sol,xx_v,ierr)

         if(ln_timing) then
            time1 = shympi_wtime()
            solver_time = solver_time + time1 - time2
         endif

         !! REDISTRIBUTE SOLUTION CORRECTLY OVER MPI PROCESSES
         call shympi_scatter_2d_nodes_d(rvec2d,rvec2d_glob)

         if(ln_timing) comm_solver_time = comm_solver_time + shympi_wtime() - time1

         deallocate(rvec2d_glob)
      else
         if(ln_timing) time1 = shympi_wtime()
         call rebuild_2d_nodes(rvec2d)
         call shympi_scatter_2d_nodes_d(rvec2d)
         if(ln_timing) comm_solver_time = comm_solver_time + shympi_wtime() - time1
      end if

      return

   end subroutine sol_petsc_debug 

   subroutine sol_petsc

      use shympi
      use system_matrix
      use timing

      implicit none

      include 'femtime.h'

      integer                   :: i,j
      double precision, pointer :: xx_v(:)
      PetscInt          :: ierr, its
      Vec               :: residual
      double precision, allocatable, dimension(:) :: total_raux2d
  
      character(len=40) :: fname, fname2,fname3
      character*20 formato
      double precision time1,time2

      if(ln_timing) time1 = shympi_wtime()

      debug = .false.  

!      if (bmpi)  call shympi_exchange_coo(c2coo,index_coo,n2zero)

      call MatZeroEntries(mat_solver,ierr)

      do i=1,n2zero 
        call MatSetValue(mat_solver,i2coo_pet(i)-1,j2coo_pet(i)-1,c2coo(i),ADD_VALUES,ierr)
      end do  

!      do j=0,n_threads-1
!        call shympi_barrier
!        if(my_id.eq.j) then
!
!          do i=1,n2zero 
!            if(i2coo_pet(i).ne.j2coo_pet(i)) then
!            call MatSetValue(mat_solver,i2coo_pet(i)-1,j2coo_pet(i)-1,c2coo(i),ADD_VALUES,ierr)
!            end if
!          end do  
!
!        end if
!        call MatAssemblyBegin(mat_solver,MAT_FLUSH_ASSEMBLY,ierr)
!        call MatAssemblyEnd(mat_solver,MAT_FLUSH_ASSEMBLY,ierr)
!      end do
!
!      do i=1,n2zero
!         if(i2coo_pet(i).eq.j2coo_pet(i)) then
!            call MatSetValue(mat_solver,i2coo_pet(i)-1,j2coo_pet(i)-1,c2coo(i),INSERT_VALUES,ierr)
!         end if
!      end do


      call MatAssemblyBegin(mat_solver,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(mat_solver,MAT_FINAL_ASSEMBLY,ierr)

      call MPI_Barrier(PETSC_COMM_WORLD,ierr)

      do i=1,nkn_inner
        call VecSetValue(rhs,rhs_glob(i),rvec2d(i),INSERT_VALUES,ierr)
      end do    

      call MPI_Barrier(PETSC_COMM_WORLD,ierr)
        
      call VecAssemblyBegin(rhs, ierr)  
      call VecAssemblyEnd(rhs, ierr)  

      call KSPSolve(ksp, rhs, sol, ierr)
      call KSPGetIterationNumber(ksp, its, ierr)
      !! DEBUGG ivb       
      if (debug) then 
        if (mod(niter,120) == 0 .or. niter == 1) then
        call PetscViewerCreate(PETSC_COMM_WORLD,viewerr,ierr);
        write(fname,'(a,i0,a,i0,a)')'rhs_petsc_mpi',n_threads,'iter',niter,'.txt'
        call PetscViewerASCIIOpen(PETSC_COMM_WORLD,fname,viewerr,ierr)
        call PetscViewerPushFormat(viewerr,PETSC_VIEWER_ASCII_MATLAB,ierr)
        !call PetscViewerPushFormat(viewerr,PETSC_VIEWER_BINARY_MATLAB,ierr)
        call VecView(rhs,viewerr,ierr)  
        call PetscViewerDestroy(viewerr,ierr) 
        !!   
        call PetscViewerCreate(PETSC_COMM_WORLD,viewerr,ierr);
        write(fname,'(a,i0,a,i0,a)')'solution_petsc_mpi',n_threads,'iter',niter,'.txt'
        call PetscViewerASCIIOpen(PETSC_COMM_WORLD,fname,viewerr,ierr)
        call PetscViewerPushFormat(viewerr,PETSC_VIEWER_ASCII_MATLAB,ierr)
        !call PetscViewerPushFormat(viewerr,PETSC_VIEWER_BINARY_MATLAB,ierr)
        call VecView(sol,viewerr,ierr)  
        call PetscViewerDestroy(viewerr,ierr)  
        !!!   
        call PetscViewerCreate(PETSC_COMM_WORLD,matviewer,ierr);
        write(fname2,'(a,i0,a,i0,a)')'mat_mpi',n_threads,'_iter',niter,'.txt'
        call PetscViewerASCIIOpen(PETSC_COMM_WORLD,fname2,matviewer,ierr)
        call PetscViewerPushFormat(matviewer,PETSC_VIEWER_ASCII_MATLAB,ierr)
        call MatView(mat_solver,matviewer,ierr)  
        call PetscViewerDestroy(matviewer,ierr)   
        end if
      end if  

      call VecGetArrayReadF90(sol,xx_v,ierr)

      do i=1,nkn_inner
        rvec2d(i) = xx_v(i)     
        !raux_pet(i) = xx_v(i)     
      end do
      call VecRestoreArrayF90(sol,xx_v,ierr)

      if (shympi_is_master().and. bmpi) then
         if (.not. allocated(total_raux2d)) then       
            allocate(total_raux2d(nkn_global))
         end if      
      end if  

      if(ln_timing) then
         time2 = shympi_wtime()
         solver_time = solver_time + time2 - time1
      endif

      if (bmpi) then      !! MPI threads > 1    
         !! REDISTRIBUTE SOLUTION CORRECTLY OVER MPI PROCESSES
         if(shympi_is_master()) then
            call petsc_gather_2d_nodes_d(rvec2d,total_raux2d)
            call shympi_scatter_2d_nodes_d(rvec2d,total_raux2d)
         else
            call petsc_gather_2d_nodes_d(rvec2d)
            call shympi_scatter_2d_nodes_d(rvec2d)
         end if
      end if

      if(ln_timing)  comm_solver_time = comm_solver_time + shympi_wtime() - time2

      return

   end subroutine sol_petsc        

   subroutine petsc_gather_2d_nodes_d(array,reb_array)

        use shympi
        use basin
        use mpi_common_struct
        use timing

        implicit none

        integer sendbuffer,ierr,i
        integer, dimension(n_threads) :: recvbuffer, displs

        double precision array(nkn_inner)
        double precision,allocatable,dimension(:),optional :: reb_array(:)
        double precision time1

        if(ln_timing) time1 = shympi_wtime()

        if ( .not. bmpi ) return

        if(my_id .eq. 0) then
          if(.not. allocated(reb_array)) then
            allocate(reb_array(nkndi))
          end if
        end if

        sendbuffer = numberNodes(my_id+1)

        displs(1) = 0
        do i=2,n_threads
          displs(i) = displs(i-1) + numberNodes(i-1)
        end do
 
        call MPI_GATHERV(array, sendbuffer, MPI_DOUBLE_PRECISION, &
                      reb_array, numberNodes, displs, &
                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        if(ln_timing) comm_time = comm_time + shympi_wtime() - time1
      
        return

   end subroutine

   subroutine petsc_final

      implicit none
      PetscInt ierr
      
      call MatDestroy(mat_solver, ierr)  
      call VecDestroy(rhs, ierr)  
      call VecDestroy(sol, ierr)  
      call KSPDestroy(ksp, ierr)  
      call PetscFinalize(ierr)  

        

   end subroutine
!
!   subroutine mumps_init(nkn_glob, n2nz,icoo, jcoo)
!
!      implicit none
!      integer nkn_glob,n2nz
!      integer :: icoo(:),jcoo(:)
!
!      return
!
!   end subroutine
!
!   subroutine mumps_solve(iteration,acoo,rvec,rhs_map,raux)
!
!      implicit none
!      integer iteration
!      double precision :: acoo(:),rvec(:),raux(:)
!      integer :: rhs_map
!
!      return
!
!   end subroutine
!
!   subroutine mumps_finalize
!
!      implicit none
!
!      return
!
!   end subroutine
!
!   subroutine mumps_gather_2d_nodes_d(raux2d,total_raux2d)
!
!      implicit none
!      double precision :: raux2d(:)
!      double precision,optional :: total_raux2d(:)
!
!      return
!
!   end subroutine
!
end module petsc_admin
