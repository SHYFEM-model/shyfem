module timing

!----------------------------------------------------------------------------

   use para
   use mpi_common_struct

   implicit none
   private

   public timing_init
   public timing_finalize
   !public timing_start
   !public timing_stop

   logical, public,save  :: ln_timing = .false.
   double precision,public, save :: comm_time = 0
   double precision,public, save :: comm_scalar_time = 0
   double precision,public, save :: comm_solver_time = 0
   double precision,public, save :: comm_transp_time = 0
   !double precision,public, save :: comm_ts_time = 0
   double precision,public, save :: transp_time = 0
   double precision,public, save :: solver_time = 0
   double precision, public, save :: scalar_time = 0
   double precision,public, save :: io_time = 0
   double precision,public, save :: io_scalar_time = 0
   double precision,public, save :: io_transp_time = 0

!----------------------------------------------------------------------------
   contains
!----------------------------------------------------------------------------

        subroutine timing_init()

          implicit none

          ln_timing   = getpar('ln_timing') .gt. 0

          comm_time = 0
          comm_scalar_time = 0
          comm_solver_time = 0
          comm_transp_time = 0
          !comm_ts_time = 0
          transp_time = 0
          solver_time = 0
          scalar_time = 0
          io_time = 0
          io_scalar_time = 0
          io_transp_time = 0

          return

        end subroutine timing_init

!----------------------------------------------------------------------------

        subroutine timing_finalize

        implicit none
        integer i,ierr

        do i=0,n_threads-1
          call MPI_Barrier(MPI_COMM_WORLD,ierr)
          if (my_id.eq.i) then
            write(6,*)'TOTAL-TRANSPORT TIME: ',transp_time
            write(6,*)'I/O-TRANSPORT TIME: ',io_transp_time
            write(6,*)'COMMUNICATION-TRANSPORT TIME: ',comm_transp_time
            write(6,*)'TOTAL-SOLVER TIME: ',solver_time+comm_solver_time
            write(6,*)'COMMUNICATION-SOLVER TIME: ',comm_solver_time
            write(6,*)'TOTAL-SCALAR TIME: ',scalar_time
            write(6,*)'I/O-SCALAR TIME: ',io_scalar_time
            write(6,*)'COMMUNICATION-SCALAR TIME: ',comm_scalar_time
            !write(6,*)'COMMUNICATION-TS TIME: ',comm_ts_time
            write(6,*)'COMMUNICATION + DISPLACEMENT TIME: ',comm_time
            write(6,*)'I/O time: ',io_time
            write(6,*)'TRANSPORT TIME:',transp_time-io_transp_time-comm_transp_time
            write(6,*)'SOLVER TIME: ',solver_time
            write(6,*)'SCALAR TIME: ',scalar_time-io_scalar_time-comm_scalar_time
          end if
        end do

        return

        end subroutine timing_finalize

!----------------------------------------------------------------------------

end module
