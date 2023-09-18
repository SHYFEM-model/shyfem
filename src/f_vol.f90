!
! routine to write volume
!
! revision log :
!
! 28.04.2010    ggu     written from scratch
!
!******************************************************************
!---------------------------------------------------------------------------
        module f_vol
!---------------------------------------------------------------------------
        contains
!---------------------------------------------------------------------------

	subroutine wrfvla

! write of finite volume data

	use layer_thickness
	use area
	use levels
	use basin, only : nkn,nel,ngr,mbw
        use shympi
        use mpi_io_admin
        use output
        use conz_util

	implicit none

	integer k,l,lmax,id,nvar

	double precision getpar
	double precision saux(nlvdi,nkn)

	integer ia_out(4)
	save ia_out

        integer, save :: icall = 0

! start of code

        if( icall .eq. -1 ) return

! initialization

        if( icall .eq. 0 ) then

	  call init_output('itmcon','idtcon',ia_out)

	  if( .not. has_output(ia_out) ) icall = -1
	  if( icall .le. -1 ) return

	  nvar = 1
	  call open_scalar_file(ia_out,nlv,nvar,'fvl')

        end if

! normal call

        icall = icall + 1

	if( .not. next_output(ia_out) ) return

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    saux(l,k) = areakv(l,k) * hdknv(l,k)
	  end do
	end do

        id = 66       			!for finite volume
        !if(bmpi) then
          !call rebuild_3d_nodes(saux, outSaux)
        !  if(shympi_is_master().and.shympi_partition_on_elements()) then
	!    call write_scalar_file(ia_out,id,nlvdi,outSaux)
        !  end if
        !else
	  call write_scalar_file(ia_out,id,nlvdi,saux)
        !end if

	end

!******************************************************************

!---------------------------------------------------------------------------
        end module f_vol
!---------------------------------------------------------------------------
