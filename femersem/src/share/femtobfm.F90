	subroutine feminit_bio

! !USES:
   use standalone
   use mem  
IMPLICIT NONE

     call init_standalone

    end subroutine feminit_bio
!
!***************************************
        subroutine do_BFM_ECO(b1cn,nbfmv1,b2cn,nbfmv2,b3cn,nbfmv3)
!
! !USES:
   use standalone
   use mem

    IMPLICIT NONE

   include '../../../fem3d/param.h'
   include '../../../fem3d/bfm_common.h'

	integer nbfmv1,nbfmv2,nbfmv3
	integer icall
	data icall /0/
	save icall
 	real b1cn(nlvdim,nkndim,nbfmv1)
 	real b2cn(nlvdim,nkndim,nbfmv2)
 	real b3cn(nlvdim,nkndim,nbfmv3)

	if(icall.eq.0)goto 1
 	call hydro_to_bfm(b1cn,nbfmv1,b2cn,nbfmv2,b3cn,nbfmv3)
1	continue
	icall=1
  	call time_manager	
	call timestepping
 	call bfm_to_hydro(b1cn,nbfmv1,b2cn,nbfmv2,b3cn,nbfmv3)
        end subroutine do_BFM_ECO

!***************************************
