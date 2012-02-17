	subroutine feminit_bio(b1cn,nbfmv1,b2cn,nbfmv2,b2cn_a,b2cn_b,b2cn_c,b2cn_d,b3cn,nbfmv3,b3cn_a,b3cn_b,b3cn_c)

! !USES:
   use standalone
   use mem  
IMPLICIT NONE
	include '../../../fem3d/param.h'
        include '../../../fem3d/bfm_common.h'

	integer nbfmv1          !number of solutes transported var 
        integer nbfmv2          !number of fitozoo transported var
        integer nbfmv3          !number of essudates transported var




        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nlvdi,nlv
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /level/ nlvdi,nlv
	 
        real b1cn(nlvdim,nkndim,nbfmv1)
        real b2cn(nlvdim,nkndim,nbfmv2)
        real b2cn_a(nlvdim,nkndim,nbfmv2)
        real b2cn_b(nlvdim,nkndim,nbfmv2)
        real b2cn_c(nlvdim,nkndim,nbfmv2)
        real b2cn_d(nlvdim,nkndim,nbfmv2)

        real b3cn(nlvdim,nkndim,nbfmv3)
        real b3cn_a(nlvdim,nkndim,nbfmv3)
        real b3cn_b(nlvdim,nkndim,nbfmv3)
        real b3cn_c(nlvdim,nkndim,nbfmv3)

! 	function x restart 

	logical has_restart



	integer k

        call init_standalone

	if( .not. has_restart(6) ) then
	 do k =1,nkn
	  call bfm_to_hydro(k,b1cn,nbfmv1,b2cn,nbfmv2,b2cn_a,b2cn_b,b2cn_c,b2cn_d,b3cn,nbfmv3,b3cn_a,b3cn_b,b3cn_c)
    	 end do
	else
	 do k =1,nkn
          call rst_bfm_to_hydro(k,b1cn,nbfmv1,b2cn,nbfmv2,b2cn_a,b2cn_b,b2cn_c,b2cn_d,b3cn,nbfmv3,b3cn_a,b3cn_b,b3cn_c)
         end do
	endif
	
	end subroutine feminit_bio
!


!***************************************



        subroutine do_BFM_ECO(k,b1cn,nbfmv1,b2cn,nbfmv2,b2cn_a,b2cn_b,b2cn_c,b2cn_d,b3cn,nbfmv3,b3cn_a,b3cn_b,b3cn_c)
!
! !USES:
   use standalone
   use mem

    IMPLICIT NONE

   include '../../../fem3d/param.h'
   include '../../../fem3d/bfm_common.h'

	integer nbfmv1          !number of solutes transported var 
        integer nbfmv2          !number of fitozoo transported var
        integer nbfmv3          !number of essudates transported var

	integer icall
	data icall /0/
	save icall

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nlvdi,nlv
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /level/ nlvdi,nlv

	integer node 
	common /node/node
	
	real b1cn(nlvdim,nkndim,nbfmv1)
        real b2cn(nlvdim,nkndim,nbfmv2)
        real b2cn_a(nlvdim,nkndim,nbfmv2)
        real b2cn_b(nlvdim,nkndim,nbfmv2)
        real b2cn_c(nlvdim,nkndim,nbfmv2)
        real b2cn_d(nlvdim,nkndim,nbfmv2)

        real b3cn(nlvdim,nkndim,nbfmv3)
        real b3cn_a(nlvdim,nkndim,nbfmv3)
        real b3cn_b(nlvdim,nkndim,nbfmv3)
	real b3cn_c(nlvdim,nkndim,nbfmv3)

	
	integer k



  	 node =k

   	 call hydro_to_bfm(k,b1cn,nbfmv1,b2cn,nbfmv2,b2cn_a,b2cn_b,b2cn_c,b2cn_d,b3cn,nbfmv3,b3cn_a,b3cn_b,b3cn_c)

!   	call time_manager
	call timestepping
 	call bfm_to_hydro(k,b1cn,nbfmv1,b2cn,nbfmv2,b2cn_a,b2cn_b,b2cn_c,b2cn_d,b3cn,nbfmv3,b3cn_a,b3cn_b,b3cn_c)

        end subroutine do_BFM_ECO

!***************************************
