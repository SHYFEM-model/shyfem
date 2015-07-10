
	module mod_nohyd

	implicit none

!        double precision qpnv(nlvdim,nkndim)
!        common /qpnv/qpnv
!        double precision qpov(nlvdim,nkndim)
!        common /qpov/qpov
!        save /qpnv/,/qpov/

        integer, private, save :: nkn_nohyd = 0
        integer, private, save :: nlv_nohyd = 0
        
        double precision, allocatable, save :: qpnv(:,:)
        double precision, allocatable, save :: qpov(:,:)

        contains

*******************************************************************

	subroutine mod_nohyd_init(nkn,nlv)

	integer nkn
        integer nlv

        if( nkn == nkn_nohyd .and. nlv == nlv_nohyd ) return

        if( nkn > 0 .or. nlv > 0 ) then
          if( nkn == 0 .or. nlv == 0 ) then
            write(6,*) 'nkn,nlv: ',nkn,nlv
            stop 'error stop mod_nohyd_init: incompatible params'
          end if
        end if

	if( nkn_nohyd > 0 ) then
         deallocate(qpnv)
	 deallocate(qpov)
        end if

         nkn_nohyd = nkn 
         nlv_nohyd = nlv       
        
        if( nkn == 0 ) return
        
         allocate (qpnv(nlv,nkn))
         allocate (qpov(nlv,nkn))
        
        end subroutine mod_nohyd_init 

!*****************************************************

        end module mod_nohyd

