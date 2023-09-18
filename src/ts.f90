	
! subroutine getts(l,k,t,s)             accessor routine to get T/S
!-------------------------------------------------------------------
        module ts
!-------------------------------------------------------------------

        implicit none

	integer, private, save :: nkn_ts = 0
	integer, private, save :: nlv_ts = 0

        double precision, allocatable, save :: rhov(:,:)
        double precision, allocatable, save :: saltv(:,:)
        double precision, allocatable, save :: tempv(:,:)

        double precision, allocatable, save :: sobsv(:,:)
        double precision, allocatable, save :: tobsv(:,:)
        double precision, allocatable, save :: rtauv(:,:)

        double precision, allocatable, save :: bpresv(:,:)
        double precision, allocatable, save :: bpresxv(:,:)
        double precision, allocatable, save :: bpresyv(:,:)

!-------------------------------------------------------------------
        contains
!-------------------------------------------------------------------

        subroutine mod_ts_init(nkn,nlv,nknex)

        integer nkn
        integer nlv
        integer, optional :: nknex

        if( nkn == nkn_ts .and. nlv == nlv_ts) return

	if( nkn > 0 .or. nlv > 0 ) then
	  if( nkn == 0 .or. nlv == 0 ) then
	    write(6,*) 'nkn,nlv: ',nkn,nlv
	    stop 'error stop mod_ts_init: incompatible parameters'
	  end if
	end if

        if( nkn_ts > 0 ) then
          deallocate(rhov)
          deallocate(saltv)
          deallocate(tempv)
          deallocate(sobsv)
          deallocate(tobsv)
          deallocate(rtauv)
          deallocate(bpresv)
          deallocate(bpresxv)
          deallocate(bpresyv)
        end if

        nkn_ts = nkn
	nlv_ts = nlv

        if( nkn == 0 ) return
          if(present(nknex)) then
            allocate(saltv(nlv,nknex))
            allocate(tempv(nlv,nknex))
          else
            allocate(saltv(nlv,nkn))
            allocate(tempv(nlv,nkn))
          end if

          allocate(rhov(nlv,nkn))
          allocate(sobsv(nlv,nkn))
          allocate(tobsv(nlv,nkn))
          allocate(rtauv(nlv,nkn))
          allocate(bpresv(nlv,nkn))
          allocate(bpresxv(nlv,nkn))
          allocate(bpresyv(nlv,nkn))

        rhov    = 0.d0
        saltv   = 0.d0
        tempv   = 0.d0
        sobsv   = 0.d0
        tobsv   = 0.d0
        rtauv   = 0.d0
        bpresv  = 0.d0
        bpresxv = 0.d0
        bpresyv = 0.d0

        end subroutine mod_ts_init

!*******************************************************************	

	subroutine getts(l,k,t,s)

! accessor routine to get T/S

        implicit none

        integer k,l
        double precision t,s

	include 'param.h'


        t = tempv(l,k)
        s = saltv(l,k)

        end

!************************************************************

        end module ts

