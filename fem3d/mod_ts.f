	
        module mod_ts

        implicit none

	integer, private, save :: nkn_ts = 0
	integer, private, save :: nlv_ts = 0

        real, allocatable, save :: rhov(:,:)
        real, allocatable, save :: saltv(:,:)
        real, allocatable, save :: tempv(:,:)

        real, allocatable, save :: sobsv(:,:)
        real, allocatable, save :: tobsv(:,:)
        real, allocatable, save :: rtauv(:,:)

        real, allocatable, save :: bpresv(:,:)
        real, allocatable, save :: bpresxv(:,:)
        real, allocatable, save :: bpresyv(:,:)

        contains

!************************************************************

        subroutine mod_ts_init(nkn,nlv)

        integer nkn
        integer nlv

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

          allocate(rhov(nlv,nkn))
          allocate(saltv(nlv,nkn))
          allocate(tempv(nlv,nkn))
          allocate(sobsv(nlv,nkn))
          allocate(tobsv(nlv,nkn))
          allocate(rtauv(nlv,nkn))
          allocate(bpresv(nlv,nkn))
          allocate(bpresxv(nlv,nkn))
          allocate(bpresyv(nlv,nkn))

	rhov = 0.
	saltv = 0.
	tempv = 0.
	sobsv = 0.
	tobsv = 0.
	rtauv = 0.
	bpresv = 0.
	bpresxv = 0.
	bpresyv = 0.

        end subroutine mod_ts_init

!************************************************************

        end module mod_ts

