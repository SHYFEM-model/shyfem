
c*****************************************************************

	subroutine transfer_hlv

	use levels, only : copy_hlv

	implicit none

	include 'param.h'

        integer l

	nlv = nlvdi
        call copy_hlv(nlv,hlv)

        write(6,*) 'hlv copied: ',nlv
        write(6,*) (hlv(l),l=1,nlv)

	end

c*****************************************************************

