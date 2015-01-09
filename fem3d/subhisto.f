c
c $Id: subhisto.f,v 1.1 2004/08/26 15:12:55 georg Exp $
c
c routines dealing with histogram
c
c contents :
c
c subroutine histo_init(nbin,bin0,dbin,rbin)
c subroutine histo_insert(value)
c subroutine histo_final(ic)
c 
c revision log :
c
c 16.03.2004    ggu     routines written from scratch
c
c****************************************************************

        subroutine histo_init(nbin,bin0,dbin,rbin)

        implicit none

        integer nbin            !total number of bins
        real bin0               !first bin (limit)
        real dbin               !regular bin size (0 => use rbin)
        real rbin(nbin)         !bin size limits (upper)

	include 'histo.h'

        integer i

        if( nbin .gt. ndim_histo ) stop 'error stop histo_init: nbin'

        if( dbin .gt. 0 ) then
          do i=1,nbin
            abin(i) = bin0 + (i-1) * dbin
          end do
        else
          do i=1,nbin
            abin(i) = rbin(i)
          end do
        end if

        ncbin = nbin
        do i=1,nbin+1
          icount(i) = 0
        end do

        end

c****************************************************************

        subroutine histo_insert(value)

        implicit none

        real value

	include 'histo.h'

        integer i

        do i=1,ncbin
          if( value .le. abin(i) ) then
            icount(i) = icount(i) + 1
            return
          end if
        end do

        i = ncbin+1
        icount(i) = icount(i) + 1

        end

c****************************************************************

        subroutine histo_final(ic)

        implicit none

        integer ic(1)

	include 'histo.h'

        integer i

        do i=1,ncbin+1
          ic(i) = icount(i)
        end do

        end

c****************************************************************

