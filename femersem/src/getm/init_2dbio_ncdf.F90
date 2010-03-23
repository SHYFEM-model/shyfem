!$Id: $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE:   init_2dbio_ncdf
! 
! !INTERFACE:
    subroutine init_2dbio_ncdf( mode,filename,name,ncben,status, var)
!
! !USES:
    use exceptions, only: getm_error
    use domain,     only: iextr,jextr
    use domain,     only: imin,imax,jmin,jmax,kmax
    use domain,     only: ioff,joff
    use ncdf_topo,  only: ncdf_read_2d
!
! !INPUT PARAMETERS: 
    implicit none

    integer ,        intent(IN)          :: mode
    character(len=*),intent(IN)          :: filename
    character(len=*),intent(IN)          :: name
!
! !OUTPUT PARAMETERS:
    integer,         intent(INOUT)       :: ncben
    integer,         intent(INOUT)       :: status
    REALTYPE,intent(INOUT)               :: var(E2DFIELD)
! !LOCAL  VARIABLES:
    integer                              :: error
    integer                              :: i
    integer                              :: il,ih,jl,jh,ilocl,jlocl,iloch,jloch
    integer                              :: dimlen
    integer                              :: name_id
    integer,dimension(4)                 :: dimidsT
    character(len=80)                    :: text

    include "netcdf.inc"
!
! !DESCRIPTION:  read 2d-field with much less checks as done for bathymetry.
!          
!
!   01-07-2006 Piet Ruardij
! 
!EOP
!-------------------------------------------------------------------------
!BOC




    
    if ( mode == 1 .or.mode == 11 ) then
       if (ncben /=0 ) then
         status=nf_close(ncben)
         if (status .ne. NF_NOERR) then
          text= "Error closing "//trim(filename)//"."
          goto 100
         endif
       endif
       status = nf_open(filename,nf_nowrite,ncben)
       if (status .ne. NF_NOERR) then
          text= "Error opening "//trim(filename)//"."
          goto 100
       endif
    endif

    if ( mode  > 0 ) then

!   Look for name
       status = nf_inq_varid(ncben,name,name_id)
       if (status .ne. NF_NOERR) then
          write(text,'(''Could not find name in '',A,''.'')') trim(filename)
          status=1
          return
       endif
   
!   Is name a matrix?
       status = nf_inq_varndims(ncben,name_id,dimlen)
       if (status .ne. NF_NOERR) then
          write(text,'(''Could not get '',A,'' of '',A,''.'')') '''dimlen''',trim(filename)
          goto 100
       endif

       if (dimlen.lt.2) then
          text="name must have 2 dimensions."
          goto 100
       endif

!   Is the size of name consistent?
       status = nf_inq_vardimid(ncben,name_id,dimidsT)
       if (status .ne. NF_NOERR) then
           write(text,'(''Could not get '',A,'' of '',A,''.'')') 'dimensions',trim(filename)
            goto 100
       endif


       do i=1,dimlen
          status = nf_inq_dimlen(ncben,dimidsT(i),dimlen)
          if (status .ne. NF_NOERR) then
            write(text,'(''Could not get dimlength'',i2,'' of '',A,''.'')') i,trim(filename)
          endif
    
          select case (i) 
            case(1)
                if (dimlen.ne.iextr) error=i
                 !   Get i-dimension for dynamic allocation
            case(2)
                 if (dimlen.ne.jextr) error=i
            case default 
               if (dimlen.ne.1) error=i
         end select 
         if ( error .ne.0 ) then
           write(text,'(''Length of dimension'',i2,'' in'',a,''inconsistent.'')') &
                                       error,name
           goto 100
         endif
       enddo

!  GLOBAL index range for variable to be read
      il    = max(imin+ioff,1);   ih    = min(imax+ioff,iextr)
      jl    = max(jmin+joff,1);   jh    = min(jmax+joff,jextr)
   
!  LOCAL index range for variable to be read
!  (different from GLOBAL range only for parallel runs)
      ilocl = max(imin-ioff,1);   jlocl = max(jmin-joff,1)
      iloch = ih-il+ilocl;        jloch = jh-jl+jlocl;

!  Read bathymetry
      call ncdf_read_2d(ncben,name_id,var(ilocl:iloch,jlocl:jloch),il,ih,jl,jh)
   endif

   if ( mode .gt. 10.or. mode.lt.0 ) status=nf_close(ncben)


    status=0
    return

    100 call getm_error("ncdf_check_grid()", text)   

    end
!EOC
!-----------------------------------------------------------------------
! Copyright by the GOTM-team and BFM-team under the GNU Public License 
! www.gnu.org
!-----------------------------------------------------------------------
