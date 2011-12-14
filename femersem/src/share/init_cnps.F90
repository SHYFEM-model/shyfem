#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !SUBROUTINE: Initialise variable components
!
! !INTERFACE:
   subroutine init_cnps(c,n,p,s,nc,pc,sc)
!
! !DESCRIPTION:
!  This subroutine initialises the other internal components
!  of biogeochemical variables
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
    REALTYPE,dimension(:),intent(in)           :: c
    REALTYPE,intent(in),optional               :: nc,pc,sc
!
! !OUTPUT PARAMETERS:
    REALTYPE,dimension(:),intent(out),optional :: n
    REALTYPE,dimension(:),intent(out),optional :: p
    REALTYPE,dimension(:),intent(out),optional :: s
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!
!LOCAL VARIABLES:
    REALTYPE                     :: nc_ratio,pc_ratio,sc_ratio
!
!EOP
!-----------------------------------------------------------------------
!BOC

    if (present(nc)) then
      nc_ratio = nc
    else
      nc_ratio = 0.0126 ! Redfield
    end if

    if (present(pc)) then
      pc_ratio = pc
    else
      pc_ratio = 0.7862e-3 ! Redfield
    end if

    if (present(sc)) then
      sc_ratio = sc
    else
      sc_ratio = 0.0145 ! Redfield
    end if

    if (present(n)) n = nc_ratio*c
    if (present(p)) p = pc_ratio*c
    if (present(s)) s = sc_ratio*c

  end subroutine init_cnps
!EOC



!-----------------------------------------------------------------------

