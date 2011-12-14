!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Specialized dummy routine to couple BFM with GOTM
!
! !INTERFACE:
   function D2toD1(BoxNumberX,BoxNumberY) result(BoxNumber)
!
! !DESCRIPTION:
!  This dummy routine, originally conceived to resolve the mapping 
!  between BFM 1D structure and 3D OGCM, is trivial in GOTM/GETM
!  Simply returns the upper level of the benthic system (1)
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi 
!
!EOP
!-----------------------------------------------------------------------
!BOC
   integer :: BoxNumberX,BoxNumberY
   integer :: BoxNumber

   BoxNumber = 1

   return
   end function D2toD1
!EOC
!-----------------------------------------------------------------------


