!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Specialized dummy routine to couple BFM with GOTM
!
! !INTERFACE:
   function D3toD1(BoxNumberX,BoxNumberY,BoxNumberZ) result(BoxNumber)
!
! !DESCRIPTION:
!  This dummy routine, originally conceived to resolve the mapping 
!  between BFM 1D structure and 3D OGCM, is trivial in GOTM/GETM
!  Simply returns the bottom level of the pelagic system
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
   integer :: BoxNumberX,BoxNumberY,BoxNumberZ
   integer :: BoxNumber

   BoxNumber = BoxNumberZ

   return
   end function D3toD1
!EOC
!-----------------------------------------------------------------------


