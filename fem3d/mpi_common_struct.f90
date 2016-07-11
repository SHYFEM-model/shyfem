
!#######################################################################!
!*******************  start module commonStructures  ******************!
!***********************************************************************!
!#######################################################################!
!------------------------------------------------------------------------
! commonStructure - common structures are declared in this module
! allowing to eliminate the deprecate common block in order to use
! the more flexible allocatable structures 
!------------------------------------------------------------------------

module commonStructures

   !REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: fxv, fyv
   !REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: fcorv
   !REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: rfricv
   !INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: nen3v
   !INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: iarv
   !INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: ilinkv,linkv,lenkv
   !REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: hm3v
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: myieltv
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: numberElements
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: numberNodes
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: numberLevels
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: allNodesAssign
   !REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: uprv,vprv
   !REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: hdenv
   !REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: hdknv
   !REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: hdkov
   !REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: tempv
   !REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: saltv
   !REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: const3d
   !REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: rtauv
   !REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: wprv
   !REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: wlnv
   !REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: wlov
   !REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: wdifhv
   !INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: iopbnd

   integer, public, save :: maxelements, maxlevelsproc, totalnodes

end module commonStructures

!------------------------------------------------------------------------
! end of commonStructures module
!------------------------------------------------------------------------
