
!-------------------------------------------------------------
! header file for ETS file format
!-------------------------------------------------------------
!
! ftype		file id
! maxvers	newest version of file
! maxcomp	compatible function calls down to here
!
! ndim		number of possible entries (open files)
! nitdim	number of integer values to be stored
! nchdim	number of string values to be stored
!
! etsitem	number of maximum entries in table
!
! etsvar	integer parameters of open files
! etschar	string parameters of open files
!
! etsvar(0,n)   iunit
! etsvar(1,n)   nvers
! etsvar(2,n)   nkn
! etsvar(3,n)   not used
! etsvar(4,n)   nlv
! etsvar(5,n)   nvar
! etsvar(6,n)   date
! etsvar(7,n)   time
!
! etschar(1,n)  title
! etschar(2,n)  femver
!
!-------------------------------------------------------------
! parameters
!-------------------------------------------------------------

        integer ftype,maxvers,maxcomp
        parameter(ftype=163,maxvers=1,maxcomp=1)

        integer ndim,nitdim,nchdim
        parameter(ndim=10,nitdim=7,nchdim=2)

!-------------------------------------------------------------
! common
!-------------------------------------------------------------

        integer etsitem
        common /etsitm/etsitem

        integer etsvar(0:nitdim,ndim)
        common /etsvar/etsvar

        character*80 etschar(nchdim,ndim)
        common /etschar/etschar

!-------------------------------------------------------------
! save
!-------------------------------------------------------------

        save /etsitm/
        save /etsvar/
        save /etschar/

!-------------------------------------------------------------
! end of header
!-------------------------------------------------------------

