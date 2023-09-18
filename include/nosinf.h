
!-------------------------------------------------------------
! header file for nos file format
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
! nositem	number of maximum entries in table
!
! nosvar	integer parameters of open files
! noschar	string parameters of open files
!
! nosvar(0,n)   iunit
! nosvar(1,n)   nvers
! nosvar(2,n)   nkn
! nosvar(3,n)   nel
! nosvar(4,n)   nlv
! nosvar(5,n)   nvar
! nosvar(6,n)   date
! nosvar(7,n)   time
!
! noschar(1,n)  title
! noschar(2,n)  femver
!
!-------------------------------------------------------------
! parameters
!-------------------------------------------------------------

        integer ftype,maxvers,maxcomp
        parameter(ftype=161,maxvers=5,maxcomp=3)

        integer ndim,nitdim,nchdim
        parameter(ndim=30,nitdim=7,nchdim=2)

!-------------------------------------------------------------
! common
!-------------------------------------------------------------

        integer nositem
        common /nositm/nositem

        integer nosvar(0:nitdim,ndim)
        common /nosvar/nosvar

        character*80 noschar(nchdim,ndim)
        common /noschar/noschar

!-------------------------------------------------------------
! save
!-------------------------------------------------------------

        save /nositm/
        save /nosvar/
        save /noschar/

!-------------------------------------------------------------
! end of header
!-------------------------------------------------------------

