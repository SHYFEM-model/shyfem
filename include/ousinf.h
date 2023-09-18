
!-------------------------------------------------------------
! header file for ous file format
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
! ousitem	number of maximum entries in table
!
! ousvar	integer parameters of open files
! ouschar	string parameters of open files
!
! ousvar(0,n)	iunit
! ousvar(1,n)	nvers
! ousvar(2,n)	nkn
! ousvar(3,n)	nel
! ousvar(4,n)	nlv
! ousvar(5,n)	not used
! ousvar(6,n)	date
! ousvar(7,n)	time
!
! ousreal(1,n)	href
! ousreal(2,n)	hzmin
!
! ouschar(1,n)	title
! ouschar(2,n)	femver
!
!-------------------------------------------------------------
! parameters
!-------------------------------------------------------------

        integer ftype,maxvers,maxcomp
        parameter(ftype=27,maxvers=2,maxcomp=1)

        integer ndim,nitdim,nrldim,nchdim
        parameter(ndim=10,nitdim=7,nrldim=2,nchdim=2)

!-------------------------------------------------------------
! common
!-------------------------------------------------------------

        integer ousitem
        common /ousitm/ousitem

        integer ousvar(0:nitdim,ndim)
        common /ousvar/ousvar

        double precision ousreal(nrldim,ndim)
        common /ousreal/ousreal

        character*80 ouschar(nchdim,ndim)
        common /ouschar/ouschar

!-------------------------------------------------------------
! save
!-------------------------------------------------------------

        save /ousitm/
        save /ousvar/
        save /ousreal/
        save /ouschar/

!-------------------------------------------------------------
! end of header
!-------------------------------------------------------------

