
c-------------------------------------------------------------
c header file for ETS file format
c-------------------------------------------------------------
c
c ftype		file id
c maxvers	newest version of file
c maxcomp	compatible function calls down to here
c
c ndim		number of possible entries (open files)
c nitdim	number of integer values to be stored
c nchdim	number of string values to be stored
c
c etsitem	number of maximum entries in table
c
c etsvar	integer parameters of open files
c etschar	string parameters of open files
c
c etsvar(0,n)   iunit
c etsvar(1,n)   nvers
c etsvar(2,n)   nkn
c etsvar(3,n)   not used
c etsvar(4,n)   nlv
c etsvar(5,n)   nvar
c etsvar(6,n)   date
c etsvar(7,n)   time
c
c etschar(1,n)  title
c etschar(2,n)  femver
c
c-------------------------------------------------------------
c parameters
c-------------------------------------------------------------

        integer ftype,maxvers,maxcomp
        parameter(ftype=163,maxvers=1,maxcomp=1)

        integer ndim,nitdim,nchdim
        parameter(ndim=10,nitdim=7,nchdim=2)

c-------------------------------------------------------------
c common
c-------------------------------------------------------------

        integer etsitem
        common /etsitm/etsitem

        integer etsvar(0:nitdim,ndim)
        common /etsvar/etsvar

        character*80 etschar(nchdim,ndim)
        common /etschar/etschar

c-------------------------------------------------------------
c save
c-------------------------------------------------------------

        save /etsitm/
        save /etsvar/
        save /etschar/

c-------------------------------------------------------------
c end of header
c-------------------------------------------------------------

