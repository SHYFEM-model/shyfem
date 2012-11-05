
c-------------------------------------------------------------
c header file for nos file format
c-------------------------------------------------------------

c maxvers	newest version of file
c maxcomp	compatible function calls down to here

c-------------------------------------------------------------
c parameters
c-------------------------------------------------------------

        integer ftype,maxvers,maxcomp
        parameter(ftype=161,maxvers=4,maxcomp=3)

        integer ndim,nitdim
        parameter(ndim=30,nitdim=7)

c-------------------------------------------------------------
c common
c-------------------------------------------------------------

        integer nositem,nosvar(0:nitdim,ndim)
        common /nosvar/nositem,nosvar

        character*80 nostitle(ndim)
        common /nostit/nostitle

c-------------------------------------------------------------
c save
c-------------------------------------------------------------

        save /nosvar/
        save /nostit/

c-------------------------------------------------------------
c end of header
c-------------------------------------------------------------

