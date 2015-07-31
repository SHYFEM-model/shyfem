
!------------------------------------------------------------------------
! param.h - parameter file for SHYFEM
!------------------------------------------------------------------------

	integer nkndim			!maximum number of nodes
	integer neldim			!maximum number of elements
	integer nlvdim			!maximum number of vertical levels

	integer mbwdim			!maximum bandwidth
	integer ngrdim			!maximum grade of nodes

	integer nbcdim			!maximum number of open boundaries
	integer nrbdim			!maximum number of open boundary nodes
	integer nb3dim			!maximum storage for boundary info

	integer nardim			!maximum area code
	integer nexdim			!maximum number of extra points
	integer nfxdim			!maximum number of flux points
	integer ncsdim			!maximum concentration variables

	integer nbdydim			!maximum particles for lagrange model

	parameter ( nkndim = 1 )
	parameter ( neldim = 1 )
	parameter ( nlvdim = 1 )

	parameter ( mbwdim = 1 )
	parameter ( ngrdim = 1 )

	parameter ( nbcdim = 1 )
	parameter ( nrbdim = 1 )
	parameter ( nb3dim = 1 )

	parameter ( nardim = 1 )
	parameter ( nexdim = 1 )
	parameter ( nfxdim = 1 )
	parameter ( ncsdim = 1 )

	parameter ( nbdydim = 100000 )

!------------------------------------------------------------------------
! do not change anything beyond this line
!------------------------------------------------------------------------

	integer nlkdim			!dimension for side index
        parameter (nlkdim=3*neldim+2*nkndim)

!------------------------------------------------------------------------
! end of parameter file
!------------------------------------------------------------------------

