
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

	parameter ( nkndim = 11000 )
	parameter ( neldim = 22000 )
	parameter ( nlvdim = 8 )

	parameter ( mbwdim = 300 )
	parameter ( ngrdim = 12 )

	parameter ( nbcdim = 50 )
	parameter ( nrbdim = 1000 )
	parameter ( nb3dim = 30000 )

	parameter ( nardim = 200 )
	parameter ( nexdim = 100 )
	parameter ( nfxdim = 5000 )
	parameter ( ncsdim = 9 )

	parameter ( nbdydim = 500000 )

!------------------------------------------------------------------------
! do not change anything beyond this line
!------------------------------------------------------------------------

	integer nlkdim			!dimension for side index
        parameter (nlkdim=3*neldim+2*nkndim)

!------------------------------------------------------------------------
! end of parameter file
!------------------------------------------------------------------------

