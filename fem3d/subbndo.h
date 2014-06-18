
c----------------------------------------------------------------------
c data structures for open boundary conditions
c----------------------------------------------------------------------

	integer kbcdim
	integer kopdim
	parameter ( kbcdim = nrbdim )	!total number of open boundary nodes
	parameter ( kopdim = ngrdim )	!maximum number of nodes close to OB

	integer nbndo			!total number of OB nodes
	integer ndebug			!unit number for debug messages

	integer nopnod(kbcdim)		!number of internal nodes close to OB
	integer ibcnod(kbcdim)		!number of boundary
	integer kbcnod(kbcdim)		!number of boundary node
	integer itynod(kbcdim)		!type of boundary

	integer nopnodes(kopdim,kbcdim)	!nodes close to OB

	real xynorm(2,kbcdim)		!normal direction for OB node
	real wopnodes(kopdim,kbcdim)	!weights of nodes close to OB

	common /ibndoc/ nbndo,ndebug,nopnod,ibcnod,kbcnod,itynod,nopnodes
	common /rbndoc/ xynorm,wopnodes

	save /ibndoc/, /rbndoc/

c----------------------------------------------------------------------
c next array is global and can be used to check for open boundary nodes
c----------------------------------------------------------------------
c
c       integer iopbnd(nkndim)          !if >0 pointer into array irv
c                                       !if <0 internal boundary (= -ibc)
c
c----------------------------------------------------------------------
c end of file
c----------------------------------------------------------------------

c possible tests:
c
c iopbnd(k) .ne. 0 		boundary node (external or internal)
c iopbnd(k) .gt. 0 		external boundary node (ibtyp = 1,2)
c iopbnd(k) .lt. 0 		internal boundary node (ibtyp = 3)
 
