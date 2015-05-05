
c----------------------------------------------------------------------
c statement functions to test for nodes with boundary condition
c----------------------------------------------------------------------
c
c iopbnd(k) = 0		no open BC
c iopbnd(k) > 0		external open BC (ibtyp=1,2)
c iopbnd(k) < 0		internal open BC (ibtyp=3)
c
c----------------------------------------------------------------------

	integer iopbnd(nkndim)
	common /iopbnd/iopbnd
	save /iopbnd/

	integer k_n
	logical is_boundary, is_external_boundary, is_internal_boundary
	logical is_inner

	is_boundary(k_n) = iopbnd(k_n) .ne. 0
	is_external_boundary(k_n) = iopbnd(k_n) .gt. 0
	is_internal_boundary(k_n) = iopbnd(k_n) .lt. 0

	is_inner(k_n) = iopbnd(k_n) .eq. 0

c----------------------------------------------------------------------
c end of file
c----------------------------------------------------------------------

