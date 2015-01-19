
        integer ierv(2,nrbdim)
        common /ierv/ierv
	save /ierv/

        real rhv(nrbdim), rlv(nrbdim)
        common /rhv/rhv, /rlv/rlv
        real rrv(nrbdim)
        integer irv(nrbdim)
        common /rrv/rrv, /irv/irv
	save /rhv/,/rlv/,/rrv/,/irv/

        integer iopbnd(nkndim)
        common /iopbnd/iopbnd
	save /iopbnd/

