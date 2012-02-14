
c extra data structure for grd files

	integer nlidim
	integer nlndim

	parameter ( nlidim = 100 )	!maximum number of lines
	parameter ( nlndim = nkndim )	!maximum number of nodes in lines

	integer iplv(nlidim)		!external line numbers

	integer iarnv(nkndim)		!types for nodes
	integer iarlv(nkndim)		!types for lines

	real hllv(nlidim)		!depth of lines

	integer ipntlv(0:nlidim)	!pointer into line structure
	integer inodlv(nlndim)		!nodes of lines

	common /iplv/iplv
	common /iarnv/iarnv
	common /iarlv/iarlv
	common /hllv/hllv
	common /ipntlv/ipntlv
	common /inodlv/inodlv

	save /iplv/,/iarnv/,/iarlv/,/hllv/,/ipntlv/,/inodlv/

