
c extra data structure for grd files

	integer nlidim
	integer nlndim

	parameter ( nlidim = 100 )	!maximum number of lines
	parameter ( nlndim = nkndim )	!maximum number of nodes in lines

	integer nli			!number of lines read
	integer iplv(nlidim)		!external line numbers

	!integer iarnv(nkndim)		!types for nodes
	integer iarlv(nkndim)		!types for lines

	real hllv(nlidim)		!depth of lines

	integer ipntlv(0:nlidim)	!pointer into line structure
	integer inodlv(nlndim)		!nodes of lines

	common /nli/nli
	common /iplv/iplv
	!common /iarnv/iarnv
	common /iarlv/iarlv
	common /hllv/hllv
	common /ipntlv/ipntlv
	common /inodlv/inodlv

	!save /nli/,/iplv/,/iarnv/,/iarlv/,/hllv/,/ipntlv/,/inodlv/
	save /nli/,/iplv/,/iarlv/,/hllv/,/ipntlv/,/inodlv/

