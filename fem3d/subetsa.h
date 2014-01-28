
	integer nets			!total number of extra nodes
	integer nkets(nexdim)		!node numbers
	real xets(nexdim)		!x-coordinates of nodes
	real yets(nexdim)		!y-coordinates of nodes
	character*80 chets(nexdim)	!description of nodes

	integer ilets(nexdim)		!layers for nodes (needed for header)
	real hets(nexdim)		!depth for nodes (needed for header)

	common /ets_int/ nets,nkets,ilets
	common /ets_real/ xets,yets,hets
	common /ets_char/ chets

	save /ets_int/,/ets_real/,/ets_char/

