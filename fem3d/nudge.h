
	integer ndgdim			!max size of variables
	parameter (ndgdim = 50)

	integer ndgdatdim		!max size for BC file
	parameter (ndgdatdim = 10*ndgdim)

	real andg_data(ndgdatdim)	!data of observations
	integer ndg_nodelist(ndgdim)	!nodes of obs
	integer ndg_use(ndgdim)		!use observations
	character*40 ndg_names(ndgdim)	!name of stations

	real andg_dist(nkndim)		!distance to be used in algorithm
	real andg_weight(nkndim)	!weight to be used in algorithm
	real andg_obs(nkndim)		!observations to be used in algorithm

	integer ndg_nodes(nkndim)	!nodes of influence
	integer ndg_area(nkndim)	!area of influence

	integer nvar
	real tramp

	common /i_nudge/ nvar,ndg_nodes,ndg_area,ndg_nodelist,ndg_use
	save /i_nudge/

	common /r_nudge/ tramp,andg_data,andg_dist,andg_weight,andg_obs
	save /r_nudge/

	common /c_nudge/ ndg_names
	save /c_nudge/

