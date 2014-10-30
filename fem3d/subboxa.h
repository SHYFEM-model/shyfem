
	integer nfxboxdim
	parameter(nfxboxdim=2000)	!total number of nodes in sections

        integer nscboxdim
        parameter(nscboxdim=140)        !maximum number of sections

        integer nbxdim
        parameter(nbxdim=50)            !maximum number of boxes

c	nbox		total number of boxes
c	nsect		total number of sections
c	nbc_ob		total number of open boundaries
c	kfluxm		total number of nodes in section
c	kflux(i)	node numbers defining sections
c
c	iflux(3,i)	aux values for nodes in sections
c	iboxes(ie)	index from elements to boxes
c	ikboxes(k)	index from nodes to boxes (<0 if node on box boundary)
c	isects(4,is)	description of section (n,ipnt,box1,box2)
c	iscbnd(4,ibc)	description of OB (n,type,box1,box2)
c
c	i runs on list of nodes
c	k runs on nodes (1,nkn)
c	ie runs on elements (1,nel)
c	is runs on sections (1,nsect)
c	ib runs on boxes (1,nbox)
c	ibc runs on open boundaries (1,nbc)

        integer nbox,nsect,nbc_ob,kfluxm,kflux(nfxboxdim)
        common /kfluxboxc/ nbox,nsect,nbc_ob,kfluxm,kflux
	save /kfluxboxc/

        integer iflux(3,nfxboxdim)
        common /ifluxbox/iflux
	save /ifluxbox/

	integer iboxes(neldim)
	common /iboxes/iboxes
	save /iboxes/

	integer ikboxes(nkndim)
	common /ikboxes/ikboxes
	save /ikboxes/

	integer isects(4,nscboxdim)	!internal section description
	common /isects/isects
	save /isects/

	integer iscbnd(4,nbcdim)	!boundary condition descriptions
	common /iscbnd/iscbnd
	save /iscbnd/

