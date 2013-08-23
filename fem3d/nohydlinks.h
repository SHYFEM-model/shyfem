

c maxdimnh = maximum dimension of the hypothetic general matrix row that has to be solved for 
c            non hydrostatic contibution
c kpl(0:maxdimnh) = array with nodes assigned in the position on row of general matrix where 
c                   there is a contribution. 
c lpl(0:maxdimnh) = array with layers assigned in the position on row of general matrix where 
c                   there is a contribution. 
c rownh(0:maxdimnh) = array of mask values for row of general (1 if in the position on row of 
c                     general matrix there is a contribution, 0 otherwise
c kpl1(0:maxdimnh) = array with maximum dimension corresponding to the number of positions on 
c                    row of the general matrix where there is a contribution. The array stores 
c                    the corresponding nodes.
c lpl1(0:maxdimnh) = array with maximum dimension corresponding to the number of positions on row 
c                    of the general matrix where there is a contribution. The array stores the 
c                    corresponding layers.
c ipl(0:nkndim) = auxiliary array that stores for each node k, the sum of lmax for all the nodes 
c       till k. It help in the assignment of pointer position in the final general array for the
c       solution of the general matrix
c iop(nkndim,nlvdim,ngrdim) =  given the row corresponding to (k,l) of the general matrix, 
c                              for each nn that 
c                corresponds to the contribution of the node k1 at layer l, connected with k,
c                iop provides the relative position on the row, respect to the diagonal
c iopup(nkndim,nlvdim) = given the row corresponding to (k,l) of the general matrix, iopup provides
c              the sum of all the positions with non zero values on the previous rows of
c              the general matrix until row (k,l)
c kop(nkndim,nlvdim,ngrdim) = for each row (k,l), array with nodes k1 of column (k1,l) contributing
c                             at non zero values in general matrix
c lop(nkndim,nlvdim,ngrdim) = for each row (k,l), array with layer l of column (k1,l) contributing 
c                             at non zero values in general matrix
c iopp(neldim,9,nlvdim) = pointer that identifies the position in the general array that is built to
c                         solve the non hydrostatic contribution, given the considered element ie, 
c                         the 9 contribution of the combination of its nodes (kn(1:3)=nen3v(ie,1:3)
c                         kn(1:3)<-->kn(1:3)
c ioii((ngrdim3*maxdimnh) = node position k on row of the general matrix corresponding to the pointer
c                          iopp
c iojj((ngrdim3*maxdimnh) = node position k1 on column of the general matrix corresponding to the 
c	                   pointer iopp
c ioll((ngrdim3*maxdimnh) = layer position l on row of the general matrix corresponding to the 
c			   pointer iopp
c ioii1((ngrdim3*maxdimnh) = position on row of the general matrix corresponding to the pointer
c                          iopp
c iojj1((ngrdim3*maxdimnh) = position on the column of the general  matrix corresponding to the 
c	                   pointer iopp


	integer csrdimnh
	parameter ( csrdimnh = 9 * neldim * nlvdim )

	integer matdimmax
	common /matdimmax/matdimmax
	save /matdimmax/

	integer nnzeronh
	common /nnzeronh/nnzeronh
	save /nnzeronh/

	integer maxdimnh
	parameter( maxdimnh = nlvdim * nkndim ) 

	integer ngrdim3
	parameter (ngrdim3 = ngrdim + 3 )

	integer kpl(0:maxdimnh)  
	integer lpl(0:maxdimnh)
	integer rownh(0:maxdimnh)
	integer kpl1(0:maxdimnh) 
	integer lpl1(0:maxdimnh)
	integer ipl(0:nkndim)
	integer iop(ngrdim,nlvdim,nkndim) 
        integer iopup(0:nlvdim,nkndim)
	integer iopp(neldim,9,nlvdim)	
	integer ikop(ngrdim3,nlvdim,nkndim)
	integer ilop(ngrdim3,nlvdim,nkndim)
	integer ioii(ngrdim3*maxdimnh)
	integer iojj(ngrdim3*maxdimnh)
	integer ioii1(ngrdim3*maxdimnh)
	integer iojj1(ngrdim3*maxdimnh)
	integer ioll(ngrdim3*maxdimnh)
	double precision conh(ngrdim3*maxdimnh)
	double precision vecnot(ngrdim3*maxdimnh)

	common /conh/conh
	common /vecnot/vecnot
	save /conh/
	save /vecnot/

	common /iopp/iopp,/ioii/ioii,/iojj/iojj,/ioll/ioll
	common /ioii1/ioii1,/iojj1/iojj1

	save /iopp/,/ioii/,/iojj/,/ioll/
	save /ioii1/,/iojj1/


	!double precision rvecnh(maxdimnh) 
	real*8 rvecnh(maxdimnh) 
	common /rvecnh/rvecnh
	!double precision rauxnh(maxdimnh)!controllare dimensione
	real*8 rauxnh(maxdimnh)!controllare dimensione
	common /rauxnh/rauxnh

