
	integer maxlnk
	parameter( maxlnk = 100 )

        integer ilinkv(nkndim+1)
        common /ilinkv/ ilinkv
        integer lenkv(nlkdim)
        common /lenkv/ lenkv
        integer lenkiiv(nlkdim)
        common /lenkiiv/ lenkiiv
        integer linkv(nlkdim)
        common /linkv/ linkv

        !integer ilinkv(1)
        !common /ilinkv/ ilinkv
        !integer lenkv(1)
        !common /lenkv/ lenkv
        !integer linkv(1)
        !common /linkv/ linkv

	integer lnk_nodes(maxlnk)
	common /lnk_nodes/ lnk_nodes

	integer lnk_elems(maxlnk)
	common /lnk_elems/ lnk_elems

	save /ilinkv/,/lenkv/,/linkv/,/lenkiiv/
	save /lnk_nodes/,/lnk_elems/

