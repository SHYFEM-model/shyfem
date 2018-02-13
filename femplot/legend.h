
        integer legdim			!maximum number of legend entries
        parameter(legdim=200)

        integer nleg,nlegdi,iplotleg
        common /nleg/ nleg,nlegdi,iplotleg

        real xleg(2,legdim), yleg(2,legdim)
        common /xyleg/ xleg,yleg

        character*80 legleg(legdim)
        common /legleg/ legleg

        character*4 whatleg(legdim)
        common /whatleg/ whatleg

        integer legsiz(legdim)
        common /legsiz/ legsiz

        real aleg(legdim), cleg(legdim)
        common /aleg/ aleg
        common /cleg/ cleg

	save /nleg/,/xyleg/,/legsiz/
	save /legleg/
	save /whatleg/
	save /aleg/,/cleg/

