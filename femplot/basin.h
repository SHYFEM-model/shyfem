
c include file if basin (BAS) is read

        !integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        !common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        !real grav,fcor,dcor,dirn,rowass,roluft
        !common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft

        integer nkn,nel,ngr,mbw
        common /nbasin/ nkn,nel,ngr,mbw
        save /nbasin/

        integer nkndi,neldi,ngrdi,mbwdi
	parameter ( nkndi = nkndim )
	parameter ( neldi = neldim )
	parameter ( ngrdi = ngrdim )
	parameter ( mbwdi = mbwdim )

        real dcorbas,dirnbas
        common /bkonst/ dcorbas,dirnbas
	save /bkonst/

        character*80 descrr
        common /descrr/ descrr
	save /descrr/

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv
        real hm3v(3,neldim)
        common /hm3v/hm3v
	save /xgv/,/ygv/,/hm3v/

        integer nen3v(3,neldim)
        common /nen3v/nen3v
        integer ipev(neldim), ipv(nkndim)
        common /ipev/ipev, /ipv/ipv
        integer iarv(neldim)
        common /iarv/iarv
        integer iarnv(nkndim)
        common /iarnv/iarnv
	save /nen3v/,/ipev/,/ipv/,/iarv/,/iarnv/

