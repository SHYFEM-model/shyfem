
c include file if basin (BAS) is read

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real grav,fcor,dcor,dirn,rowass,roluft
        common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft

        character*80 descrr
        common /descrr/ descrr

        real xgv(nkndim), ygv(nkndim)
        real hm3v(3,neldim)
        common /xgv/xgv, /ygv/ygv
        common /hm3v/hm3v

        integer nen3v(3,neldim)
        integer ipev(neldim), ipv(nkndim)
        integer iarv(neldim)
        common /nen3v/nen3v
        common /ipev/ipev, /ipv/ipv
        common /iarv/iarv

