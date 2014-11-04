
	!integer nmtdim			!size of meteo regular data
	!parameter (nmtdim = 10000)

	!real mdata(nmtdim)

        real wxv(1),wyv(1)
        common /wxv/wxv,/wyv/wyv
        real ppv(1)
        common /ppv/ppv

        real metrad(1),methum(1)
        real mettair(1),metcc(1)
        common /metrad/metrad, /methum/methum
        common /mettair/mettair, /metcc/metcc

        real metrain(1)
        common /metrain/metrain

        real tauxnv(1),tauynv(1)
        common /tauxnv/tauxnv,/tauynv/tauynv

        real metwbt(1),metws(1)
        common /metwbt/metwbt, /metws/metws

        real windcd(1)          !wave drag coefficient
        common /windcd/windcd

	real evapv(1)		!evaporation
	common /evapv/evapv

c metrain and evapv are in [m/s]
c metrain is read from file in [mm/day] and converted to [m/s]

