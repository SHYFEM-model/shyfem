
        real wxv(nkndim),wyv(nkndim)
        common /wxv/wxv,/wyv/wyv
        real ppv(nkndim)
        common /ppv/ppv

        real metrad(nkndim),methum(nkndim)
        real mettair(nkndim),metcc(nkndim)
        common /metrad/metrad, /methum/methum
        common /mettair/mettair, /metcc/metcc

        real metrain(nkndim)
        common /metrain/metrain

        real tauxnv(nkndim),tauynv(nkndim)
        common /tauxnv/tauxnv,/tauynv/tauynv

        real metwbt(nkndim),metws(nkndim)
        common /metwbt/metwbt, /metws/metws

        real windcd(nkndim)          !wave drag coefficient
        common /windcd/windcd

        real metice(nkndim)          !ice cover
        common /metice/metice
	save /metice/

	real evapv(nkndim)		!evaporation
	common /evapv/evapv

c metrain and evapv are in [m/s]
c metrain is read from file in [mm/day] and converted to [m/s]

