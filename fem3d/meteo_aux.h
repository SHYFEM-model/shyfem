
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

	save /metwbt/,/metws/,/metrain/
	save /ppv/,/wxv/,/wyv/,/tauxnv/,/tauynv/
	save /metrad/,/methum/,/mettair/,/metcc/

