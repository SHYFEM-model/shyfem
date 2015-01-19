
        character*80 boundn(nbcdim)
        character*80 conzn(nbcdim)
        character*80 saltn(nbcdim)
        character*80 tempn(nbcdim)
        character*80 bio2dn(nbcdim)
        character*80 sed2dn(nbcdim)
        character*80 mud2dn(nbcdim)
        character*80 lam2dn(nbcdim)
        character*80 dmf2dn(nbcdim)
        character*80 tox3dn(nbcdim)
        character*80 bfm1bc(nbcdim)
        character*80 bfm2bc(nbcdim)
        character*80 bfm3bc(nbcdim)
        character*80 vel3dn(nbcdim)

        common /boundn/ boundn
        common /conzn/ conzn
        common /saltn/ saltn
        common /tempn/ tempn
        common /bio2dn/ bio2dn
        common /sed2dn/ sed2dn
        common /mud2dn/ mud2dn
        common /lam2dn/ lam2dn  !!!!!!!!!!!!!!!!! BUG
        common /dmf2dn/ dmf2dn
        common /tox3dn/ tox3dn
        common /bfm1bc/ bfm1bc
        common /bfm2bc/ bfm2bc
        common /bfm3bc/ bfm3bc
        common /vel3dn/ vel3dn

	save /boundn/,/conzn/,/saltn/,/tempn/
	save /bio2dn/,/sed2dn/,/mud2dn/,/lam2dn/,/dmf2dn/,/tox3dn/
	save /bfm1bc/,/bfm2bc/,/bfm3bc/
	save /vel3dn/

