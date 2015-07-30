
	integer nbc_dim
	parameter ( nbc_dim = 100 )

        character*80 boundn(nbc_dim)
        character*80 conzn(nbc_dim)
        character*80 saltn(nbc_dim)
        character*80 tempn(nbc_dim)
        character*80 bio2dn(nbc_dim)
        character*80 sed2dn(nbc_dim)
        character*80 mud2dn(nbc_dim)
        character*80 lam2dn(nbc_dim)
        character*80 dmf2dn(nbc_dim)
        character*80 tox3dn(nbc_dim)
        character*80 bfm1bc(nbc_dim)
        character*80 bfm2bc(nbc_dim)
        character*80 bfm3bc(nbc_dim)
        character*80 vel3dn(nbc_dim)

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

