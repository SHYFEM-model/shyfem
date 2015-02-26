
	logical bdebug_out
	parameter (bdebug_out = .false.)

	integer nunit

	integer iformat,iwave
        common /supout1/ iformat,iwave

	integer nunit_wave,nunit_ous,nunit_nos,nunit_fvl
     +			,nunit_eos,nunit_fem
        common /supout2/ nunit_wave,nunit_ous,nunit_nos,nunit_fvl
     +                  ,nunit_eos,nunit_fem

	real regp(7)
	common /supout3/regp

	save /supout1/,/supout2/,/supout3/

