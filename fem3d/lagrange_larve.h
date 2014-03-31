c setting for larval simulation 

	logical sedeff
	parameter (sedeff=.true.)
	logical sedgrowth
	parameter (sedgrowth=.false.)
	logical sedbug
	parameter (sedbug=.false.)
	logical typsink 
	parameter (typsink=.true.)
	logical foodyn 
	parameter (foodyn=.false.)
	logical labflag
	parameter (labflag=.false.)

c        integer nsedstate
c        parameter( nsedstate = 3 )

c        real larsed(nlvdim,nkndim,nsedstate) ! sediment composition %
c        common /larsed/larsed


c	integer iform 	! funcform of growth dependance 
c	integer lunit 	! unit espress in [um] or [cm]
c	integer itemp 	! funcform of temperature dependance 

	integer icalib 	! combination of funcform dependance 

	real itime  	! output frequency fort.99	
	real rlinit  	! initial larval length [um]
	real Wwinit  	! initial larval weight []

	real p_Lmax    	! maximum length 
	real p_Limp  	! settlement length
	real p_Kl    	! maximum length 
	real p_Kll    	! maximum length 
        integer p_PLD


c 10-20 Von Bert + Lassiter; ==11 VonBert + beta;  ==12 VonBert + lin
c 20-30 Logistic + lassiter; ==21 logisti + beta;  ==22 Logistic + lin
c 30-40 Lineare  + Lassiter; ==31 lineare + beta;  ==32 linear + lin
c 15 = VB + las
c 25 = LG + las
	parameter ( icalib = 15)

	parameter (rlinit = 60) 
	parameter (Wwinit = 0.0001) 


	parameter (p_Lmax = 240)
	parameter (p_Limp  = 220)
	parameter (p_Kl  = 70)
	parameter (p_Kll  = 100)
	parameter (p_PLD = 45*86400)


