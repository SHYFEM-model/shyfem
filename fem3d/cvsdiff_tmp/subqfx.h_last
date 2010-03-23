
c-----------------------------------------------------------
c integer variables
c-----------------------------------------------------------

c	ifunit			unit of flux file
c	itfold,itfnew,itfact	time for old/new/actual level
c	itperiod		periodic heat flux condition 
c				(0 if none, default)
c				(31536000 for one year - 365 days)
c	irhumid			relative humidity given 
c				(0 for wet bulb given)
c				(1 for relative humidity, default)

        integer ifunit,itfold,itfnew,itfact,itperiod,irhumid
        common /qflxi/ ifunit,itfold,itfnew,itfact,itperiod,irhumid

c-----------------------------------------------------------
c values of old time level
c-----------------------------------------------------------

        real qsold,taold,tbold,uwold,ccold
     +			,urold,pold,eold,rold,qold
        common /qflxro/ qsold,taold,tbold,uwold,ccold
     +			,urold,pold,eold,rold,qold

c-----------------------------------------------------------
c values of new time level
c-----------------------------------------------------------

        real qsnew,tanew,tbnew,uwnew,ccnew
     +			,urnew,pnew,enew,rnew,qnew
        common /qflxrn/ qsnew,tanew,tbnew,uwnew,ccnew
     +			,urnew,pnew,enew,rnew,qnew

c-----------------------------------------------------------
c values of actual time level
c-----------------------------------------------------------

        real qsact,taact,tbact,uwact,ccact
     +			,uract,pact,eact,ract,qact
        common /qflxra/ qsact,taact,tbact,uwact,ccact
     +			,uract,pact,eact,ract,qact

c-----------------------------------------------------------
c save common blocks
c-----------------------------------------------------------

        save /qflxi/, /qflxro/, /qflxrn/, /qflxra/

c-----------------------------------------------------------
c end of header
c-----------------------------------------------------------

