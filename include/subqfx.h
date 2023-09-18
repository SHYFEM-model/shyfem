
!-----------------------------------------------------------
! integer variables
!-----------------------------------------------------------

!	ifunit			unit of flux file
!	itfold,itfnew,itfact	time for old/new/actual level
!	itperiod		periodic heat flux condition 
!				(0 if none, default)
!				(31536000 for one year - 365 days)
!	irhumid			relative humidity given 
!				(0 for wet bulb given)
!				(1 for relative humidity, default)

        integer ifunit,itfold,itfnew,itfact,itperiod,irhumid
        common /qflxi/ ifunit,itfold,itfnew,itfact,itperiod,irhumid

!-----------------------------------------------------------
! values of old time level
!-----------------------------------------------------------

        double precision qsold,taold,tbold,uwold,ccold,urold,pold,eold,rold,qold
        common /qflxro/ qsold,taold,tbold,uwold,ccold,urold,pold,eold,rold,qold

!-----------------------------------------------------------
! values of new time level
!-----------------------------------------------------------

        double precision qsnew,tanew,tbnew,uwnew,ccnew,urnew,pnew,enew,rnew,qnew
        common /qflxrn/ qsnew,tanew,tbnew,uwnew,ccnew,urnew,pnew,enew,rnew,qnew

!-----------------------------------------------------------
! values of actual time level
!-----------------------------------------------------------

        double precision qsact,taact,tbact,uwact,ccact,uract,pact,eact,ract,qact
        common /qflxra/ qsact,taact,tbact,uwact,ccact,uract,pact,eact,ract,qact

!-----------------------------------------------------------
! save common blocks
!-----------------------------------------------------------

        save /qflxi/, /qflxro/, /qflxrn/, /qflxra/

!-----------------------------------------------------------
! end of header
!-----------------------------------------------------------

