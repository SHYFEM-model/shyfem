
!--------------------------------------------------------------------------
!
! Copyright (C) 1999-2007 The GOTM group
!
! This file is part of SHYFEM. (m)
!
! This directory contains the code for GOTM, version 4.0.0
!
! Please see this web page for more information: http://gotm.net/
!
! GOTM is distributed under the GNU public license
!
!--------------------------------------------------------------------------

! !REVISION HISTORY:
!
!  03May99  ppm (grads)
!

      double precision kappa,z0b,
     &              cm0,z0s,ce1,ce2,ce3minus,ce3plus,
     &              k_min,epsmin,L_min,
     &              sig_k,sig_e,Prandtl0,cde,cdL,cmucst,
     &              sl,e1,e2,e3,a1,a2,b1,b2,c1,
     &              galp,qeghmin,qeghmax,qeghcrit,
     &              alfa,c2,c3,klimiw,rich_cr,numiw,nuhiw,numshear   
      integer MaxN,Stab,MYLength,TKEMeth,LengthMeth,Iwmodel
      logical fluxcond,lengthlim  
      logical qesmooth
      parameter(MaxN=250)

!     Different TKE Models
      integer TKEloceq,TKE_keps,TKEMY  
      parameter(TKEloceq=1)
      parameter(TKE_keps=2)
      parameter(TKEMY=3)

!     Different Lengthscale Models
      integer Parabola,Triangle,Xing,RobertOuellet,
     &        Blackadar,BougeaultAndre,Ispramix,
     &        DissEq,LengthEq,DistParabola  
      parameter(Parabola=1)
      parameter(Triangle=2)
      parameter(Xing=3)
      parameter(RobertOuellet=4)
      parameter(Blackadar=5)
      parameter(BougeaultAndre=6)
      parameter(Ispramix=7)
      parameter(DissEq=8)
      parameter(LengthEq=9)
      parameter(DistParabola=10)

! Different 'stability functions'
      integer KanClay, BurBaum, CanutoA, CanutoB, KanClayQe,
     &        BurBaumQe, CanutoAQe, CanutoBQe, Constan,
     &        MunkAnd, SchumGerz, FluxRich
      parameter(KanClay=1)
      parameter(BurBaum=2)
      parameter(CanutoA=3)
      parameter(CanutoB=4)
      parameter(KanClayQe=5)
      parameter(BurBaumQe=6)
      parameter(CanutoAQe=7)
      parameter(CanutoBQe=8)
      parameter(Constan=9)
      parameter(MunkAnd=10)
      parameter(SchumGerz=11)
      parameter(FluxRich=12)


      common /turbconstants/  kappa,z0b,
     &                    cm0,z0s,ce1,ce2,ce3minus,ce3plus,
     &                    k_min,epsmin,L_min,
     &                    sig_k,sig_e,Prandtl0,cde,cdL,
     &                    sl,e1,e2,e3,a1,a2,b1,b2,c1,
     &                    galp,cmucst,qeghmin,qeghmax,
     &                    qeghcrit,
     &                    alfa,c2,c3,klimiw,rich_cr,numiw,nuhiw,
     &                    numshear,Stab,TKEMeth,LengthMeth,IwModel,
     &                    MYLength,
     &                    fluxcond,lengthlim,  
     &                    qesmooth 

