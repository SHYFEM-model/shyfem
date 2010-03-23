
        real  ! param fito1
     :  mumax1,tmax1,tott1,coeff1,iopt1,
     :  ksp1,krf1,kmf1

        real  ! param fito2
     :   mumax2,tmax2,tott2,coeff2,
     :  ksdoc2,ksp2,kmf2

        real  ! parametri zoo
     :   kgr1,kfz1,alfa,eff1,eff2,kmz,kexz

        real  ! param other
     :  kdcd,kdcp,ksodec,used,krear,kestf,kestd,kestw,sal,
     :  rpc1,rpc2,rpcz,roc1,roc2,rocd

        common /param/
     :  mumax1,tmax1,tott1,coeff1,iopt1,
     :  ksp1,krf1,kmf1,
     :  mumax2,tmax2,tott2,coeff2,
     :  ksdoc2,ksp2,kmf2,
     :  kgr1,kfz1,alfa,eff1,eff2,kmz,kexz,
     :  kdcd,kdcp,ksodec,used,krear,kestf,kestd,kestw,sal,
     :  rpc1,rpc2,rpcz,roc1,roc2,rocd

         integer po4,phy1,bac,zoo,oxy,detc,detp,doc,dop

         common /variabili/
     :   po4,phy1,bac,zoo,oxy,detc,
     :   detp,doc,dop

         integer nvd,nvnd,nvt

         parameter(nvd=9)
         parameter(nvnd=1)
         parameter(nvt=nvd+nvnd)

