#!/bin/sh
#
#-------------------------------------------

ElabExt()
{
  numb=$1
  what=$2

  sim=mm_hyd_$numb.ext

  if [ -f $sim ]; then
    shyelab -split $sim
  fi

  TS $what
  gps -jpg -dpi 200 out.ps
  mv out.jpg ${what}_$numb.jpg
  mv out.ps ${what}_$numb.ps
}

ElabSalt()
{
  numb=$1
  what=$2
  apn=$3
  slide=$4

  [ -z "$slide" ] && slide=1
  sim=mm_hyd_$numb.ts.shy

  shyplot -freq $freq -varname salinity $sim $apn
  gps -split plot.ps
  gps -eps -jpg -dpi 200 plot.$slide.ps
  mv plot.$slide.jpg ${what}_$numb.jpg
  mv plot.$slide.eps ${what}_$numb.eps
}

ElabHydro()
{
  numb=$1
  what=$2
  apn=$3

  sim=mm_hyd_$numb.hydro.shy

  shyplot -freq $freq -varname velocity $sim $apn
  gps -split -eps -jpg -dpi 200 plot.ps
  mv plot.1.jpg ${what}_$numb.jpg
  mv plot.1.eps ${what}_$numb.eps
}

TS()
{
  var=$1

  tx="time"
  if [ $var = "zeta" ]; then
    title="Water level"
    ty="Water level [m]"
  elif [ $var = "zcomp" ]; then
    title="Water level"
    ty="Water level [m]"
  elif [ $var = "speed" ]; then
    title="Current speed"
    ty="Current speed [m/s]"
  elif [ $var = "salt" ]; then
    title="Salinity"
    ty="Salinity [psu]"
  elif [ $var = "scomp" ]; then
    title="Salinity"
    ty="Salinity [psu]"
  else
    echo "no match: $var"
    exit 1
  fi
  gp -t "$title" -tx "$tx" -ty "$ty" $var.2d.*
}

#-------------------------------------------

level=0
freq=32

if [ 1 = 2 ]; then

ElabExt 01 zeta
ElabExt 02 zeta
ElabExt 02 speed
ElabExt 03 zeta
ElabExt 03 speed

ElabHydro 03 velcol apnbath.str
ElabHydro 03 vel apn_vel.str

ElabExt 11 zeta
ElabExt 11 speed
ElabHydro 11 velcol apnbath.str
ElabHydro 11 vel apn_vel.str

ElabExt 11 zeta
cp zeta.2d.3 zcomp.2d.00
ElabExt 12a zeta
cp zeta.2d.3 zcomp.2d.30+
ElabExt 12b zeta
cp zeta.2d.3 zcomp.2d.30-
ElabExt 12 zcomp

ElabHydro 11 velriver apn_velriver.str
ElabHydro 13 velriver apn_velriver.str
ElabHydro 11 velcolriver apn_velcolriver.str
ElabHydro 13 velcolriver apn_velcolriver.str

ElabHydro 21 velcol apnbath.str
ElabHydro 21 vel apn_vel.str
ElabHydro 22a velcol apnbath.str
ElabHydro 22a vel apn_vel.str
ElabHydro 22b velcol apnbath.str
ElabHydro 22b vel apn_vel.str
ElabHydro 23 velcol apnbath.str
ElabHydro 23 vel apn_vel.str

freq=12
ElabSalt 31 saltmap apnbath.str 7

ElabExt 31 salt
cp salt.2d.3 scomp.2d.35
ElabExt 32a salt
cp salt.2d.3 scomp.2d.20
ElabExt 32b salt
cp salt.2d.3 scomp.2d.50
ElabExt 32 scomp
fi

ElabExt 41 salt
ElabExt 41 zeta
ElabExt 41 speed

freq=32
ElabHydro 41 velcol apnbath.str

#-------------------------------------------

