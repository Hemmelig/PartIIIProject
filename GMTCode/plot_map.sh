#!/bin/csh
 
# gmt settings 
 gmtset COLOR_MODEL rgb
 gmtset LABEL_FONT_SIZE 18 LABEL_OFFSET 0.2c ANNOT_FONT Helvetica ANNOT_FONT_SIZE 10 ANNOT_OFFSET 0.05c HEADER_FONT_SIZE 14 HEADER_FONT Helvetica-Bold  HEADER_OFFSET 0.0c
 gmtset BASEMAP_TYPE fancy FRAME_PEN 0.75p TICK_LENGTH 0.15c
 
 gmtset OUTPUT_DEGREE_FORMAT +D
 
fig=map_events # figure name: $fig.pdf

[ -e $fig.ps ] && rm $fig.ps
 
topo=./etopo2.grd
#cpt=../Resource/GMT_globe.cpt
cpt=topo.cpt
makecpt -Crelief -T-8000/5000/100 -D -Z > $cpt
#
cpt_eq=quakes.cpt  
makecpt -Cseis -T400/650/50 -D -Z > $cpt_eq # make color bar for source depth
 
scale=7i     # size: width
 
lon1=150
lon2=200
lat1=-40
lat2=-10
 
Rvalue=$lon1/$lon2/$lat1/$lat2
 
grd=loc.grd
norm=loc.norm
grdsample -G$grd $topo -I2m -R$Rvalue -V -fg
grdgradient $grd -G$norm -A0 -M -Nt1 -V
 
grdimage -R$Rvalue -JM$scale -Ba10f10WSne:."Jan 2016 - Apr 2017": -X0.75i -Y2i -C$cpt $grd -I$norm -K  -N1 > $fig.ps
#pscoast -R -J -O -K  -W0.2p/100 -Di  >> $fig.ps
 
inp=./datagcmt.txt
awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' $inp | psmeca -R  -J -Sm0.15i/-1 -O -K -W0.2p/0 -T0 -Z$cpt_eq  >> $fig.ps
 
psscale -C$cpt_eq -D7.5i/2i/-3.5i/0.5 -Ba50g50:"Depth (km)": -O  -K >> $fig.ps
pwd | psxy -R -J -O -H >> $fig.ps

#psconvert $fig.ps -A -E1000 -P -Tg
ps2pdf -A $fig.ps
evince $fig.pdf 
