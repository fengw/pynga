#! /bin/csh

set MAGS = ( 6.50 7.00)	# Magnitude
set FLEN = ( 18.0 44.0)  # fault length
set FWID = ( 18.0 23.0)  # fault width

set NSTK = ( 3 5 )	# number of stations along strike
set DSTK = ( 6 8 ) # increment of stations along strike
set NDIP = ( 2 5 ) # number of stations at edge of strike, along dip
set DDIP = ( 8 8 ) # increment of stations along dip (starting with -1*ddip, up to ndip*ddip)

set ELON = -122.0
set ELAT = 36.5
set XAZIM = 0  # fault azimuth

set GENERIC_DIST = 1.0
set GENERIC_VS30 = 865

set NN_PROF_INC = ( 4 5)	# increment in North direction for stations beyond edge of fault
set NN_STAT = 4	# number of stations beyond edge of fault (along strike)
set EE_PROF_DISTS = ( -100 -90 -80 -70 -60 -50 -40 -30 -25 -20 -15 -12 -8 -5 -3 -1 0 \
1 3 5 8 12 15 20 25 30 40 50 60 70 80 90 100) # station sequence for the 'along strike' stations

set m = 0
foreach mag ( $MAGS )
@ m ++

set TC_STATS = rv01-m${mag}_stats.ne
set LL_STATS = rv01-m${mag}_stats.ll
set STATLIST_BBP = rv01-m${mag}_stats.stl

echo $EE_PROF_DISTS | gawk -v fl=$FLEN[$m] -v ninc=$NN_PROF_INC[$m] -v nsta=$NN_STAT \
-v ns=$NSTK[$m] -v ds=$DSTK[$m] -v nd=$NDIP[$m] -v dd=$DDIP[$m] 'BEGIN{np=0;}{ \
for(i=1;i<=NF;i++){np++;pd[i]=$i;}} \
END { \
for(j=0;j<ns;j++){ \
xx=(int(ns/2)-j)*ds; \
if(xx>=0.0)xnam=sprintf("p%.3d",xx); \
else xnam=sprintf("m%.3d",-xx); \
for(i=1;i<=np;i++){ \
if(pd[i]>=0.0)ynam=sprintf("p%.3d",pd[i]); \
else ynam=sprintf("m%.3d",-pd[i]); \
printf "%8.2f %8.2f %s%s\n",xx,pd[i],xnam,ynam;}} \
fl2=0.5*fl; \
for(j=3;j<=(2*nd+3);j++){ \
yy=-(nd-(j-1))*dd; \
if(yy>=0.0)ynam=sprintf("p%.3d",yy); \
else ynam=sprintf("m%.3d",-yy); \
for(i=0;i<=nsta;i++){ \
xx=fl2+i*ninc; \
if(xx>=0.0)xnam=sprintf("p%.3d",xx); \
else xnam=sprintf("m%.3d",-xx); \
printf "%8.2f %8.2f %s%s\n",xx,yy,xnam,ynam;} \
for(i=0;i<=nsta;i++){ \
xx=-(fl2+i*ninc); \
if(xx>=0.0)xnam=sprintf("p%.3d",xx); \
else xnam=sprintf("m%.3d",-xx); \
printf "%8.2f %8.2f %s%s\n",xx,yy,xnam,ynam;} \
}}' > $TC_STATS

./xy2ll mlon=$ELON mlat=$ELAT xazim=$XAZIM < $TC_STATS > $LL_STATS

gawk -v dx=$GENERIC_DIST -v vs=$GENERIC_VS30 -v title="Mag= $mag" \
'BEGIN{ \
printf "# %s\n",title; \
printf "# Date: %s\n",strftime(); \
printf "# \n"; \
printf "# Data fields are TAB-separated\n"; \
printf "# \n"; \
printf "# Column 1: station longitude\n"; \
printf "# Column 2: station latitude\n"; \
printf "# Column 3: station name/code\n"; \
printf "# Column 4: dummy distance to fault plane (set to 1 km)\n"; \
printf "# Column 5: Vs30 (m/s) (set to 865 m/s)\n"; \
printf "#\n";}{ \
printf "%.5f\t%.5f\t%s\t%.2f\t%.0f\n",$1,$2,$3,dx,vs;}' $LL_STATS > $STATLIST_BBP

end
