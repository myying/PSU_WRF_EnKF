#!/bin/bash
#This function calculates i/j_parent_start so that the nested domains are centerred over
#  the storm of interest. Storm location is from tcvitals data.
#Input:  icdate = start date of a wrf run
#Output: a text file containing i_parent_start and j_parent_start

. $CONFIG_FILE

icdate=$1
outfile=$2

#find tcvitals data that covers icdate
minute_off=`echo "(${icdate:8:2}*60+${icdate:10:2})%$LBC_INTERVAL" |bc`
fgdate=`advance_time $icdate -$minute_off`

for i in 0 1; do
  tcdate=`advance_time $fgdate $((i*$LBC_INTERVAL))`
  tcvitals_data=$TCVITALS_DIR/${tcdate:0:4}/${tcdate}.${STORM_ID}-tcvitals.dat
  if [ ! -f $tcvitals_data ]; then echo "$tcvitals_data not found!"; exit; fi
  latstr=`head -n1 $tcvitals_data |awk '{print $6}'`
  lonstr=`head -n1 $tcvitals_data |awk '{print $7}'`
  if [ ${latstr:3:1} == "N" ]; then
    lat[$i]=`echo "${latstr:0:3}/10" |bc -l`
  else
    lat[$i]=`echo "-${latstr:0:3}/10" |bc -l`
  fi
  if [ ${lonstr:4:1} == "E" ]; then
    lon[$i]=`echo "${lonstr:0:4}/10" |bc -l`
  else
    lon[$i]=`echo "-${lonstr:0:4}/10" |bc -l`
  fi
done

iclat=$(echo "${lat[0]}+(${lat[1]} - ${lat[0]})*$minute_off/$LBC_INTERVAL" |bc -l)
iclon=$(echo "${lon[0]}+(${lon[1]} - ${lon[0]})*$minute_off/$LBC_INTERVAL" |bc -l)

cat << EOF > ll2ij.ncl
begin
   opt = True
   opt@MAP_PROJ          = 3
   opt@TRUELAT1          = $TRUELAT1
   opt@TRUELAT2          = $TRUELAT2
   opt@STAND_LON         = $STAND_LON
   opt@REF_LAT           = $REF_LAT
   opt@REF_LON           = $REF_LON
   opt@KNOWNI            = `echo "${E_WE[0]}/2" |bc -l`
   opt@KNOWNJ            = `echo "${E_SN[0]}/2" |bc -l`
   opt@DX                = ${DX[0]}
   opt@DY                = ${DY[0]}
   loc=wrf_ll_to_ij($iclon,$iclat,opt)
   asciiwrite("tcij",round(loc,3))
end
EOF
ncl ll2ij.ncl
tci=`cat tcij |head -n1`
tcj=`cat tcij |tail -n1`
cat tcij
rm -f ll2ij.ncl tcij

istart=$(echo "$tci - ((${E_WE[1]}-1)/${GRID_RATIO[1]}-1)/2 + 1" |bc)
jstart=$(echo "$tcj - ((${E_SN[1]}-1)/${GRID_RATIO[1]}-1)/2 + 1" |bc)

i_parent_start="1 $istart"
j_parent_start="1 $jstart"
if [ $MAX_DOM -gt 2 ]; then
  for n in `seq 3 $MAX_DOM`; do
    i_parent_start="$i_parent_start ${I_PARENT_START[$n-1]}"
    j_parent_start="$j_parent_start ${J_PARENT_START[$n-1]}"
  done
fi

echo $i_parent_start > $outfile
echo $j_parent_start >> $outfile

