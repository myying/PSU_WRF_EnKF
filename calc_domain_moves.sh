#!/bin/bash
#This function generates domain_moves for namelist.input so that the nested domains follows 
#  storm center (from tcvitals data). Storm center lat/lon are calculated at each time (including
#  start/end time for wrf and any lbc (0 6 12 18Z) time in between. The storm is assumed to move
#  linearly between any two given location points.
#      
#  Only move the smallest domain, and use corral_dist to force move its parent domains
#
#Input: date1 = wrf run start date
#       date2 = wrf run end date
#Output: a text file that contains namelist options related to domain moves.

. $CONFIG_FILE

date1=$1
date2=$2
outfile=$3

# Find tcvitals data that covers the period
off1=`echo "(${date1:8:2}*60+${date1:10:2})%$LBC_INTERVAL" |bc`
off2=`echo "(${date2:8:2}*60+${date2:10:2})%$LBC_INTERVAL" |bc`
fgdate1=`advance_time $date1 -$off1`
fgdate2=`advance_time $date2 $((LBC_INTERVAL-$off2))`
n=$(echo "`diff_time $fgdate1 $fgdate2`/$LBC_INTERVAL" |bc)

# Get storm lat lon at each time node
for i in `seq 0 $n`; do
  tcdate=`advance_time $fgdate1 $((i*$LBC_INTERVAL))`
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

# Interpolate if time node is not at lbc time
iclat1=$(echo "${lat[0]} + (${lat[1]} - ${lat[0]})*$off1/$LBC_INTERVAL" |bc -l)
iclon1=$(echo "${lon[0]} + (${lon[1]} - ${lon[0]})*$off1/$LBC_INTERVAL" |bc -l)
iclat2=$(echo "${lat[$n-1]} + (${lat[$n]} - ${lat[$n-1]})*$off2/$LBC_INTERVAL" |bc -l)
iclon2=$(echo "${lon[$n-1]} + (${lon[$n]} - ${lon[$n-1]})*$off2/$LBC_INTERVAL" |bc -l)
if [ $off2 -gt 0 ]; then
  n2=$((n-1))
else
  n2=$((n-2))
fi
latlist=""
lonlist=""
datelist=""
m=0
for i in `seq 1 $n2`; do 
  latlist="$latlist ${lat[$i]}"
  lonlist="$lonlist ${lon[$i]}"
  datelist="$datelist `advance_time $fgdate1 $((i*$LBC_INTERVAL))`"
  m=$((m+1))
done
latlist="$iclat1 $latlist $iclat2"
lonlist="$iclon1 $lonlist $iclon2"
datelist="$date1 $datelist $date2"
m=$((m+1))
nlat=($latlist)
nlon=($lonlist)
ndate=($datelist)

# Get domain d01 i/j at storm center each time
for i in `seq 0 $m`; do
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
   loc=wrf_ll_to_ij(${nlon[$i]},${nlat[$i]},opt)
   asciiwrite("ij",round(loc,3))
end
EOF
  ncl ll2ij.ncl
  ni[$i]=`cat ij |head -n1`
  nj[$i]=`cat ij |tail -n1`
  rm -f ll2ij.ncl ij
done


# Generate domain_moves options
num_moves=0
dt=0
move_id=""
move_interval=""
move_cd_x=""
move_cd_y=""
corral_dist="8, 8,"
for i in `seq 3 $MAX_DOM`; do
  corral_dist="$corral_dist $(echo "`min ${I_PARENT_START[$i-1]} ${J_PARENT_START[$i-1]}`-2" |bc),"
done

for i in `seq 0 $((m-1))`; do  #For each time
  echo "${ndate[$i]} (${ni[$i]} ${nj[$i]}) -> ${ndate[$i+1]} (${ni[$i+1]} ${nj[$i+1]})"
  dmin=`diff_time ${ndate[$i]} ${ndate[$i+1]}`
  di=`echo "(${ni[$i+1]} - ${ni[$i]})" |bc`
  dj=`echo "(${nj[$i+1]} - ${nj[$i]})" |bc`
  s=$(echo "`abs $di` + `abs $dj`" |bc)
  echo $dmin $s $di $dj
  if [ $s -gt 0 ]; then  # If displacement is non-zero, move the domain
    num_moves=$((num_moves+$s))
    for t in `seq 1 $s`; do
      if [ $t -lt $s ]; then
        dt=`echo "$dt + $dmin/$s" |bc`
      else
        dt=`echo "$dt + $dmin - ($s-1)*$dmin/$s" |bc`
      fi
      move_id="$move_id $MAX_DOM,"
      move_interval="$move_interval $dt,"
      if [ `abs $di` -gt 0 ] && [ `abs $dj` -gt 0 ]; then
        if [ `abs $di` -ge `abs $dj` ]; then
          ds=$(echo "$s/`abs $dj`" |bc)
          if [ `echo "$t % $ds" |bc` == 0 ] && [ $t -le $(echo "`abs $dj` * $ds" |bc) ]; then
            move_cd_x="$move_cd_x 0,"
            move_cd_y="$move_cd_y $(echo "$dj/`abs $dj`" |bc),"
          else
            move_cd_x="$move_cd_x $(echo "$di/`abs $di`" |bc),"
            move_cd_y="$move_cd_y 0,"
          fi
        else
          ds=$(echo "$s/`abs $di`" |bc)
          if [ `echo "$t % $ds" |bc` == 0 ] && [ $t -le $(echo "`abs $di` * $ds" |bc) ]; then
            move_cd_x="$move_cd_x $(echo "$di/`abs $di`" |bc),"
            move_cd_y="$move_cd_y 0,"
          else
            move_cd_x="$move_cd_x 0,"
            move_cd_y="$move_cd_y $(echo "$dj/`abs $dj`" |bc),"
          fi
        fi
      elif [ `abs $di` -gt 0 ]; then
        move_cd_x="$move_cd_x $(echo "$di/`abs $di`" |bc),"
        move_cd_y="$move_cd_y 0,"
      else
        move_cd_x="$move_cd_x 0,"
        move_cd_y="$move_cd_y $(echo "$dj/`abs $dj`" |bc),"
      fi
    done
  fi
done

# Output
echo "num_moves     = $num_moves," > $outfile
if [ $num_moves -gt 0 ]; then
  echo "move_id       = $move_id" >> $outfile
  echo "move_interval = $move_interval" >> $outfile
  echo "move_cd_x     = $move_cd_x" >> $outfile
  echo "move_cd_y     = $move_cd_y" >> $outfile
fi
echo "corral_dist   = $corral_dist" >> $outfile

