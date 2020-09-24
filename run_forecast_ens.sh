#!/bin/bash
source ~/.bashrc

#load configuration files, functions, parameters
cd $WORK/PSU_WRF_EnKF
export CONFIG_FILE=$WORK/PSU_WRF_EnKF/config/Patricia/EnSRF_s1_fcst
. $CONFIG_FILE
. util.sh

export DATE=201510212100
export DATE_FORECAST_END=201510230000

rundir=$WORK_DIR/run/$DATE/wrf_ens_fcst
if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi

cd $rundir
if [[ `cat stat` == "complete" ]]; then exit; fi

#Check dependency
wait_for_module ../icbc
if [ $DATE -gt $DATE_START ]; then
  if $RUN_ENKF; then wait_for_module ../enkf; fi
fi

echo running > stat

tid=0
nt=$((total_ntasks/$wrf_ntasks))
for NE in `seq 1 $NUM_ENS`; do
  id=`expr $NE + 1000 |cut -c2-`
  if [[ ! -d $id ]]; then mkdir $id; fi
  touch $id/rsl.error.0000
  if [[ `tail -n5 $id/rsl.error.0000 |grep SUCCESS` ]]; then continue; fi

  cd $id

  ln -fs $WORK/code/WRFV3.9_vortexfollow/run/* .
  rm -f namelist.*

  for n in `seq 1 2`; do
    dm=d`expr $n + 100 |cut -c2-`
    ln -fs ../../../../fc/$DATE/wrfinput_${dm}_$id wrfinput_$dm
  done
  ln -fs ../../../../fc/wrfbdy_d01 wrfbdy_d01

  ###nestdown for d03
  export run_minutes=0
  export start_date=$DATE
  i_parent_start="1 "
  j_parent_start="1 "
  outfile=$WORK_DIR/fc/$DATE/wrfinput_d02_$id
  i_parent_start="$i_parent_start $(ncdump -h $outfile |grep :I_PARENT_START |awk '{print $3}')"
  j_parent_start="$j_parent_start $(ncdump -h $outfile |grep :J_PARENT_START |awk '{print $3}')"
  i_parent_start="$i_parent_start ${I_PARENT_START[2]}"
  j_parent_start="$j_parent_start ${J_PARENT_START[2]}"
  echo $i_parent_start > ij_parent_start
  echo $j_parent_start >> ij_parent_start
  for n in 3; do
    dm=d`expr $n + 100 |cut -c2-`
    parent_dm=d`expr ${PARENT_ID[$n-1]} + 100 |cut -c2-`
    if [[ ! -d $dm ]]; then mkdir -p $dm; fi
    touch $dm/rsl.error.0000
    if [[ `tail -n5 $dm/rsl.error.0000 |grep SUCCESS` ]]; then continue; fi
    cd $dm
    $SCRIPT_DIR/namelist_wrf.sh ndown $n > namelist.input
    ln -fs $WRF_DIR/run/ndown.exe .
    ln -fs ../wrfinput_${parent_dm} wrfout_d01_`wrf_time_string $DATE`
    ln -fs ../../../../../rc/$DATE/wrfinput_${dm} wrfndi_d02
    cp -L wrfndi_d02 wrfinput_d02
    $SCRIPT_DIR/job_submit.sh $wps_ntasks 0 $HOSTPPN ./ndown.exe >& ndown.log &
    watch_log rsl.error.0000 SUCCESS 1 `pwd`
    cd ..
  done
  ln -fs d03/wrfinput_d02 wrfinput_d03

  ####Running model
  export start_date=$DATE
  export run_minutes=`diff_time $DATE $DATE_FORECAST_END`
  export restart_interval=1440
  $SCRIPT_DIR/namelist_wrf.sh wrf_fcst > namelist.input
  for n in `seq 1 $MAX_DOM`; do
    dm=d`expr $n + 100 |cut -c2-`
    echo "+:h:0:H_DIABATIC" > my_output_${dm}.txt
  done

  $SCRIPT_DIR/job_submit.sh $wrf_ntasks $((tid*$wrf_ntasks)) $HOSTPPN ./wrf.exe '08:00:00' >& wrf.log &

  tid=$((tid+1))
  if [[ $tid == $nt ]]; then
    tid=0
    wait
  fi
  cd ..
done
