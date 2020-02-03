#!/bin/bash
source ~/.bashrc

#load configuration files, functions, parameters
cd $WORK/PSU_WRF_EnKF
export CONFIG_FILE=$WORK/PSU_WRF_EnKF/config/Patricia/control_cheyenne
. $CONFIG_FILE
. util.sh

export DATE=201510211200

rundir=$WORK_DIR/run/$DATE/wrf
if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi

cd $rundir
if [[ `cat stat` == "complete" ]]; then exit; fi

#Check dependency
wait_for_module ../icbc
if [ $DATE -gt $DATE_START ]; then
  if $RUN_ENKF; then wait_for_module ../enkf; fi
fi

echo running > stat

start_date=$DATE
run_minutes=$LBC_INTERVAL


####FINISH RUNNING FIRST 6h
###Prepare icbc
#ln -fs $WORK/code/WRFV3.9_vortexfollow/run/* .
#rm -f namelist.*

#for n in `seq 1 $MAX_DOM`; do
#  dm=d`expr $n + 100 |cut -c2-`
#  if $RUN_ENKF && [ $n == 2 ]; then
#    ncea ../../../fc/$DATE/wrfinput_${dm}_??? wrfinput_$dm
#  else
#    ncea ../../../fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE`_??? wrfinput_$dm
#  fi
#done
#cp ../../../rc/$DATE_START/wrfinput_d03 .
#ln -fs ../../../fc/wrfbdy_d01 wrfbdy_d01

###Nestdown inner domains
#for n in 3; do
#  dm=d`expr $n + 100 |cut -c2-`
#  parent_dm=d`expr ${PARENT_ID[$n-1]} + 100 |cut -c2-`
#  if [[ `tail -n5 $dm/rsl.error.0000 |grep SUCCESS` ]]; then continue; fi
#  if [[ ! -d $dm ]]; then mkdir -p $dm; fi
#  cd $dm
#  ln -fs ../wrfinput_d0? .
#  export run_minutes=0
#  export start_date=$DATE
#  if $FOLLOW_STORM; then
#    cp $WORK_DIR/rc/$DATE/ij_parent_start .
#  fi
#  $SCRIPT_DIR/namelist_wrf.sh ndown $n > namelist.input
#  rm -f wrfinput_d0?
#  ln -fs $WRF_DIR/run/ndown.exe .
#  ln -fs ../wrfinput_${parent_dm} wrfout_d01_`wrf_time_string $DATE`
#  ln -fs ../wrfinput_$dm wrfndi_d02
#  $SCRIPT_DIR/job_submit.sh $wps_ntasks 0 $HOSTPPN ./ndown.exe >& ndown.log
#  cd ..
#  watch_log $dm/rsl.error.0000 SUCCESS 1 $rundir
#  mv $dm/wrfinput_d02 wrfinput_${dm}
#done
next_date=`advance_time $start_date $LBC_INTERVAL`
start_date=$next_date

while [[ $next_date -le $DATE_END ]]; do  #time loop
  next_date=`advance_time $start_date $LBC_INTERVAL`
  echo $start_date

  export start_date
  export run_minutes
  export time_step_ratio=1
  export sst_update=0

  touch rsl.error.0000_$start_date
  if [[ `tail -n5 rsl.error.0000_$start_date |grep SUCCESS` ]]; then
    start_date=$next_date
    continue
  fi
  for n in `seq 1 $MAX_DOM`; do
    dm=d`expr $n + 100 |cut -c2-`
  done

  ####prepare sst
  mkdir -p sst_icbc
  cd sst_icbc

  ###get child domain location
  i_parent_start="1 "
  j_parent_start="1 "
  for n in `seq 2 $MAX_DOM`; do
    dm=d`expr $n + 100 |cut -c2-`
    outfile=../wrfrst_${dm}_`wrf_time_string $start_date`
    i_parent_start="$i_parent_start $(ncdump -h $outfile |grep :I_PARENT_START |awk '{print $3}')"
    j_parent_start="$j_parent_start $(ncdump -h $outfile |grep :J_PARENT_START |awk '{print $3}')"
  done
  echo $i_parent_start > ij_parent_start
  echo $j_parent_start >> ij_parent_start

  $SCRIPT_DIR/namelist_wps.sh > namelist.wps

  ln -fs $WORK/data/Patricia/icbc/sst/FILE* .

  ln -fs $WPS_DIR/geogrid/src/geogrid.exe .
  $SCRIPT_DIR/job_submit.sh $wps_ntasks 0 $HOSTPPN ./geogrid.exe >& geogrid.log
  watch_log geogrid.log Successful 10 `pwd`

  ln -fs $WPS_DIR/metgrid/METGRID.TBL.ARW METGRID.TBL
  ln -fs $WPS_DIR/metgrid/src/metgrid.exe .
  $SCRIPT_DIR/job_submit.sh $wps_ntasks 0 $HOSTPPN ./metgrid.exe >& metgrid.log
  watch_log metgrid.log Successful 10 `pwd`

  $SCRIPT_DIR/namelist_wrf.sh real > namelist.input
  ln -fs $WORK/code/WRFV3.9_serial/main/real.exe .
  ./real.exe >& rsl.error.0000
  watch_log rsl.error.0000 SUCCESS 10 `pwd`
  for n in `seq 1 $MAX_DOM`; do
    dm=d`expr $n + 100 |cut -c2-`
    mv wrfinput_$dm wrfinput_${dm}_$start_date
  done

  cd ..
  cp sst_icbc/ij_parent_start .
  export restart=true
  $SCRIPT_DIR/namelist_wrf.sh wrf_fcst > namelist.input

  for n in `seq 1 $MAX_DOM`; do
    dm=d`expr $n + 100 |cut -c2-`
    $SCRIPT_DIR/diagnose/update_var.py sst_icbc/wrfinput_${dm}_$start_date wrfrst_${dm}_`wrf_time_string $start_date` SST
    $SCRIPT_DIR/diagnose/update_var.py sst_icbc/wrfinput_${dm}_$start_date wrfrst_${dm}_`wrf_time_string $start_date` TSK
  done

  $SCRIPT_DIR/job_submit.sh $wrf_ntasks 0 $HOSTPPN ./wrf.exe >& wrf.log
  watch_log rsl.error.0000 SUCCESS 10 $WORK_DIR/run/$DATE/wrf
  mv rsl.error.0000 rsl.error.0000_$start_date

  start_date=$next_date

done


