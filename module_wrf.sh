#!/bin/bash
. $CONFIG_FILE

rundir=$WORK_DIR/run/$DATE/wrf
if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi

cd $rundir
if [[ `cat stat` == "complete" ]]; then exit; fi

#Check dependency
wait_for_module ../icbc
if $RUN_ENKF; then  
  wait_for_module ../enkf
  if [ $JOB_SUBMIT_MODE == 1 ]; then  wait_for_module ../wrf_ens; fi
fi
if $RUN_4DVAR; then  wait_for_module ../4dvar ../wrf_window; fi

echo running > stat

export start_date=$DATE
export run_minutes=$run_minutes_forecast


for i in 1; do
  touch rsl.error.0000
  if [[ `tail -n5 rsl.error.0000 |grep SUCCESS` ]]; then continue; fi

  ln -fs $WRF_DIR/run/* .
  rm -f namelist.*

#  ln -fs ../../../fc/$DATE/wrfinput_d01_mean wrfinput_d01  ####
  for n in `seq 1 $MAX_DOM`; do
    dm=d`expr $n + 100 |cut -c2-`
    if $RUN_4DVAR; then
      ln -fs ../../../fc/$DATE/wrfinput_${dm}_`wrf_time_string $DATE` wrfinput_$dm
    else
      if $RUN_ENKF; then
        ln -fs ../../../fc/$DATE/wrfinput_${dm}_mean wrfinput_$dm
      fi
    fi
  done
  #ln -fs ../../../rc/$DATE_START/wrfbdy_d01 .
  #ln -fs ../../../fc/$DATE/wrfbdy_d01 .
  ln -fs ../../../fc/wrfbdy_d01 .
  if [[ $SST_UPDATE == 1 ]]; then
    ln -fs ../../../rc/$DATE_START/wrflowinp_d?? .
  fi

  if $FOLLOW_STORM; then
    cp $WORK_DIR/rc/$DATE/ij_parent_start_4dvar ij_parent_start
    cp $WORK_DIR/rc/$DATE/domain_moves_4dvar domain_moves
  fi
  if $MULTI_PHYS_ENS; then
    $SCRIPT_DIR/multi_physics_reset.sh >& multi_physics_reset.log
  fi
  $SCRIPT_DIR/namelist_wrf.sh wrf > namelist.input

  $SCRIPT_DIR/job_submit.sh $wrf_single_ntasks 0 $HOSTPPN ./wrf.exe >& wrf.log
done

#Check output
watch_log rsl.error.0000 SUCCESS 1 $rundir

mv wrfout* $WORK_DIR/output/$DATE

echo complete > stat

