#!/bin/bash

. $CONFIG_FILE

rundir=$WORK_DIR/run/$DATE/wrf_window1
if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi

cd $rundir
if [[ `cat stat` == "complete" ]]; then exit; fi

#Check dependency
wait_for_module ../icbc 
if $RUN_ENKF; then wait_for_module ../../$PREVDATE/wrf_ens; fi

echo "  Running WRF through obs window (fg -> fg02)..."
echo running > stat

indate=`advance_time $DATE $OBS_WIN_MIN`
outdate=`advance_time $DATE $OBS_WIN_MAX`

#1. Update BC with new wrfinput
for i in 1; do
  touch update_wrf_bc.log
  if [[ `tail -n5 update_wrf_bc.log |grep successfully` ]]; then continue; fi

  dd=`diff_time $DATE_START $LBDATE`
  n_1=$((dd/$LBC_INTERVAL+1))
 
  cat > parame.in << EOF
&control_param
 wrf_3dvar_output_file = 'wrfinput_d01_update'
 wrf_bdy_file          = 'wrfbdy_d01_update'
 wrf_bdy_file_real     = 'wrfbdy_d01_real'
 wrf_input_from_si     = 'wrfinput_d01_real'
 wrf_input_from_si_randmean = 'wrfinput_d01_real'
 wrf_3dvar_random_draw = 'wrfinput_d01_real'
 cycling = .true.
 debug   = .true. 
 low_bdy_only = .false. 
 perturb_bdy = .false.
 n_1 = $n_1
/
EOF
  ln -fs $WRF_BC_DIR/update_wrf_bc.exe .
  ln -fs ../../../fc/wrfbdy_d01_window wrfbdy_d01_real
  cp -L wrfbdy_d01_real wrfbdy_d01_update
  ln -fs ../../../rc/$DATE/wrfinput_d01 wrfinput_d01_real
  rm -f wrfinput_d01_update

  if $RUN_ENKF; then
    ln -fs ../../../fc/$PREVDATE/wrfinput_d01_`wrf_time_string $indate`_mean wrfinput_d01_update
  else
    ln -fs ../../../fc/$PREVDATE/wrfinput_d01_`wrf_time_string $indate` wrfinput_d01_update
  fi

  ./update_wrf_bc.exe >& update_wrf_bc.log 
  watch_log update_wrf_bc.log successfully 1 $rundir
  mv wrfbdy_d01_update $WORK_DIR/fc/wrfbdy_d01_window
done

#2. Run wrf.exe
for i in 1; do
  touch rsl.error.0000
  if [[ `tail -n5 rsl.error.0000 |grep SUCCESS` ]]; then continue; fi

  ln -fs $WRF_DIR/run/* .
  rm -f namelist.*

  for n in `seq 1 $MAX_DOM`; do
    dm=d`expr $n + 100 |cut -c2-`
    if $RUN_ENKF; then
      ln -fs ../../../fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $indate`_mean wrfinput_${dm}
    else
      ln -fs ../../../fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $indate` wrfinput_${dm}
    fi
  done
  ln -fs ../../../fc/wrfbdy_d01_window wrfbdy_d01
  if [[ $SST_UPDATE == 1 ]]; then
    ln -fs ../../../rc/$DATE_START/wrflowinp_d?? .
  fi

  export start_date=$indate
  export run_minutes=$(echo "$OBS_WIN_MAX - $OBS_WIN_MIN" |bc)
  export inputout_interval=$run_minutes
  export inputout_begin=0
  export inputout_end=$run_minutes

  if $MULTI_PHYS_ENS; then
    $SCRIPT_DIR/multi_physics_reset.sh >& multi_physics_reset.log
  fi

  if $FOLLOW_STORM; then
    cp $WORK_DIR/rc/$DATE/ij_parent_start_4dvar ij_parent_start
  fi
  $SCRIPT_DIR/namelist_wrf.sh wrf > namelist.input
  $SCRIPT_DIR/job_submit.sh $wrf_single_ntasks 0 $HOSTPPN ./wrf.exe >& wrf.log
done

#Check output
watch_log rsl.error.0000 SUCCESS 1 $rundir

##wrfinput for next cycle
outfile=wrfinput_d01_`wrf_time_string $outdate`
watch_file $outfile 1 $rundir
if $RUN_ENKF; then
  mv $outfile $WORK_DIR/fc/$PREVDATE/wrfinput_d01_`wrf_time_string $outdate`_mean
else
  mv $outfile $WORK_DIR/fc/$PREVDATE/wrfinput_d01_`wrf_time_string $outdate`
fi
if [ $MAX_DOM -gt 1 ]; then
  for n in `seq 2 $MAX_DOM`; do
    dm=d`expr $n + 100 |cut -c2-`
    outfile=wrfout_${dm}_`wrf_time_string $outdate`
    watch_file $outfile 1 $rundir
    if $RUN_ENKF; then
      mv $outfile $WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $outdate`_mean
    else
      mv $outfile $WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $outdate`
    fi
  done
fi

if $CLEAN; then rm -f wrfout*; fi

echo complete > stat

