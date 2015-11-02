#!/bin/bash
. $CONFIG_FILE

rundir=$WORK_DIR/run/$DATE/wrf_window
if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi

cd $rundir
if [[ `cat stat` == "complete" ]]; then exit; fi

#Check dependency
wait_for_module ../icbc 
if [ $DATE -ge $DATE_CYCLE_START ]; then wait_for_module ../4dvar; fi

echo running > stat

if [ $DATE == $DATE_START ]; then
  echo "  Running initial WRF deterministic forecast to first cycle win_min..."
  export start_date=$DATE
  export run_minutes=$(echo "$run_minutes_cycle+$OBS_WIN_MIN" |bc)
else
  if $RUN_ENKF; then
    echo "  Running WRF deterministic forecast (4DVar analysis -> correct analysis time)..."
    export start_date=`advance_time $DATE $OBS_WIN_MIN`
    export run_minutes=`abs $OBS_WIN_MIN`
  else
    echo "  Running WRF through obs window again (4DVar analysis -> next cycle win_min)..."
    export start_date=`advance_time $DATE $OBS_WIN_MIN`
    export run_minutes=$run_minutes_cycle
  fi
fi
export inputout_interval=$run_minutes
export inputout_begin=0
export inputout_end=$run_minutes

#1. Update BC with 4dvar analysis if DA is performed
if [ $DATE -ge $DATE_CYCLE_START ]; then

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
  ln -fs ../../../fc/wrfbdy_d01 wrfbdy_d01_real
  cp -L wrfbdy_d01_real wrfbdy_d01_update
  ln -fs ../../../rc/$DATE/wrfinput_d01 wrfinput_d01_real
  rm -f wrfinput_d01_update

  ln -fs ../../../fc/$DATE/wrfinput_d01_`wrf_time_string $start_date` wrfinput_d01_update

  ./update_wrf_bc.exe >& update_wrf_bc.log 

  watch_log update_wrf_bc.log successfully 1 $rundir

  mv wrfbdy_d01_update $WORK_DIR/fc/wrfbdy_d01
done
fi

#2. Run wrf.exe
for i in 1; do
  touch rsl.error.0000
  if [[ `tail -n5 rsl.error.0000 |grep SUCCESS` ]]; then continue; fi

  ln -fs $WRF_DIR/run/* .
  rm -f namelist.*

  for n in `seq 1 $MAX_DOM`; do
    dm=d`expr $n + 100 |cut -c2-`
    if [ $DATE == $DATE_START ]; then
      ln -fs ../../../rc/$DATE/wrfinput_$dm .
    else
      ln -fs ../../../fc/$DATE/wrfinput_${dm}_`wrf_time_string $start_date` wrfinput_$dm
    fi
  done
  if [ $DATE == $DATE_START ]; then
    ln -fs ../../../rc/$DATE/wrfbdy_d01 .
  else
    ln -fs ../../../fc/wrfbdy_d01 .
  fi
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

if $RUN_ENKF; then
##wrfinput at $DATE for enkf replace mean (hybrid)
  outfile=wrfinput_d01_`wrf_time_string $DATE`
  watch_file $outfile 1 $rundir
  mv $outfile $WORK_DIR/fc/$DATE/wrfinput_d01_`wrf_time_string $DATE`
  if [ $MAX_DOM -gt 1 ]; then
    for n in `seq 2 $MAX_DOM`; do
      dm=d`expr $n + 100 |cut -c2-`
      outfile=wrfout_${dm}_`wrf_time_string $DATE`
      watch_file $outfile 1 $rundir
      cp $outfile $WORK_DIR/fc/$DATE/wrfinput_${dm}_`wrf_time_string $DATE`
    done
  fi
else
##save wrfinput as next cycle fg (4dvar only)
  for i in $OBS_WIN_MIN; do
    outdate=`advance_time $NEXTDATE $i`
    outfile=wrfinput_d01_`wrf_time_string $outdate`
    watch_file $outfile 1 $rundir
    mv $outfile $WORK_DIR/fc/$DATE/wrfinput_d01_`wrf_time_string $outdate`
    if [ $MAX_DOM -gt 1 ]; then
      for n in `seq 2 $MAX_DOM`; do
        dm=d`expr $n + 100 |cut -c2-`
        outfile=wrfout_${dm}_`wrf_time_string $outdate`
        watch_file $outfile 1 $rundir
        cp $outfile $WORK_DIR/fc/$DATE/wrfinput_${dm}_`wrf_time_string $outdate`
      done
    fi
  done
##save wrfout as deterministic forecast input (4dvar only)
  for n in `seq 1 $MAX_DOM`; do
    dm=d`expr $n + 100 |cut -c2-`
    outfile=wrfout_${dm}_`wrf_time_string $DATE`
    watch_file $outfile 1 $rundir
    cp $outfile $WORK_DIR/fc/$DATE/wrfinput_${dm}_`wrf_time_string $DATE`
  done
fi

if $CLEAN; then rm -f wrfout*; fi

echo complete > stat

