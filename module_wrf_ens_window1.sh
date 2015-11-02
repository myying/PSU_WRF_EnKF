#!/bin/bash
. $CONFIG_FILE
rundir=$WORK_DIR/run/$DATE/wrf_ens_window1
if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi

cd $rundir
if [[ `cat stat` == "complete" ]]; then exit; fi

#Check dependency
wait_for_module ../icbc ../../$DATE_START/perturb_ic
wait_for_module ../../$PREVDATE/wrf_ens
if [[ $JOB_SUBMIT_MODE == 1 ]]; then wait_for_module ../wrf_window1; fi

echo running > stat

export start_date=`advance_time $DATE $OBS_WIN_MIN`
export run_minutes=`echo "$OBS_WIN_MAX - $OBS_WIN_MIN" |bc`
export inputout_interval=$MINUTES_PER_SLOT
export inputout_begin=0
export inputout_end=$run_minutes
export wrfout_interval=$inputout_interval

#1 Perturb BCs with random_samples
echo "  Perturbing ensemble BCs..."
dd=`diff_time $DATE_START $LBDATE`
n_1=$((dd/$LBC_INTERVAL+1))

tid=0
nt=$total_ntasks
for NE in `seq 1 $NUM_ENS`; do
  id=`expr $NE + 1000 |cut -c2-`
  if [[ ! -d $id ]]; then mkdir $id; fi
  touch $id/update_wrf_bc.log
  if [[ `tail -n5 $id/update_wrf_bc.log |grep successfully` ]]; then continue; fi
  cd $id

  cat > parame.in << EOF
&control_param
 wrf_3dvar_output_file = 'wrfinput_d01_update'
 wrf_bdy_file          = 'wrfbdy_d01_update'
 wrf_bdy_file_real     = 'wrfbdy_d01_real'
 wrf_input_from_si     = 'wrfinput_d01_real'
 wrf_input_from_si_randmean = 'random_mean'
 wrf_3dvar_random_draw = 'random_draw'
 cycling = .true.
 debug   = .true. 
 low_bdy_only = .false. 
 perturb_bdy = .true.
 n_1 = $n_1
/
EOF

  ln -fs $WRF_BC_DIR/update_wrf_bc.exe .

  if [ $DATE == $DATE_CYCLE_START ]; then
    ln -fs ../../../../fc/wrfbdy_d01 wrfbdy_d01_real
  else
    ln -fs ../../../../fc/wrfbdy_d01_window_$id wrfbdy_d01_real
  fi
  cp -L wrfbdy_d01_real wrfbdy_d01_update

  ln -fs ../../../../rc/$DATE/wrfinput_d01 wrfinput_d01_real
  rm -f wrfinput_d01_update
  ln -fs ../../../../fc/$PREVDATE/wrfinput_d01_`wrf_time_string $start_date`_$id wrfinput_d01_update

  ln -fs ../../../../fc/$DATE_START/wrfinput_d01_`expr $((RANDOM%($NUM_ENS-1)+1)) + 1000 |cut -c2-` random_draw
  ln -fs ../../../../rc/$DATE_START/wrfinput_d01 random_mean

#  $SCRIPT_DIR/job_submit.sh 1 $tid $HOSTPPN ./update_wrf_bc.exe >& update_wrf_bc.log &
  ./update_wrf_bc.exe >& update_wrf_bc.log

  tid=$((tid+1))
  if [[ $tid == $nt ]]; then
    tid=0
    wait
  fi
  cd ..
done
wait

for NE in `seq 1 $NUM_ENS`; do
  id=`expr $NE + 1000 |cut -c2-`
  watch_log $id/update_wrf_bc.log successfully 1 $rundir
  mv $id/wrfbdy_d01_update $WORK_DIR/fc/wrfbdy_d01_window_$id
done


#Run wrf forecast for each member
echo "  Running WRF ensemble forecast through the obs window (fg -> fg02)..."

tid=0
nt=$((total_ntasks/$wrf_ntasks))
for NE in `seq 1 $NUM_ENS`; do
  id=`expr $NE + 1000 |cut -c2-`
  touch $id/rsl.error.0000
  if [[ `tail -n5 $id/rsl.error.0000 |grep SUCCESS` ]]; then continue; fi

  cd $id
  ln -fs $WRF_DIR/run/* .
  rm -f namelist.*

  for n in `seq 1 $MAX_DOM`; do
    dm=d`expr $n + 100 |cut -c2-`
    ln -fs ../../../../fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $start_date`_$id wrfinput_$dm
  done
  ln -fs ../../../../fc/wrfbdy_d01_window_$id wrfbdy_d01

  if [[ $SST_UPDATE == 1 ]]; then
    ln -fs ../../../../rc/$DATE_START/wrflowinp_d?? .
  fi

  if $FOLLOW_STORM; then
    cp $WORK_DIR/rc/$DATE/ij_parent_start_4dvar ij_parent_start
  fi
  $SCRIPT_DIR/namelist_wrf.sh wrf > namelist.input
  $SCRIPT_DIR/job_submit.sh $wrf_ntasks $((tid*$wrf_ntasks)) $HOSTPPN ./wrf.exe >& wrf.log &
  tid=$((tid+1))
  if [[ $tid == $nt ]]; then
    tid=0
    wait
  fi
  cd ..
done
wait

#Check outputs
for NE in `seq 1 $NUM_ENS`; do
  id=`expr $NE + 1000 |cut -c2-`
  watch_log $id/rsl.error.0000 SUCCESS 1 $rundir
  e_offset=`seq $OBS_WIN_MIN $MINUTES_PER_SLOT $OBS_WIN_MAX`
  for i in $e_offset; do
    outdate=`advance_time $DATE $i`
    if [ $outdate -eq $start_date ]; then
      for n in `seq 1 $MAX_DOM`; do
        dm=d`expr $n + 100 |cut -c2-`
        outfile=$WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $start_date`_$id
        ln -fs $outfile $WORK_DIR/fc/$DATE/wrfinput_${dm}_`wrf_time_string $start_date`_window_$id
      done
    else
      outfile=$id/wrfinput_d01_`wrf_time_string $outdate`
      mv $outfile $WORK_DIR/fc/$DATE/wrfinput_d01_`wrf_time_string $outdate`_window_$id
      if [ $MAX_DOM -gt 1 ]; then
        for n in `seq 2 $MAX_DOM`; do
          dm=d`expr $n + 100 |cut -c2-`
          outfile=$id/wrfout_${dm}_`wrf_time_string $outdate`
          mv $outfile $WORK_DIR/fc/$DATE/wrfinput_${dm}_`wrf_time_string $outdate`_window_$id
        done
      fi
    fi
  done
done

if $CLEAN; then 
  for NE in `seq 1 $NUM_ENS`; do
    id=`expr $NE + 1000 |cut -c2-`
    rm -f $rundir/$id/wrfout*
  done
fi

echo complete > stat
