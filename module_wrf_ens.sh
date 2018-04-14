#!/bin/bash
. $CONFIG_FILE
rundir=$WORK_DIR/run/$DATE/wrf_ens
if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi

cd $rundir
if [[ `cat stat` == "complete" ]]; then exit; fi

#Check dependency
wait_for_module ../icbc 
if [ $DATE == $DATE_START ]; then wait_for_module ../perturb_ic; fi
if [ $DATE -gt $DATE_START ]; then
  if $RUN_ENKF; then wait_for_module ../enkf; fi
fi

echo running > stat

export start_date=$DATE
export run_minutes=$run_minutes_cycle 
if $RUN_4DVAR; then
  export inputout_interval=`abs $OBS_WIN_MIN`
  export inputout_begin=`echo "$run_minutes + $OBS_WIN_MIN" |bc`
  export inputout_end=$run_minutes
  export wrfout_interval=`gcd $inputout_begin $inputout_end`
else
  export inputout_interval=$run_minutes
  export inputout_begin=0
  export inputout_end=$run_minutes
fi

echo "  Running WRF ensemble forecast..."

tid=0
nt=$((total_ntasks/$wrf_ntasks))
for r in 1 1.5; do
  export time_step_ratio=$r
  for NE in `seq 1 $NUM_ENS`; do
    id=`expr $NE + 1000 |cut -c2-`
    if [[ ! -d $id ]]; then mkdir $id; fi
    touch $id/rsl.error.0000
    if [[ `tail -n5 $id/rsl.error.0000 |grep SUCCESS` ]]; then continue; fi
  
    cd $id

    ####Update and perturb BC
    #if [ $DATE == $LBDATE ] && [ $r == 1 ]; then
      #dd=`diff_time $DATE_START $LBDATE`
      #n_1=$((dd/$LBC_INTERVAL+1))
      #cat > parame.in << EOF
#&control_param
 #wrf_3dvar_output_file = 'wrfinput_d01_update'
 #wrf_bdy_file          = 'wrfbdy_d01_update'
 #wrf_bdy_file_real     = 'wrfbdy_d01_real'
 #wrf_input_from_si     = 'wrfinput_d01_real'
 #wrf_input_from_si_randmean = 'random_mean'
 #wrf_3dvar_random_draw = 'random_draw'
 #cycling = .true.
 #low_bdy_only = .false. 
 #perturb_bdy = .true.
 #n_1 = $n_1
#/
#EOF
      #ln -fs $WRF_BC_DIR/update_wrf_bc.exe .

      #ln -fs ../../../../fc/wrfbdy_d01_$id wrfbdy_d01_real
      #cp -L wrfbdy_d01_real wrfbdy_d01_update
      #ln -fs ../../../../rc/$DATE_START/wrfinput_d01_`wrf_time_string $DATE` wrfinput_d01_real
      #ln -fs ../../../../fc/$DATE/wrfinput_d01_$id wrfinput_d01_update

      #randnum=`expr $((RANDOM%99+1)) + 1000 |cut -c2-`
      #ln -fs ../../../../fc/$DATE_START/wrfinput_d01_$randnum random_draw
      #ln -fs ../../../../fc/$DATE_START/wrfinput_d01 random_mean
      #echo $n_1 $randnum >> ../../../../fc/rand_$id

      #./update_wrf_bc.exe >& update_wrf_bc.log 
      #watch_log update_wrf_bc.log successfully 1 $rundir
      #mv wrfbdy_d01_update $WORK_DIR/fc/wrfbdy_d01_$id
    #fi
    if [ $DATE == $LBDATE ]; then
      export sst_update=1
    else
      export sst_update=0
    fi

    ####Running model
    ln -fs $WRF_DIR/run/* .
    rm -f namelist.*
  
    for n in `seq 1 $MAX_DOM`; do
      dm=d`expr $n + 100 |cut -c2-`
      if $RUN_ENKF; then
        ln -fs ../../../../fc/$DATE/wrfinput_${dm}_$id wrfinput_$dm
      else
        ln -fs ../../../../fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE`_$id wrfinput_$dm
      fi
    done
    #ln -fs ../../../../fc/wrfbdy_d01_$id wrfbdy_d01
		ln -fs ../../../../fc/wrfbdy_d01 wrfbdy_d01

    if $MULTI_PHYS_ENS; then
      if [ $DATE -gt $DATE_START ]; then
        $SCRIPT_DIR/multi_physics_set.sh $id >& multi_physics_set.log
      fi
    fi

    if [[ $SST_UPDATE == 1 ]]; then
      ln -fs ../../../../fc/wrflowinp_d?? .
    fi

    if $FOLLOW_STORM; then
      cp $WORK_DIR/rc/$DATE/ij_parent_start .
      cp $WORK_DIR/rc/$DATE/domain_moves .
    fi
    $SCRIPT_DIR/namelist_wrf.sh wrf > namelist.input
    $SCRIPT_DIR/job_submit.sh $wrf_ntasks $((tid*$wrf_ntasks)) $HOSTPPN ./wrf.exe >& wrf.log &
#    mpirun.lsf ./wrf.exe

    tid=$((tid+1))
    if [[ $tid == $nt ]]; then
      tid=0
      wait
    fi
    cd ..
  done
  wait
done

#Check outputs
for NE in `seq 1 $NUM_ENS`; do
  id=`expr $NE + 1000 |cut -c2-`
  watch_log $id/rsl.error.0000 SUCCESS 1 $rundir
  if $RUN_4DVAR; then
    e_offset="$OBS_WIN_MIN 0"
  else
    e_offset="0"
  fi
  for i in $e_offset; do
    outdate=`advance_time $NEXTDATE $i`
    outfile=$id/wrfinput_d01_`wrf_time_string $outdate`
    mv $outfile $WORK_DIR/fc/$DATE/wrfinput_d01_`wrf_time_string $outdate`_$id
    if [ $MAX_DOM -gt 1 ]; then
      for n in `seq 1 $MAX_DOM`; do
        dm=d`expr $n + 100 |cut -c2-`
        outfile=$id/wrfout_${dm}_`wrf_time_string $outdate`
        mv $outfile $WORK_DIR/fc/$DATE/wrfinput_${dm}_`wrf_time_string $outdate`_$id
      done
    fi
  done
done

for NE in `seq 1 $NUM_ENS`; do
  id=`expr $NE + 1000 |cut -c2-`
  if $RUN_ENKF; then
    mv $id/wrfout_d01_`wrf_time_string $DATE` $WORK_DIR/output/$DATE/wrfout_d01_`wrf_time_string $DATE`_$id
  fi
  mv $id/wrfout_d01_`wrf_time_string $NEXTDATE` $WORK_DIR/output/$DATE/wrfout_d01_`wrf_time_string $NEXTDATE`_$id
done

#Calculate ensemble mean for 4DVar fg (next cycle)
if $RUN_4DVAR; then
  echo "  Calculating ensemble mean..."
  cd $rundir
  for i in $OBS_WIN_MIN; do
    outdate=`advance_time $NEXTDATE $i`
    if [[ $outdate -ne $DATE ]]; then
      for n in `seq 1 $MAX_DOM`; do
        dm=d`expr $n + 100 |cut -c2-`
        if [[ ! -d mean_${dm}_$outdate ]]; then mkdir -p mean_${dm}_$outdate; fi
        cd mean_${dm}_$outdate
        for NE in `seq 1 $NUM_ENS`; do
          id=`expr $NE + 1000 |cut -c2-`
          ln -fs ../../../../fc/$DATE/wrfinput_${dm}_`wrf_time_string $outdate`_$id fort.`expr 80010 + $NE`
        done
        cp -L fort.80011 fort.`expr 80011 + $NUM_ENS`
        $SCRIPT_DIR/namelist_enkf.sh $n > namelist.enkf
        ln -fs $ENKF_DIR/ensemble_mean.exe .
        $SCRIPT_DIR/job_submit.sh $enkf_ntasks 0 $enkf_ppn ./ensemble_mean.exe >& ensemble_mean.log 
#        mpirun.lsf ./ensemble_mean.exe >& ensemble_mean.log

        watch_log ensemble_mean.log Successful 1 $rundir
        mv fort.`expr 80011 + $NUM_ENS` $WORK_DIR/fc/$DATE/wrfinput_${dm}_`wrf_time_string $outdate`_mean
        cd ..
      done
    else
      ln -sf $WORK_DIR/fc/$DATE/wrfinput_d01_mean $WORK_DIR/fc/$DATE/wrfinput_d01_`wrf_time_string $outdate`_mean
    fi
  done
fi

if $CLEAN; then 
  for NE in `seq 1 $NUM_ENS`; do
    id=`expr $NE + 1000 |cut -c2-`
    rm -f $rundir/$id/wrfout*
  done
fi

echo complete > stat
rm ???/rsl.* &

