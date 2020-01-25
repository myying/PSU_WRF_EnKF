#!/bin/bash
. $CONFIG_FILE
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

export start_date=$DATE
export run_minutes=$run_minutes_cycle
export minute_off=`echo "(${DATE:8:2}*60+${DATE:10:2})%$LBC_INTERVAL" |bc`
export LBDATE=`advance_time $DATE -$minute_off`

echo "  Running WRFforecast..."

for r in 1; do
  export time_step_ratio=$r
    touch rsl.error.0000
    if [[ `tail -n5 rsl.error.0000 |grep SUCCESS` ]]; then continue; fi

    if [[ $DATE == $LBDATE ]]; then
      export sst_update=1
    else
      export sst_update=0
    fi

    ####Running model
    ln -fs $WRF_DIR/run/* .
    rm -f namelist.*

    for n in `seq 1 $MAX_DOM`; do
      dm=d`expr $n + 100 |cut -c2-`
      #if $RUN_ENKF; then
      #  ln -fs ../../../fc/$DATE/wrfinput_${dm}_mean wrfinput_$dm
      #else
      #  ln -fs ../../../fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE`_mean wrfinput_$dm
      #fi
      ncea ../../../fc/$DATE/wrfinput_${dm}_??? wrfinput_$dm
    done
		ln -fs ../../../fc/wrfbdy_d01 wrfbdy_d01

    if [[ $SST_UPDATE == 1 ]]; then
      ln -fs ../../../rc/$DATE_START/wrflowinp_d?? .
    fi

    if $FOLLOW_STORM; then
      $SCRIPT_DIR/calc_ij_parent_start.sh $DATE ij_parent_start >& follow_storm.log
      watch_file ij_parent_start 1 $rundir
      $SCRIPT_DIR/calc_domain_moves.sh $DATE `advance_time $DATE_END -$TCV_INTERVAL` domain_moves >& follow_storm.log
      watch_file domain_moves 1 $rundir
    fi
    $SCRIPT_DIR/namelist_wrf.sh wrf > namelist.input
    $SCRIPT_DIR/job_submit.sh $wrf_ntasks 0 $HOSTPPN ./wrf.exe >& wrf.log
done

#Check outputs
  watch_log rsl.error.0000 SUCCESS 1 $rundir
  #if $RUN_4DVAR; then
  #  e_offset="$OBS_WIN_MIN 0"
  #else
  #  e_offset="0"
  #fi
  #for i in $e_offset; do
  #  outdate=`advance_time $NEXTDATE $i`
  #  outfile=wrfinput_d01_`wrf_time_string $outdate`
  #  mv $outfile $WORK_DIR/fc/$DATE/wrfinput_d01_`wrf_time_string $outdate`_mean
  #  if [ $MAX_DOM -gt 1 ]; then
  #    for n in `seq 1 $MAX_DOM`; do
  #      dm=d`expr $n + 100 |cut -c2-`
  #      outfile=wrfout_${dm}_`wrf_time_string $outdate`
  #      mv $outfile $WORK_DIR/fc/$DATE/wrfinput_${dm}_`wrf_time_string $outdate`_mean
  #    done
  #  fi
  #done

echo complete > stat
