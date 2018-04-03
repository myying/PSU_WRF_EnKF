#!/bin/bash
#BSUB -P UPSU0001
#BSUB -J run_obsproc
#BSUB -W 2:00
#BSUB -q small
#BSUB -n 16
#BSUB -R "span[ptile=16]"
#BSUB -o log/%J.out
#BSUB -e log/%J.err
source /glade/u/apps/opt/lmod/4.2.1/init/bash
source ~/.bashrc

#load configuration files, functions, parameters
cd $WORK/PSU_WRF_EnKF
export CONFIG_FILE=$WORK/PSU_WRF_EnKF/config/EnKF_OSSE/obsproc
. $CONFIG_FILE
. util.sh

if [[ ! -d $WORK_DIR ]]; then mkdir -p $WORK_DIR; fi
cd $WORK_DIR

if [ $JOB_SUBMIT_MODE == 1 ]; then
  if [[ $HOSTTYPE == "stampede" ]]; then
    export total_ntasks=$SLURM_NTASKS
  fi
  if [[ $HOSTTYPE == "jet" ]]; then
    export total_ntasks=$PBS_NP
  fi
  if [[ $HOSTTYPE == "yellowstone" ]]; then
    export total_ntasks=$LSB_MAX_NUM_PROCESSORS
  fi
else
  export total_ntasks=9999999
fi

#start cycling
date
export DATE=$DATE_START
export PREVDATE=$DATE
export NEXTDATE=$DATE

while [[ $NEXTDATE -le $DATE_CYCLE_END ]]; do  #CYCLE LOOP

  #TIME CALCULATION
  if [[ $DATE == $DATE_START ]]; then
    export run_minutes_cycle=`diff_time $DATE $DATE_CYCLE_START`
  else
    export run_minutes_cycle=$CYCLE_PERIOD
  fi
  export NEXTDATE=`advance_time $DATE $run_minutes_cycle`
  if $FORECAST_TO_END; then
    export run_minutes_forecast=`diff_time $DATE $DATE_END`
  else
    export run_minutes_forecast=`max $run_minutes_cycle $FORECAST_MINUTES`
  fi
  #LBDATE = Closest to DATE when FG is available
  export minute_off=`echo "(${DATE:8:2}*60+${DATE:10:2})%$LBC_INTERVAL" |bc`
  export LBDATE=`advance_time $DATE -$minute_off` 

  echo "----------------------------------------------------------------------"
  echo "CURRENT CYCLE: `wrf_time_string $DATE` => `wrf_time_string $NEXTDATE`"
  mkdir -p {run,rc,fc,output,obs}/$DATE

  #CLEAR ERROR TAGS
  for d in `ls run/$DATE/`; do
    if [[ `cat run/$DATE/$d/stat` != "complete" ]]; then
      echo waiting > run/$DATE/$d/stat
    fi
  done

  #RUN COMPONENTS---------------------------------------

  # Data assimilation for each cycle
  if [ $DATE -ge $DATE_CYCLE_START ] && [ $DATE -le $DATE_CYCLE_END ]; then
    # Processing observations
    if $RUN_ENKF || $RUN_4DVAR; then
      $SCRIPT_DIR/module_obsproc.sh &
    fi
  fi
  wait

  #CHECK ERRORS
  for d in `ls -t run/$DATE/`; do
    if [[ `cat run/$DATE/$d/stat` == "error" ]]; then
      echo CYCLING STOP DUE TO FAILED COMPONENT: $d
      exit 1
    fi
  done

  #ADVANCE TO NEXT CYCLE
  export PREVDATE=$DATE
  export DATE=$NEXTDATE
done
echo CYCLING COMPLETE

