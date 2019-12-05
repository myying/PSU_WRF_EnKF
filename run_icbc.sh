#!/bin/bash
#SBATCH -J run_cycle
#SBATCH -n 32 -N 1
#SBATCH -p development
#SBATCH -t 2:00:00
#SBATCH -o out
source ~/.bashrc

#load configuration files, functions, parameters
cd $WORK/PSU_WRF_EnKF
export CONFIG_FILE=$WORK/PSU_WRF_EnKF/config/Patricia_test
. $CONFIG_FILE
. util.sh

if [[ ! -d $WORK_DIR ]]; then mkdir -p $WORK_DIR; fi
cd $WORK_DIR

if [ $JOB_SUBMIT_MODE == 1 ]; then
  if [[ $HOSTTYPE == "stampede2" ]]; then
    export total_ntasks=$SLURM_NTASKS
  fi
  if [[ $HOSTTYPE == "jet" ]]; then
    export total_ntasks=$PBS_NP
  fi
  if [[ $HOSTTYPE == "cheyenne" ]]; then
    export total_ntasks=$NCPUS
  fi
else
  export total_ntasks=9999999
fi
echo total_ntasks=$total_ntasks

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
  #LBDATE: Closest to DATE when LBC is available in wrfbdy/wrflowinp
  export minute_off=`echo "(${DATE:8:2}*60+${DATE:10:2})%$LBC_INTERVAL" |bc`
  export LBDATE=`advance_time $DATE -$minute_off`

  echo "----------------------------------------------------------------------"
  echo "CURRENT CYCLE: `wrf_time_string $DATE` => `wrf_time_string $NEXTDATE`"
  mkdir -p {run,rc,fc,output,obs}/$DATE

  #CLEAR ERROR TAGS
  for d in `ls run/$DATE/`; do
    if [[ `cat run/$DATE/$d/stat` == "error" ]]; then
      echo waiting > run/$DATE/$d/stat
    fi
  done

  #RUN COMPONENTS---------------------------------------

  # ICBC
  $SCRIPT_DIR/module_icbc.sh

  #CHECK ERRORS
  for d in `ls -t run/$DATE/`; do
    if [[ `cat run/$DATE/$d/stat` == "error" ]]; then
      echo CYCLING STOP DUE TO FAILED COMPONENT: $d
      date
      exit 1
    fi
  done

  #ADVANCE TO NEXT CYCLE
  export PREVDATE=$DATE
  export DATE=$NEXTDATE
done
echo CYCLING COMPLETE
date

