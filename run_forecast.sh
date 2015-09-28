#!/bin/bash
#####header for jet######
#PBS -A hfip-psu
#PBS -N BAMEX_EnKF
#PBS -l walltime=8:00:00
#PBS -q batch
#PBS -l partition=ujet:tjet
#PBS -l procs=720
#PBS -j oe
#PBS -o ./log
#PBS -d .

#####header for stampede######
##SBATCH -J run_cycle
##SBATCH -n 336
##SBATCH -p normal
##SBATCH -t 2:00:00

source ~/.bashrc

#load configuration files, functions, parameters
cd $WORK/DA
export CONFIG_FILE=$WORK/DA/config/dynamo_osse1
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

#  # ICBC
#  $SCRIPT_DIR/module_icbc.sh &
#  # Ensemble initialization and forecast
#  if $RUN_ENKF; then
#    if [ $DATE == $DATE_START ]; then
#      $SCRIPT_DIR/module_perturb_ic.sh &
#    fi
#    if [ $NEXTDATE -le $DATE_CYCLE_END ]; then
#      $SCRIPT_DIR/module_wrf_ens.sh &
#    fi
#  fi
#  # First deterministic run for 4DVar
#  if $RUN_4DVAR && ! $RUN_ENKF && [ $DATE == $DATE_START ]; then
#    $SCRIPT_DIR/module_wrf_window.sh &
#  fi

  # Data assimilation for each cycle
  if [ $DATE -ge $DATE_CYCLE_START ] && [ $DATE -le $DATE_CYCLE_END ]; then
    # Processing observations
#    if $RUN_ENKF || $RUN_4DVAR; then
#      $SCRIPT_DIR/module_obsproc.sh &
#    fi
#    # EnKF
#    if $RUN_ENKF; then
#      $SCRIPT_DIR/module_enkf.sh &
#    fi
#    # 4DVar
#    if $RUN_4DVAR; then
#      $SCRIPT_DIR/module_wrf_window1.sh &
#      $SCRIPT_DIR/module_4dvar.sh &
#      $SCRIPT_DIR/module_wrf_window.sh &
#    fi
#    # EnVar need an extra ensemble run through the obs window to get perturbations
#    if $RUN_ENVAR; then
#      $SCRIPT_DIR/module_wrf_ens_window1.sh &
#    fi
     $SCRIPT_DIR/module_wrf.sh &
  fi
#wait
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

