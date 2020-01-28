#!/bin/bash
source ~/.bashrc

#load configuration files, functions, parameters
cd $WORK/PSU_WRF_EnKF
export CONFIG_FILE=$WORK/PSU_WRF_EnKF/config/Patricia/control_cheyenne
. $CONFIG_FILE
. util.sh

if [[ ! -d $WORK_DIR ]]; then mkdir -p $WORK_DIR; fi
cd $WORK_DIR

export total_ntasks=9999999

date
export DATE=201510210600
export run_minutes_cycle=`diff_time $DATE $DATE_END`
$SCRIPT_DIR/module_wrf.sh &
