#!/bin/bash
. $WORK/PSU_WRF_EnKF/util.sh

for t in 201110{16..18}*; do
  if [ ! -f $t/wrfout_d01_`wrf_time_string $t`_mean ]; then
    ncra -h $t/wrfout_d01_`wrf_time_string $t`_{001..060} $t/wrfout_d01_`wrf_time_string $t`_mean
  fi
  t1=`advance_time $t 180`
  if [ $t/wrfout_d01_`wrf_time_string $t1`_mean ]; then
    ncra -h $t/wrfout_d01_`wrf_time_string $t1`_{001..060} $t/wrfout_d01_`wrf_time_string $t1`_mean
  fi
done
