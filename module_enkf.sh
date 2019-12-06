#!/bin/bash
. $CONFIG_FILE
rundir=$WORK_DIR/run/$DATE/enkf
if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi

cd $rundir
if [[ `cat stat` == "complete" ]]; then exit; fi

wait_for_module ../../$PREVDATE/wrf_ens ../obsproc ../icbc
if [[ $JOB_SUBMIT_MODE == 1 ]]; then
  wait_for_module ../icbc
  if $RUN_4DVAR; then  wait_for_module ../4dvar ../wrf_window; fi
fi

#Run EnKF
echo running > stat
echo "  Running EnKF..."

domlist=`seq 1 $MAX_DOM`

###preparing files
for n in $domlist; do
  dm=d`expr $n + 100 |cut -c2-`
  if [[ ! -d $dm ]]; then mkdir -p $dm; fi
  if [ -f $dm/${DATE}.finish_flag ]; then continue; fi
  cd $dm

  #link prior members
  for NE in `seq 1 $NUM_ENS`; do
    id=`expr $NE + 1000 |cut -c2-`
    ln -fs ../../../../fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE`_$id fort.`expr 80010 + $NE`
  done
  cp -L fort.80011 fort.`expr 80011 + $NUM_ENS` ##prior mean

  #make several copies
  for NE in `seq 1 $((NUM_ENS+1))`; do
    if [ $NUM_SCALES == 1 ]; then
      ln -fs fort.`expr 80010 + $NE` fort.`expr 50010 + $NE`
      ln -fs fort.`expr 80010 + $NE` fort.`expr 90010 + $NE`
    else
      cp -L fort.`expr 80010 + $NE` fort.`expr 50010 + $NE`
      cp -L fort.`expr 80010 + $NE` fort.`expr 60010 + $NE`
      cp -L fort.`expr 80010 + $NE` fort.`expr 90010 + $NE`
    fi
    cp -L  fort.`expr 80010 + $NE` fort.`expr 70010 + $NE`
  done

  ln -fs $ENKF_DIR/enkf.mpi .
  ln -fs $WRF_DIR/run/LANDUSE.TBL .

  if [ $NUM_SCALES -gt 1 ]; then
    ln -fs $ENKF_DIR/scale_decompose.exe .
    ln -fs $ENKF_DIR/alignment.exe .
  fi

  ##link observations
  #LITTLE_R format from obsproc
  ln -fs $WORK_DIR/obs/$DATE/obs_gts_`wrf_time_string $DATE`.3DVAR obs_3dvar_${DATE}00

  #airborne radar superobs
  if [ -f $DATA_DIR/airborne_radar/SO/${DATE}_all.so_ass ]; then
    ln -fs $DATA_DIR/airborne_radar/SO/${DATE}_all.so_ass airborne_${DATE}_so
  else
    echo > airborne_${DATE}_so
  fi

  #radiance obs
  ln -fs $WORK/code/CRTM/crtm_wrf/coefficients
  if [ -f $DATA_DIR/radiance/SO/radiance_d03_${DATE}_so ]; then
    ln -fs $DATA_DIR/radiance/SO/radiance_d03_${DATE}_so radiance_${DATE}_so
  else
    echo > radiance_${DATE}_so
  fi

  cd ..
done

if [ $NUM_SCALES == 1 ]; then
  for n in $domlist; do
    dm=d`expr $n + 100 |cut -c2-`
    if [[ ! -d $dm ]]; then mkdir -p $dm; fi
    if [ -f $dm/${DATE}.finish_flag ]; then continue; fi
    cd $dm
    $SCRIPT_DIR/namelist_enkf.sh $n 1 > namelist.enkf
    $SCRIPT_DIR/job_submit.sh $enkf_ntasks 0 $enkf_ppn ./enkf.mpi >& enkf.log
    watch_log enkf.log Successful 5 $rundir
    cd ..
  done
else
  ###runing enkf.mpi multiscale scheme
  for s in `seq 1 $NUM_SCALES`; do
    dm=d01
    if [[ ! -d $dm ]]; then mkdir -p $dm; fi
    cd $dm
    if [[ ! -d scale$s ]]; then mkdir -p scale$s; fi
    for n in `seq 2 $MAX_DOM`; do
      dm1=d`expr $n + 100 |cut -c2-`
      mkdir -p ../$dm1/scale$s
    done
    if [ -f scale$s/${DATE}.finish_flag ]; then continue; fi
    $SCRIPT_DIR/namelist_enkf.sh 1 $s > namelist.enkf
    $SCRIPT_DIR/job_submit.sh $NUM_ENS 0 $enkf_ppn ./scale_decompose.exe >& scale_decompose.log
    watch_log scale_decompose.log Successful 5 $rundir
    $SCRIPT_DIR/job_submit.sh $enkf_ntasks 0 $enkf_ppn ./enkf.mpi >& enkf.log
    watch_log enkf.log Successful 5 $rundir
    $SCRIPT_DIR/job_submit.sh $NUM_ENS 0 $enkf_ppn ./alignment.exe >& alignment.log
    watch_log alignment.log Successful 5 $rundir
    ##save copy
    mv ${DATE}.finish_flag scale$s/.
    cp fort.5* fort.6* fort.7* fort.9* scale$s/.
    for n in `seq 2 $MAX_DOM`; do
      dm1=d`expr $n + 100 |cut -c2-`
      cp ../$dm1/fort.9* ../$dm1/scale$s/.
    done
    cd ..
  done
  for n in `seq 2 $MAX_DOM`; do
    dm=d`expr $n + 100 |cut -c2-`
    cd $dm
    $SCRIPT_DIR/namelist_enkf.sh $n 1 > namelist.enkf
    $SCRIPT_DIR/job_submit.sh $enkf_ntasks 0 $enkf_ppn ./enkf.mpi >& enkf.log
    watch_log enkf.log Successful 5 $rundir
    cd ..
  done
fi

#Check output
for n in $domlist; do
  dm=d`expr $n + 100 |cut -c2-`
  cd $dm

  #Replace mean
  #1. replacing mean with 4DVar analysis (recentering) if running hybrid DA
  #if $RUN_4DVAR; then
  #  wait_for_module ../../4dvar ../../wrf_window
  #  if [[ ! -d replace_mean ]]; then mkdir -p replace_mean; fi
  #  cd replace_mean
  #  echo "  Replacing ens mean with 4DVar analysis for domain $dm"
  #  for NE in `seq 1 $((NUM_ENS+1))`; do
  #    mv ../fort.`expr 90010 + $NE` fort.`expr 80010 + $NE`
  #    cp fort.`expr 80010 + $NE` fort.`expr 90010 + $NE`
  #  done
  #  ln -fs $WORK_DIR/fc/$DATE/wrfinput_${dm}_`wrf_time_string $DATE` fort.70010
  #  ln -fs $ENKF_DIR/replace_mean.exe .
  #  ./replace_mean.exe $NUM_ENS >& replace_mean.log
  #  watch_log replace_mean.log Successful 1 $rundir
  #  for NE in `seq 1 $((NUM_ENS+1))`; do
  #    mv fort.`expr 90010 + $NE` ../
  #  done
  #  cd ..
  #fi

  ###2. replacing mean with first guess (GFS/FNL) reanalysis
  #if [[ ! -d replace_mean ]]; then mkdir -p replace_mean; fi
  #cd replace_mean
  #echo "  Replacing ens mean with $REPLACE_MEAN_WITH for domain $dm"
  #for NE in `seq 1 $((NUM_ENS+1))`; do
  #  mv ../fort.`expr 90010 + $NE` fort.`expr 80010 + $NE`
  #  cp fort.`expr 80010 + $NE` fort.`expr 90010 + $NE`
  #done
  #if [[ $REPLACE_MEAN_WITH == "forecast" ]]; then
  #  ln -fs $WORK_DIR/fc/$DATE/wrfinput_$dm fort.70010
  #elif [[ $REPLACE_MEAN_WITH == "gfs" ]]; then
  #  ln -fs $WORK_DIR/rc/$DATE/wrfinput_$dm fort.70010
  #fi
  #ln -fs $ENKF_DIR/replace_mean.exe .
  #./replace_mean.exe $NUM_ENS >& replace_mean.log
  #watch_log replace_mean.log Successful 1 $rundir
  #for NE in `seq 1 $((NUM_ENS+1))`; do
  #  mv fort.`expr 90010 + $NE` ../
  #done
  #cd ..

  for NE in `seq 1 $NUM_ENS`; do
    id=`expr $NE + 1000 |cut -c2-`
    if [ $NUM_SCALES == 1 ]; then
      mv fort.`expr 70010 + $NE` $WORK_DIR/fc/$DATE/wrfinput_${dm}_$id
    else
      if [ $n == 1 ]; then
        mv fort.`expr 90010 + $NE` $WORK_DIR/fc/$DATE/wrfinput_${dm}_$id
      else
        mv fort.`expr 70010 + $NE` $WORK_DIR/fc/$DATE/wrfinput_${dm}_$id
      fi
    fi
  done
  if [ $NUM_SCALES == 1 ]; then
    cp fort.`expr 70011 + $NUM_ENS` $WORK_DIR/fc/$DATE/wrfinput_${dm}_mean
  else
    if [ $n == 1 ]; then
      cp fort.`expr 90011 + $NUM_ENS` $WORK_DIR/fc/$DATE/wrfinput_${dm}_mean
    else
      cp fort.`expr 70011 + $NUM_ENS` $WORK_DIR/fc/$DATE/wrfinput_${dm}_mean
    fi
  fi
  cd ..
done

echo complete > stat

