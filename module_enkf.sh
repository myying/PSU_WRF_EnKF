#!/bin/bash
. $CONFIG_FILE
rundir=$WORK_DIR/run/$DATE/enkf
if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi

cd $rundir
if [[ `cat stat` == "complete" ]]; then exit; fi

wait_for_module ../../$PREVDATE/wrf_ens #../obsproc
if [[ $JOB_SUBMIT_MODE == 1 ]]; then 
  wait_for_module ../icbc
  if $RUN_4DVAR; then  wait_for_module ../4dvar ../wrf_window; fi
fi

#Run EnKF
echo running > stat

domlist=`seq 1 $MAX_DOM`

for n in $domlist; do
  dm=d`expr $n + 100 |cut -c2-`
  if [[ ! -d $dm ]]; then mkdir -p $dm; fi
  if [ -f $dm/${DATE}.finish_flag ]; then continue; fi
  cd $dm

  #link priors
  for NE in `seq 1 $NUM_ENS`; do
    id=`expr $NE + 1000 |cut -c2-`
    ln -fs ../../../../fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE`_$id fort.`expr 80010 + $NE`
    cp -L fort.`expr 80010 + $NE` fort.`expr 90010 + $NE` >> link.log 2>&1 &
  done
  wait

  #ln -fs ../../../../../truth/wrfout_${dm}_`wrf_time_string $DATE` fort.70010
  cp -L fort.80011 fort.`expr 80011 + $NUM_ENS`
  cp -L fort.80011 fort.`expr 90011 + $NUM_ENS`
  ln -fs $ENKF_DIR/enkf.mpi .
  ln -fs $WRF_DIR/run/LANDUSE.TBL .

#  #link observations
#  #LITTLE_R format from obsproc
#  ln -fs $WORK_DIR/obs/$DATE/obs_gts_`wrf_time_string $DATE`.3DVAR obs_3dvar_${DATE}00
#  #airborne radar superobs
#  ln -fs $DATA_DIR/so/${DATE}_all.so_ass airborne_${DATE}_so

	ln -fs $OBS_DIR/obs_gts_`wrf_time_string $DATE`.3DVAR obs_3dvar_${DATE}00

#radiance obs
  ln -fs $WORK/code/CRTM/crtm_wrf/coefficients
  ln -fs $OBS_DIR/Met7/ch3/Tb_d01_${DATE}_so radiance_${DATE}_so

  $SCRIPT_DIR/namelist_enkf.sh $n > namelist.enkf
  cd ..
done

#run enkf.mpi
echo "  Running EnKF..."
tid=0
nn=$((($enkf_ntasks+$enkf_ppn-$enkf_ntasks%$enkf_ppn)/$enkf_ppn))
nt=$(($total_ntasks/$HOSTPPN/$nn))
for n in $domlist; do
  dm=d`expr $n + 100 |cut -c2-`
  if [ -f $dm/${DATE}.finish_flag ]; then continue; fi
  cd $dm
  $SCRIPT_DIR/job_submit.sh $enkf_ntasks $((tid*$enkf_ntasks)) $enkf_ppn ./enkf.mpi >& enkf.log &
#  mpirun.lsf ./enkf.mpi >& enkf.log
  tid=$((tid+1))
  if [[ $tid == $nt ]]; then
    tid=0
    wait
  fi
  cd ..
done
wait

#Check output
for n in $domlist; do
  dm=d`expr $n + 100 |cut -c2-`
  cd $dm
  watch_log enkf.log Successful 5 $rundir

  #Replace mean
  #1. replacing mean with 4DVar analysis (recentering) if running hybrid DA
  if $RUN_4DVAR; then
    wait_for_module ../../4dvar ../../wrf_window
    if [[ ! -d replace_mean ]]; then mkdir -p replace_mean; fi 
    cd replace_mean
    echo "  Replacing ens mean with 4DVar analysis for domain $dm"
    for NE in `seq 1 $((NUM_ENS+1))`; do
      mv ../fort.`expr 90010 + $NE` fort.`expr 80010 + $NE`
      cp fort.`expr 80010 + $NE` fort.`expr 90010 + $NE`
    done
    ln -fs $WORK_DIR/fc/$DATE/wrfinput_${dm}_`wrf_time_string $DATE` fort.70010
    ln -fs $ENKF_DIR/replace_mean.exe .
    ./replace_mean.exe $NUM_ENS >& replace_mean.log 
    watch_log replace_mean.log Successful 1 $rundir
    for NE in `seq 1 $((NUM_ENS+1))`; do
      mv fort.`expr 90010 + $NE` ../ 
    done
    cd ..
  fi

  #2. replacing mean with first guess (GFS/FNL) reanalysis
  if $REPLACE_MEAN; then
    if [[ ! -d replace_mean ]]; then mkdir -p replace_mean; fi
    cd replace_mean
    echo "  Replacing ens mean with $REPLACE_MEAN_WITH for domain $dm"
    for NE in `seq 1 $((NUM_ENS+1))`; do
      mv ../fort.`expr 90010 + $NE` fort.`expr 80010 + $NE`
      cp fort.`expr 80010 + $NE` fort.`expr 90010 + $NE`
    done
    if [[ $REPLACE_MEAN_WITH == "forecast" ]]; then
      ln -fs $WORK_DIR/fc/$DATE/wrfinput_$dm fort.70010
    elif [[ $REPLACE_MEAN_WITH == "gfs" ]]; then
      ln -fs $WORK_DIR/rc/$DATE/wrfinput_$dm fort.70010
    fi
    ln -fs $ENKF_DIR/replace_mean.exe .
    ./replace_mean.exe $NUM_ENS >& replace_mean.log
    watch_log replace_mean.log Successful 1 $rundir
    for NE in `seq 1 $((NUM_ENS+1))`; do
      mv fort.`expr 90010 + $NE` ../
    done
    cd ..
  fi

	for NE in `seq 1 $NUM_ENS`; do
		id=`expr $NE + 1000 |cut -c2-`
		mv fort.`expr 90010 + $NE` $WORK_DIR/fc/$DATE/wrfinput_${dm}_$id
		ln -fs ../../../../fc/$DATE/wrfinput_${dm}_$id fort.`expr 90010 + $NE`
	done
	cp fort.`expr 90011 + $NUM_ENS` $WORK_DIR/fc/$DATE/wrfinput_${dm}_mean
	cp fort.`expr 80011 + $NUM_ENS` $WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE`_mean
  cp fort.10000 $WORK_DIR/obs/$DATE/assimilated_obs_$dm
  cd ..
done

echo complete > stat

