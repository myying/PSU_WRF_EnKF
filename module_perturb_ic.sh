#!/bin/bash
. $CONFIG_FILE

rundir=$WORK_DIR/run/$DATE/perturb_ic
if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi

cd $rundir
if [[ `cat stat` == "complete" ]]; then exit; fi

#Check dependency
wait_for_module ../icbc

echo running > stat
echo "  Perturbing IC using WRF 3DVar..."

#Run randomcv wrfvar to generate 100 random perturbations
tid=0
nt=$((total_ntasks/$var3d_ntasks))
nm=$NUM_ENS

for i in `seq 1 $nm`; do 
  id=`expr $i + 1000 |cut -c2-`

  if [[ ! -d $id ]]; then mkdir $id; fi
  touch $id/rsl.error.0000
  if [[ `tail -n5 $id/rsl.error.0000 |grep successful` ]]; then continue; fi
  cd $id

  ln -fs $WRFDA_DIR/run/LANDUSE.TBL .
  ln -fs $WRFDA_DIR/var/build/da_wrfvar.exe .
  cp $WORK_DIR/fc/$DATE_START/wrfinput_d?? .
  ln -fs wrfinput_d01 fg
  if [[ $CV_OPTIONS == 3 ]]; then
    ln -fs $WRFDA_DIR/var/run/be.dat.cv3 be.dat
  else
    ln -fs $BE_DIR/be.dat .
  fi

  export analysis_type="RANDOMCV"
  export time_window_min=$DATE
  export time_window_max=$DATE
  $SCRIPT_DIR/namelist_wrfvar.sh > namelist.input
  export start_date=$DATE
  export run_minutes=0
  $SCRIPT_DIR/namelist_wrf.sh wrfvar 1 >> namelist.input

  $SCRIPT_DIR/job_submit.sh $var3d_ntasks $((tid*$var3d_ntasks)) $HOSTPPN ./da_wrfvar.exe >& da_wrfvar.log &

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
  watch_log $id/rsl.error.0000 successfully 1 $rundir
  mv $id/wrfvar_output $WORK_DIR/fc/$DATE/wrfinput_d01_$id
done

#recenter ensemble to first guess
if [[ ! -d replace_mean ]]; then mkdir -p replace_mean; fi
cd replace_mean
  for NE in `seq 1 $NUM_ENS`; do
    id=`expr $NE + 1000 |cut -c2-`
    ln -fs $WORK_DIR/fc/$DATE/wrfinput_d01_$id fort.`expr 80010 + $NE`
    cp fort.`expr 80010 + $NE` fort.`expr 90010 + $NE`
  done
  ln -fs $WORK_DIR/fc/$DATE_START/wrfinput_d01 fort.70010
  ln -fs $ENKF_DIR/replace_mean.exe .
  ./replace_mean.exe $NUM_ENS >& replace_mean.log
  watch_log replace_mean.log Successful 1 $rundir
  for NE in `seq 1 $NUM_ENS`; do
    id=`expr $NE + 1000 |cut -c2-`
    mv fort.`expr 90010 + $NE` $WORK_DIR/fc/$DATE/wrfinput_d01_$id
  done
cd ..

###
#for NE in `seq 1 $NUM_ENS`; do
#  id=`expr $NE + 1000 |cut -c2-`
#  cp $WORK_DIR/fc/$DATE_START/wrfinput_d01 $WORK_DIR/fc/$DATE/wrfinput_d01_$id
#done

if [ $DATE == $DATE_START ]; then
  #nest down perturbations for inner domains
  if [ $MAX_DOM -gt 1 ]; then
    tid=0
    nt=$((total_ntasks/$wrf_ntasks))
    echo "  Nestdown perturbations for inner domains..."
    for n in `seq 2 $MAX_DOM`; do
      dm=d`expr $n + 100 |cut -c2-`
      parent_dm=d`expr ${PARENT_ID[$n-1]} + 100 |cut -c2-`
      for NE in `seq 1 $NUM_ENS`; do
        id=`expr $NE + 1000 |cut -c2-`
        cd $id
        if [[ ! -d $dm ]]; then mkdir -p $dm; fi
        cd $dm
        ln -fs ../wrfinput_d0? .
        export run_minutes=0
        export start_date=$DATE
        if $FOLLOW_STORM; then
          cp $WORK_DIR/rc/$DATE/ij_parent_start .
        fi
        $SCRIPT_DIR/namelist_wrf.sh ndown $n > namelist.input
        rm -f wrfinput_d0?
        ln -fs $WRF_DIR/run/ndown.exe .
        ln -fs $WORK_DIR/fc/$DATE/wrfinput_${parent_dm}_$id wrfout_d01_`wrf_time_string $DATE`
        ln -fs ../wrfinput_$dm wrfndi_d02
        $SCRIPT_DIR/job_submit.sh $wrf_ntasks $((tid*$wrf_ntasks)) $HOSTPPN ./ndown.exe >& ndown.log &
        tid=$((tid+1))
        if [[ $tid == $nt ]]; then
          tid=0
           wait
        fi
        cd ../..
      done
      wait
      for NE in `seq 1 $NUM_ENS`; do
        id=`expr $NE + 1000 |cut -c2-`
        watch_log $id/$dm/rsl.error.0000 SUCCESS 1 $rundir
        mv $id/$dm/wrfinput_d02 $WORK_DIR/fc/$DATE/wrfinput_${dm}_$id
      done
    done
  fi

  #if using multi-physics ensemble, randomly assign physics options
  if $MULTI_PHYS_ENS; then
    for NE in `seq 1 $NUM_ENS`; do
      id=`expr $NE + 1000 |cut -c2-`
      $SCRIPT_DIR/multi_physics_randraw.sh $WORK_DIR/fc/$DATE/wrfinput_d??_$id >& multi_physics_draw.log
    done
#    $SCRIPT_DIR/multi_physics_draw.sh $NUM_ENS $WORK_DIR/fc/$DATE >& multi_physics_draw.log
  fi
fi

echo complete > stat
