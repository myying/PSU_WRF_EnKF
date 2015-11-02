#!/bin/bash
. $CONFIG_FILE
rundir=$WORK_DIR/run/$DATE/4dvar
if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi

cd $rundir
if [[ `cat stat` == "complete" ]]; then exit; fi

#Check dependency
wait_for_module ../icbc ../wrf_window1
if $RUN_ENVAR; then
  wait_for_module ../wrf_ens_window1
fi
if $RUN_ENKF; then
  wait_for_module ../../$PREVDATE/wrf_ens
fi

echo running > stat
echo "  Running 4DVar..."

export analysis_type="3D-VAR"
export time_window_min=`advance_time $DATE $OBS_WIN_MIN`
export time_window_max=`advance_time $DATE $OBS_WIN_MAX`

domlist=`seq 1 $MAX_DOM`

fgat_num=0
for i in `seq $OBS_WIN_MIN $MINUTES_PER_SLOT $OBS_WIN_MAX`; do fgat_num=$((fgat_num+1)); done

#Prepare input files: ep fg fg02 wrfbdy be.dat obs-------------------------------------------
for n in $domlist; do
  dm=d`expr $n + 100 |cut -c2-`
  if [[ ! -d $dm ]]; then mkdir -p $dm; echo running > $dm/stage0; fi
  if [[ `cat $dm/stage0` == "complete" ]]; then continue; fi
  cd $dm

  #Calculate ensemble perturbation (ep)
  if $RUN_ENKF; then
    echo "    calculating ep for domain $dm"
    if $RUN_ENVAR; then
      ep_offset=`seq $OBS_WIN_MIN $MINUTES_PER_SLOT $OBS_WIN_MAX`
    else
      ep_offset=$OBS_WIN_MIN
    fi
    j=1
    for i in $ep_offset; do
      epdate=`advance_time $DATE $i`
      if [[ ! -d ens_$epdate ]]; then mkdir -p ens_$epdate; fi
      for NE in `seq 1 $NUM_ENS`; do
        id=`expr $NE + 1000 |cut -c2-`
        if $RUN_ENVAR; then
          infile=$WORK_DIR/fc/$DATE/wrfinput_${dm}_`wrf_time_string $epdate`_window_$id
        else
          infile=$WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $epdate`_$id
        fi
        outfile=$rundir/$dm/ens_$epdate/wrfinput_${dm}.e$id
        if $MULTI_INC; then
           $MULTI_INC_DIR/da_decimation.exe -i $infile -o $outfile -thin $DECIMATION_FACTOR >& da_decimation.log
        else
           ln -fs $infile $outfile
        fi
      done

      if $RUN_ENVAR; then
        epdir=ep`expr 100 + $j |cut -c2-`
      else 
        epdir=ep
      fi
      if [[ ! -d $epdir ]]; then mkdir -p $epdir; fi
      cd $epdir
      if [ $ALPHACV_METHOD == 1 ]; then
        epexe=$WRFDA_DIR/var/da/gen_be_ep1.exe
      else
        epexe=$WRFDA_DIR/var/da/gen_be_ep2.exe
      fi
      $epexe ${epdate:0:10} $NUM_ENS ../ens_$epdate wrfinput_$dm >& ../ens_ep_$epdate.log &
      cd ..
      j=$((j+1))
    done
    wait
  fi


  #Link background error covariance (be.dat)
  if [[ $CV_OPTIONS == 3 ]]; then
    ln -fs $WRFDA_DIR/var/run/be.dat.cv3 be.dat
  else
    if $MULTI_INC; then
      ln -fs $BE_DIR/be.dat.lores be.dat
    else
      ln -fs $BE_DIR/be.dat .
    fi
  fi

  #Link Observations
  j=1
  for i in `seq $OBS_WIN_MIN $MINUTES_PER_SLOT $OBS_WIN_MAX`; do
    obsdate=`advance_time $DATE $i`
#    ln -fs $DATA_DIR/madis/${DATE:0:10}/obs_gts_`wrf_time_string $obsdate`.4DVAR ob`expr 100 + $j |cut -c2-`.ascii
    ln -fs $OBS_DIR/$DATE/obs_gts_`wrf_time_string $obsdate`.4DVAR.$OBS_TYPE ob`expr 100 + $j |cut -c2-`.ascii
    j=$((j+1))
  done

  #Prepare first guess and BC
  if $RUN_ENKF; then
    fg_file=$WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $time_window_min`_mean
    fg02_file=$WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $time_window_max`_mean
  else
    fg_file=$WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $time_window_min`
    fg02_file=$WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $time_window_max`
  fi
  if $MULTI_INC; then
    ln -fs $fg_file fg_hires
    ln -fs $fg02_file fg02_hires
    ln -fs $WORK_DIR/fc/wrfbdy_d01_window wrfbdy_d01_orig
    $MULTI_INC_DIR/da_decimation.exe -i fg_hires -o fg_thin -thin $DECIMATION_FACTOR >& da_decimation.log
    $MULTI_INC_DIR/da_decimation.exe -i fg02_hires -o fg02_thin -thin $DECIMATION_FACTOR >& da_decimation.log
    $MULTI_INC_DIR/da_bdy.exe -fg fg_hires -fg02 fg02_hires -bdy wrfbdy_d01_orig -o wrfbdy_hires >& da_bdy.log
    $MULTI_INC_DIR/da_bdy.exe -fg fg_thin -fg02 fg02_thin -bdy wrfbdy_d01_orig -o wrfbdy_thin >& da_bdy.log
  else
    ln -fs $fg_file fg
    ln -fs $fg02_file fg02
    ln -fs $fg_file wrfinput_d01
    ln -fs $WORK_DIR/fc/wrfbdy_d01_window wrfbdy_d01_orig
    $MULTI_INC_DIR/da_bdy.exe -fg fg -fg02 fg02 -bdy wrfbdy_d01_orig -o wrfbdy_d01 >& da_bdy.log
  fi

  #Other necessary files
  ln -fs $WRFDA_DIR/run/* .
  for i in `/bin/ls *DBL`; do   #use *DATA_DBL instead of *DATA 
    rm `echo $i |rev |cut -c5- |rev`
    mv $i `echo $i |rev |cut -c5- |rev`
  done

  ln -fs $WRFDA_DIR/var/da/da_wrfvar.exe

  echo complete > stage0
  cd ..
done

#Run wrfvar -----------------------------------------------------------------------------
if $MULTI_INC; then

  #Multi-incremental 4DVar
  echo "    stage1. first outer loop with full resolution"
  tid=0
  nn=$((($var4d_ntasks+$var4d_ppn-$var4d_ntasks%$var4d_ppn)/$var4d_ppn))
  nt=$(($total_ntasks/$HOSTPPN/$nn))
  for n in $domlist; do
    dm=d`expr 100 + $n |cut -c2-`
    if [ ! -f $dm/stage1 ]; then echo running > $dm/stage1; fi
    if [[ `cat $dm/stage1` == "complete" ]]; then continue; fi
    cd $dm

    export multi_inc=1
    export num_fgat_time=$fgat_num
    if $RUN_ENVAR; then
      export var4d=false
      export tot_en_size=$NUM_ENS
    else
      export var4d=true
    fi
    $SCRIPT_DIR/namelist_wrfvar.sh > namelist.input.stage1
    export start_date=$time_window_min
    export run_minutes=`echo "$OBS_WIN_MAX - $OBS_WIN_MIN" |bc`
    export time_step=${TIME_STEP[$n-1]}
    $SCRIPT_DIR/namelist_wrf.sh wrfvar $n >> namelist.input.stage1

    ln -fs namelist.input.stage1 namelist.input
    ln -fs fg_hires fg
    ln -fs fg02_hires fg02
    ln -fs fg_hires wrfinput_d01
    ln -fs wrfbdy_hires wrfbdy_d01

    $SCRIPT_DIR/job_submit.sh ${var4d_ntasks} $((tid*$var4d_ntasks)) $var4d_ppn ./da_wrfvar.exe >& da_wrfvar.log &
    tid=$((tid+1))
    if [[ $tid == $nt ]]; then
      tid=0
      wait
    fi
    cd ..
  done
  wait
  for n in $domlist; do
    dm=d`expr $n + 100 |cut -c2-`
    if [[ `cat $dm/stage1` == "complete" ]]; then continue; fi
    cd $dm
    watch_log rsl.error.0000 successfully 1 $rundir
    for f in `ls rsl*`; do mv $f stage1_outer_$f; done
    echo complete > stage1
    cd ..
  done

  #stage2. additional inner/outer loops
  for it in `seq 1 $MAX_EXT_ITS`; do

    echo "    stage2: inner loop with low resolution (iteration $it)"
    tid=0
    nn=$((($var4d_ntasks+$var4d_ppn-$var4d_ntasks%$var4d_ppn)/$var4d_ppn))
    nt=$(($total_ntasks/$HOSTPPN/$nn))
    for n in $domlist; do
      dm=d`expr 100 + $n |cut -c2-`
      if [ ! -f $dm/stage2_inner$it ]; then echo running > $dm/stage2_inner$it; fi
      if [[ `cat $dm/stage2_inner$it` == "complete" ]]; then continue; fi 
      cd $dm

      if [ $it == 1 ]; then
        export multi_inc=2
        export num_fgat_time=$fgat_num
        if $RUN_ENVAR; then
          export tot_en_size=$NUM_ENS
          export var4d=false
        else
          export var4d=true
        fi
        if $RUN_ENKF; then
          export ensdim_alpha=$NUM_ENS
        else
          export ensdim_alpha=0
        fi
        export outer_loop_restart=true
        $SCRIPT_DIR/namelist_wrfvar.sh > namelist.input.stage2
        export start_date=$time_window_min
        export run_minutes=`echo "$OBS_WIN_MAX - $OBS_WIN_MIN" |bc`
        export thin_factor=$DECIMATION_FACTOR
        export time_step=`echo "${TIME_STEP[$n-1]} * $DECIMATION_FACTOR" |bc`
        until [ `echo "($MINUTES_PER_SLOT*60) % $time_step" |bc` == 0 ]; do export time_step=$((time_step-1)); done
        $SCRIPT_DIR/namelist_wrf.sh wrfvar $n >> namelist.input.stage2
        ln -fs namelist.input.stage2 namelist.input
        ln -fs fg_thin fg
        ln -fs fg02_thin fg02
        ln -fs fg_thin wrfinput_d01
        ln -fs wrfbdy_thin wrfbdy_d01
      else
        ln -fs namelist.input.stage2 namelist.input
        ln -fs fg_thin_new fg
        ln -fs fg02_thin_new fg02
        ln -fs fg_thin_new wrfinput_d01
        ln -fs wrfbdy_thin_new wrfbdy_d01
      fi
      $SCRIPT_DIR/job_submit.sh ${var4d_ntasks} $((tid*$var4d_ntasks)) $var4d_ppn ./da_wrfvar.exe >& da_wrfvar.log &
      tid=$((tid+1))
      if [[ $tid == $nt ]]; then
        tid=0
        wait
      fi
      cd ..
    done
    wait
    for n in $domlist; do
      dm=d`expr $n + 100 |cut -c2-`
      if [[ `cat $dm/stage2_inner$it` == "complete" ]]; then continue; fi
      cd $dm
      watch_log rsl.error.0000 successfully 1 $rundir
      mv wrfvar_output fg_thin_new
      if $VAR4D_LBC; then
        mv ana02 fg02_thin_new
      else
        cp fg02_thin fg02_thin_new
      fi
      for f in `ls rsl.*`; do mv $f stage$((it+1))_inner_$f; done
      mv cost_fn cost_fn_$((it+1))
      mv grad_fn grad_fn_$((it+1))

      #increment regridding
      $MULTI_INC_DIR/da_bilin.exe -fg_lores fg_thin -an_lores fg_thin_new -fg_hires fg_hires -ns $DECIMATION_FACTOR -o fg_hires_new >& da_bilin.log
      $MULTI_INC_DIR/da_bilin.exe -fg_lores fg02_thin -an_lores fg02_thin_new -fg_hires fg02_hires -ns $DECIMATION_FACTOR -o fg02_hires_new >& da_bilin.log

      $MULTI_INC_DIR/da_bdy.exe -fg fg_hires_new -fg02 fg02_hires_new -bdy wrfbdy_d01_orig -o wrfbdy_hires >& da_bdy.log
      $MULTI_INC_DIR/da_bdy.exe -fg fg_thin_new -fg02 fg02_thin_new -bdy wrfbdy_d01_orig -o wrfbdy_thin >& da_bdy.log

      echo complete > stage2_inner$it
      cd ..
    done

    if [ $it == $MAX_EXT_ITS ]; then continue; fi
    echo "    stage2: outer loop with full resolution (iteration $it)"
    tid=0
    nn=$((($var4d_ntasks+$var4d_ppn-$var4d_ntasks%$var4d_ppn)/$var4d_ppn))
    nt=$(($total_ntasks/$HOSTPPN/$nn))
    for n in $domlist; do
      dm=d`expr 100 + $n |cut -c2-`
      if [ ! -f $dm/stage2_outer$it ]; then echo running > $dm/stage2_outer$it; fi
      if [[ `cat $dm/stage2_outer$it` == "complete" ]]; then continue; fi
      cd $dm

      rm -f gts*
      ln -fs namelist.input.stage1 namelist.input
      ln -fs fg_hires_new fg
      ln -fs fg02_hires_new fg02
      ln -fs fg_hires_new wrfinput_d01
      ln -fs wrfbdy_hires wrfbdy_d01
      $SCRIPT_DIR/job_submit.sh ${var4d_ntasks} $((tid*$var4d_ntasks)) $var4d_ppn ./da_wrfvar.exe >& da_wrfvar.log &
      tid=$((tid+1))
      if [[ $tid == $nt ]]; then
        tid=0
        wait
      fi
      cd ..
    done
    wait
    for n in $domlist; do
      dm=d`expr $n + 100 |cut -c2-`
      if [[ `cat $dm/stage2_outer$it` == "complete" ]]; then continue; fi
      cd $dm 
      watch_log rsl.error.0000 successfully 1 $rundir
      for f in `ls rsl*`; do mv $f stage$((it+1))_outer_$f; done
      echo complete > stage2_outer$it
      cd ..
    done

  done #it

  #Final analysis
  for n in `seq 1 $MAX_DOM`; do
    dm=d`expr $n + 100 |cut -c2-`
    cd $dm
    mv fg_hires_new wrfvar_output
    cd ..
  done

else  #full resolution 4DVar

  tid=0
  nn=$((($var4d_ntasks+$var4d_ppn-$var4d_ntasks%$var4d_ppn)/$var4d_ppn))
  nt=$(($total_ntasks/$HOSTPPN/$nn))
  for n in $domlist; do
    dm=d`expr 100 + $n |cut -c2-`
    cd $dm
   
    export num_fgat_time=$fgat_num
    if $RUN_ENKF; then
      export ensdim_alpha=$NUM_ENS
    else
      export ensdim_alpha=0
    fi
    if $RUN_ENVAR; then
      export var4d=false
      export tot_en_size=$NUM_ENS
    else
      export var4d=true
    fi
    $SCRIPT_DIR/namelist_wrfvar.sh > namelist.input
    export start_date=$time_window_min
    export run_minutes=`echo "$OBS_WIN_MAX - $OBS_WIN_MIN" |bc`
    export time_step=${TIME_STEP[$n-1]}
    $SCRIPT_DIR/namelist_wrf.sh wrfvar $n >> namelist.input
    $SCRIPT_DIR/job_submit.sh $var4d_ntasks $((tid*$var4d_ntasks)) $var4d_ppn ./da_wrfvar.exe >& da_wrfvar.log &
    tid=$((tid+1))
    if [[ $tid == $nt ]]; then
      tid=0
      wait
    fi
    cd ..
  done
  wait
  for n in $domlist; do
    dm=d`expr 100 + $n |cut -c2-`
    watch_log $dm/rsl.error.0000 successfully 1 $rundir
  done

fi  #Multi-inc

for n in $domlist; do
  dm=d`expr 100 + $n |cut -c2-`
  watch_file $dm/wrfvar_output 1 $rundir
  mv $dm/wrfvar_output $WORK_DIR/fc/$DATE/wrfinput_${dm}_`wrf_time_string $time_window_min`
done

echo complete > stat
