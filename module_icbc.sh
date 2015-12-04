#!/bin/bash
. $CONFIG_FILE

rundir=$WORK_DIR/run/$DATE/icbc
if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi

cd $rundir
if [[ `cat stat` == "complete" ]]; then exit; fi

#no dependency
if [[ $DATE -gt $DATE_START ]]; then
  wait_for_module ../../$DATE_START/icbc
fi

echo running > stat

#0. Calculate nested domain locations (centered over storm) and domain move steps
if $FOLLOW_STORM; then
  #Nested domain location i,j: calculate from tcvatils if first cycle, otherwise get from previous cycle outputs
  if [ $DATE == $DATE_START ]; then
    $SCRIPT_DIR/calc_ij_parent_start.sh $DATE $WORK_DIR/rc/$DATE/ij_parent_start >& follow_storm.log
    watch_file $WORK_DIR/rc/$DATE/ij_parent_start 1 $rundir
    cp $WORK_DIR/rc/$DATE/ij_parent_start $WORK_DIR/rc/$DATE/ij_parent_start_4dvar
  else
    if $RUN_ENKF; then
      i_parent_start="1 "
      j_parent_start="1 "
      for n in `seq 2 $MAX_DOM`; do
        dm=d`expr $n + 100 |cut -c2-`
        outfile=$WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE`_001
        watch_file $outfile 1 $rundir
        i_parent_start="$i_parent_start $(ncdump -h $outfile |grep :I_PARENT_START |awk '{print $3}')"
        j_parent_start="$j_parent_start $(ncdump -h $outfile |grep :J_PARENT_START |awk '{print $3}')"
      done
      echo $i_parent_start > $WORK_DIR/rc/$DATE/ij_parent_start
      echo $j_parent_start >> $WORK_DIR/rc/$DATE/ij_parent_start
    fi
    if $RUN_4DVAR; then
      i_parent_start="1 "
      j_parent_start="1 "
      for n in `seq 2 $MAX_DOM`; do
        dm=d`expr $n + 100 |cut -c2-`
        if $RUN_ENKF; then
          outfile=$WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_$(wrf_time_string `advance_time $DATE $OBS_WIN_MIN`)_mean
        else
          outfile=$WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_$(wrf_time_string `advance_time $DATE $OBS_WIN_MIN`)
        fi
        watch_file $outfile 1 $rundir
        i_parent_start="$i_parent_start $(ncdump -h $outfile |grep :I_PARENT_START |awk '{print $3}')"
        j_parent_start="$j_parent_start $(ncdump -h $outfile |grep :J_PARENT_START |awk '{print $3}')"
      done
      echo $i_parent_start > $WORK_DIR/rc/$DATE/ij_parent_start_4dvar
      echo $j_parent_start >> $WORK_DIR/rc/$DATE/ij_parent_start_4dvar
    fi
  fi

  #Domain move steps
  $SCRIPT_DIR/calc_domain_moves.sh $DATE $NEXTDATE $WORK_DIR/rc/$DATE/domain_moves >& follow_storm.log
  watch_file $WORK_DIR/rc/$DATE/domain_moves 1 $rundir
  if $RUN_4DVAR; then
    if [ $DATE == $DATE_START ]; then
      $SCRIPT_DIR/calc_domain_moves.sh $DATE `advance_time $NEXTDATE $OBS_WIN_MIN` $WORK_DIR/rc/$DATE/domain_moves_4dvar >& follow_storm.log
    else
      $SCRIPT_DIR/calc_domain_moves.sh `advance_time $DATE $OBS_WIN_MIN` `advance_time $NEXTDATE $OBS_WIN_MIN` $WORK_DIR/rc/$DATE/domain_moves_4dvar >& follow_storm.log
    fi
    watch_file $WORK_DIR/rc/$DATE/domain_moves_4dvar 1 $rundir
  fi

  ln -fs $WORK_DIR/rc/$DATE/ij_parent_start .
fi

#if CP < LBC_INTERVAL, cannot generate wrfinput and wrfbdy from LBC data
#instead, we will fetch wrfbdy from the previous cycle where LBC is available
#and wrfinput will be from the previous cycle wrf run.
if [[ $LBDATE != $DATE ]]; then echo complete > stat; exit; fi

export start_date=$DATE
if [ $DATE == $DATE_START ]; then
  export run_minutes=`diff_time $DATE_START $DATE_END`
else
  export run_minutes=$((LBC_INTERVAL*2))
fi

#export run_minutes=`max $run_minutes $run_minutes_forecast`

$SCRIPT_DIR/namelist_wps.sh > namelist.wps

#1. geogrid.exe --------------------------------------------------------------------
if [[ $DATE == $DATE_START ]]; then
  echo "  Running geogrid.exe..."
  ln -sf $WPS_DIR/geogrid/src/geogrid.exe .
  $SCRIPT_DIR/job_submit.sh $real_ntasks 0 $HOSTPPN ./geogrid.exe >& geogrid.log 
  watch_log geogrid.log Successful 1 $rundir
  mv geo_em.d0?.nc $WORK_DIR/rc/.
fi
ln -fs ../../../rc/geo_em.d0?.nc .

#2. ungrib.exe --------------------------------------------------------------------
echo "  Running ungrib.exe..."
#Link first guess files (FNL, GFS or ECWMF-interim)
fgdate=$start_date
gribfile=""
while [[ $fgdate -le `advance_time $start_date $run_minutes` ]]; do
  ccyymm=`echo $fgdate |cut -c1-6`
  dd=`echo $fgdate |cut -c7-8`
  hh=`echo $fgdate |cut -c9-10`
#  file="$FG_DIR/gfs.$ccyymm$dd$hh/`date -u -d $ccyymm$dd' '$hh':00' +%y%j%H`000000" #GFS
  file="$FG_DIR/fnl_${ccyymm}${dd}_${hh}_00"                                #FNL
  if [ -e $file ]; then 
    gribfile="$gribfile $file"
  fi
  fgdate=`advance_time $fgdate $LBC_INTERVAL`
done
$WPS_DIR/link_grib.csh $gribfile
ln -sf $WPS_DIR/ungrib/Variable_Tables/Vtable.GFS Vtable
ln -fs $WPS_DIR/ungrib/src/ungrib.exe .
./ungrib.exe >& ungrib.log
watch_log ungrib.log Successful 2 $rundir

#3. metgrid.exe --------------------------------------------------------------------
echo "  Running metgrid.exe..."
ln -fs $WPS_DIR/metgrid/METGRID.TBL.ARW METGRID.TBL
ln -fs $WPS_DIR/metgrid/src/metgrid.exe .
$SCRIPT_DIR/job_submit.sh $real_ntasks 0 $HOSTPPN ./metgrid.exe >& metgrid.log
watch_log metgrid.log Successful 2 $rundir
mv met_em* $WORK_DIR/rc/$DATE/.

#4. real.exe ----------------------------------------------------------------------
echo "  Running real.exe..."
$SCRIPT_DIR/namelist_wrf.sh real > namelist.input
ln -fs ../../../rc/$DATE/met_em* .
ln -fs $WRF_DIR/main/real.exe .
$SCRIPT_DIR/job_submit.sh $real_ntasks 0 $HOSTPPN ./real.exe >& real.log
watch_log rsl.error.0000 SUCCESS 2 $rundir
if [ $SST_UPDATE == 1 ]; then
  if [ $CYCLE_PERIOD -lt $LBC_INTERVAL ]; then
    for n in `seq 1 $MAX_DOM`; do
      dm=d`expr $n + 100 |cut -c2-`
      ncl $SCRIPT_DIR/util_linint_nc_time.ncl dmin=$CYCLE_PERIOD 'infile="wrflowinp_'$dm'"' >> lowinp.log 2>&1
      mv tmp.nc $WORK_DIR/rc/$DATE/wrflowinp_$dm
    done
  else
    cp wrflowinp_d?? $WORK_DIR/rc/$DATE
  fi
fi
cp wrfinput_d?? $WORK_DIR/rc/$DATE/.
cp wrfbdy_d01 $WORK_DIR/rc/$DATE/.
if [[ $DATE == $DATE_START ]]; then
  cp wrfinput_d?? $WORK_DIR/fc/$DATE/.
  cp wrfbdy_d01 $WORK_DIR/fc/.
  if $RUN_4DVAR; then
    cp $WORK_DIR/fc/wrfbdy_d01 $WORK_DIR/fc/wrfbdy_d01_window
  fi
fi

if $CLEAN; then rm -f *log.???? ; fi
echo complete > stat
