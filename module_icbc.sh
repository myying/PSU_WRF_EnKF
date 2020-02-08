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

echo "  Preparing IC BC..."
echo running > stat

#0. Calculate nested domain locations (centered over storm) and domain move steps
if $FOLLOW_STORM; then
  echo "    calculate domain moves"
  #Nested domain location i,j: calculate from tcvatils if first cycle, otherwise get from previous cycle outputs
  #if [ $DATE == $DATE_START ]; then
    $SCRIPT_DIR/calc_ij_parent_start.sh $DATE $WORK_DIR/rc/$DATE/ij_parent_start >& follow_storm.log
    watch_file $WORK_DIR/rc/$DATE/ij_parent_start 1 $rundir
    #cp $WORK_DIR/rc/$DATE/ij_parent_start $WORK_DIR/rc/$DATE/ij_parent_start_4dvar
  #else
  #  i_parent_start="1 "
  #  j_parent_start="1 "
  #  for n in `seq 2 $MAX_DOM`; do
  #    dm=d`expr $n + 100 |cut -c2-`
  #    outfile=$WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE`_001
  #    watch_file $outfile 1 $rundir
  #    i_parent_start="$i_parent_start $(ncdump -h $outfile |grep :I_PARENT_START |awk '{print $3}')"
  #    j_parent_start="$j_parent_start $(ncdump -h $outfile |grep :J_PARENT_START |awk '{print $3}')"
  #  done
  #  echo $i_parent_start > $WORK_DIR/rc/$DATE/ij_parent_start
  #  echo $j_parent_start >> $WORK_DIR/rc/$DATE/ij_parent_start
  ##  if $RUN_4DVAR; then
  ##    i_parent_start="1 "
  ##    j_parent_start="1 "
  ##    for n in `seq 2 $MAX_DOM`; do
  ##      dm=d`expr $n + 100 |cut -c2-`
  ##      if $RUN_ENKF; then
  ##        outfile=$WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_$(wrf_time_string `advance_time $DATE $OBS_WIN_MIN`)_mean
  ##      else
  ##        outfile=$WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_$(wrf_time_string `advance_time $DATE $OBS_WIN_MIN`)
  ##      fi
  ##      watch_file $outfile 1 $rundir
  ##      i_parent_start="$i_parent_start $(ncdump -h $outfile |grep :I_PARENT_START |awk '{print $3}')"
  ##      j_parent_start="$j_parent_start $(ncdump -h $outfile |grep :J_PARENT_START |awk '{print $3}')"
  ##    done
  ##    echo $i_parent_start > $WORK_DIR/rc/$DATE/ij_parent_start_4dvar
  ##    echo $j_parent_start >> $WORK_DIR/rc/$DATE/ij_parent_start_4dvar
  ##  fi
  #fi

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

$SCRIPT_DIR/namelist_wps.sh > namelist.wps
#1. geogrid.exe --------------------------------------------------------------------
echo "    running geogrid.exe"
ln -sf $WPS_DIR/geogrid/src/geogrid.exe .
$SCRIPT_DIR/job_submit.sh $wps_ntasks 0 $HOSTPPN ./geogrid.exe >& geogrid.log
watch_log geogrid.log Successful 10 $rundir
mv geo_em.d??.nc $WORK_DIR/rc/$DATE/.
ln -fs $WORK_DIR/rc/$DATE/geo_em.d??.nc .

#2. ungrib.exe --------------------------------------------------------------------
echo "    running ungrib.exe"
if [[ $DATE == $DATE_START ]]; then
  #Link first guess files (FNL, GFS or ECWMF-interim)
  $WPS_DIR/link_grib.csh $FG_DIR/*
  ln -sf $WPS_DIR/ungrib/Variable_Tables/Vtable.GFS Vtable
  ln -fs $WPS_DIR/ungrib/src/ungrib.exe .
  #$SCRIPT_DIR/job_submit.sh $wps_ntasks 0 $HOSTPPN ./ungrib.exe >& ungrib.log
  ./ungrib.exe >& ungrib.log
  watch_log ungrib.log Successful 10 $rundir
else
  rm -f FILE*
  ln -fs $WORK/data/Patricia/icbc/FILE* .
fi

#3. metgrid.exe --------------------------------------------------------------------
echo "    running metgrid.exe"
ln -fs $WPS_DIR/metgrid/METGRID.TBL.ARW METGRID.TBL
ln -fs $WPS_DIR/metgrid/src/metgrid.exe .
$SCRIPT_DIR/job_submit.sh $wps_ntasks 0 $HOSTPPN ./metgrid.exe >& metgrid.log
watch_log metgrid.log Successful 10 $rundir

#4. real.exe ----------------------------------------------------------------------
echo "    running real.exe"
$SCRIPT_DIR/namelist_wrf.sh real > namelist.input
ln -fs $WORK/code/WRFV3.9_serial/main/real.exe .
#$SCRIPT_DIR/job_submit.sh $wps_ntasks 0 $HOSTPPN ./real.exe >& real.log
./real.exe >& rsl.error.0000
watch_log rsl.error.0000 SUCCESS 10 $rundir

###deprecated, new sst update method is replacing SST TSK variables directly
#if [ $SST_UPDATE == 1 ]; then
#  if [ $CYCLE_PERIOD -lt $LBC_INTERVAL ]; then
#    for n in `seq 1 $MAX_DOM`; do
#      dm=d`expr $n + 100 |cut -c2-`
#      ncl $SCRIPT_DIR/util_linint_nc_time.ncl dmin=$CYCLE_PERIOD 'infile="wrflowinp_'$dm'"' >> lowinp.log 2>&1
#      mv tmp.nc $WORK_DIR/rc/$DATE/wrflowinp_$dm
#    done
#  else
#    cp wrflowinp_d?? $WORK_DIR/rc/$DATE
#  fi
#fi
cp wrfinput_d?? $WORK_DIR/rc/$DATE/.
cp wrfbdy_d01 $WORK_DIR/rc/$DATE/.
if [[ $DATE == $DATE_START ]]; then
  cp wrfinput_d?? $WORK_DIR/fc/$DATE/.
  cp wrfbdy_d01 $WORK_DIR/fc/.

  if $RUN_4DVAR; then
    cp $WORK_DIR/fc/wrfbdy_d01 $WORK_DIR/fc/wrfbdy_d01_window
  fi
fi

###get wrfinput for SST update
echo "    preparing sst update"
rm -f FILE*
ln -fs $WORK/data/Patricia/icbc/sst/FILE* .
$SCRIPT_DIR/job_submit.sh $wps_ntasks 0 $HOSTPPN ./metgrid.exe >& metgrid.log
watch_log metgrid.log Successful 10 $rundir
ln -fs $WORK/code/WRFV3.9_serial/main/real.exe .
./real.exe >& rsl.error.0000
watch_log rsl.error.0000 SUCCESS 10 $rundir
for n in `seq 1 $MAX_DOM`; do
  dm=d`expr $n + 100 |cut -c2-`
  cp wrfinput_$dm $WORK_DIR/rc/$DATE/wrfinput_${dm}_sst
done

if $CLEAN; then rm -f *log.???? GRIB* rsl.*; fi
echo complete > stat
