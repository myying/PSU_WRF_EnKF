#!/bin/bash
. $CONFIG_FILE

rundir=$WORK_DIR/run/$DATE/icbc
if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi

cd $rundir
echo complete > stat ############################

if [[ `cat stat` == "complete" ]]; then exit; fi

#no dependency
if [[ $DATE -gt $DATE_START ]]; then
  wait_for_module ../../$DATE_START/icbc
fi

echo running > stat

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

if [[ $DATE == $DATE_START ]]; then

$SCRIPT_DIR/namelist_wps.sh > namelist.wps
#1. geogrid.exe --------------------------------------------------------------------
  echo "  Running geogrid.exe..."
  ln -sf $WPS_DIR/geogrid/src/geogrid.exe .
#  $SCRIPT_DIR/job_submit.sh $real_ntasks 0 $HOSTPPN ./geogrid.exe >& geogrid.log 
./geogrid.exe >& geogrid.log 
  watch_log geogrid.log Successful 1 $rundir
  mv geo_em.d0?.nc $WORK_DIR/rc/.
ln -fs ../../../rc/geo_em.d0?.nc .

#2. ungrib.exe --------------------------------------------------------------------
echo "  Running ungrib.exe..."
#Link first guess files (FNL, GFS or ECWMF-interim)
$WPS_DIR/link_grib.csh $DATA_DIR/ERA_Interim/*grib
ln -sf $WPS_DIR/ungrib/Variable_Tables/Vtable.ERA-interim.pl Vtable
ln -fs $WPS_DIR/ungrib/src/ungrib.exe .
./ungrib.exe >& ungrib.log
watch_log ungrib.log Successful 2 $rundir

#3. metgrid.exe --------------------------------------------------------------------
echo "  Running metgrid.exe..."
ln -fs $WPS_DIR/metgrid/METGRID.TBL.ARW METGRID.TBL
ln -fs $WPS_DIR/metgrid/src/metgrid.exe .
#$SCRIPT_DIR/job_submit.sh $real_ntasks 0 $HOSTPPN ./metgrid.exe >& metgrid.log
ibrun -np 16 ./metgrid.exe >& metgrid.log
watch_log metgrid.log Successful 2 $rundir
mv met_em* $WORK_DIR/rc/$DATE/.

fi

#4. real.exe ----------------------------------------------------------------------
echo "  Running real.exe..."
export NUM_METGRID_LEVELS=38
$SCRIPT_DIR/namelist_wrf.sh real > namelist.input
ln -fs ../../../rc/$DATE_START/met_em* .
ln -fs $WRF_DIR/main/real.exe .
#$SCRIPT_DIR/job_submit.sh $real_ntasks 0 $HOSTPPN ./real.exe >& real.log
ibrun -np 16 ./real.exe >& real.log
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
#  cp wrfinput_d?? $WORK_DIR/fc/$DATE/.
#  cp wrfbdy_d01 $WORK_DIR/fc/.
#USE WRONG IC/BC FOR OSSE
  cp $SCRATCH/EnKF_OSSE/icbc_FNL/wrfinput_d01_$DATE $WORK_DIR/fc/$DATE/wrfinput_d01
  cp $SCRATCH/EnKF_OSSE/icbc_FNL/wrfbdy_d01_long $WORK_DIR/fc/.

  if $RUN_4DVAR; then
    cp $WORK_DIR/fc/wrfbdy_d01 $WORK_DIR/fc/wrfbdy_d01_window
  fi
fi

if $CLEAN; then rm -f *log.???? met_em* *FILE*; fi
echo complete > stat
