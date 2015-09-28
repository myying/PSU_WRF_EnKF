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

#if CP < LINT, cannot generate wrfinput and wrfbdy from LBC data
#instead, we will fetch wrfbdy from the previous cycle where LBC is available
#and wrfinput will be from the previous cycle wrf run.
#if [[ $LBDATE != $DATE ]]; then echo complete > stat; exit; fi

LINT=360
export start_date=$DATE
if [ $DATE == $DATE_START ]; then
  export run_minutes=`diff_time $DATE_START $DATE_END`
else
  export run_minutes=$run_minutes_forecast
fi

for i in 1; do
if [[ $LBDATE != $DATE ]]; then continue; fi
$SCRIPT_DIR/namelist_wps.sh > namelist.wps

#1. geogrid.exe --------------------------------------------------------------------
if [[ $DATE == $DATE_START ]]; then
  echo "  Running geogrid.exe..."
  ln -sf $WPS_DIR/geogrid/src/geogrid.exe .
  $SCRIPT_DIR/job_submit.sh $real_ntasks 0 $HOSTPPN ./geogrid.exe >& geogrid.log 
#  ./geogrid.exe >& geogrid.log
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
  file="$FG_DIR/${ccyymm:0:4}/fnl_${ccyymm}${dd}_${hh}_00"                                #FNL
  if [ -e $file ]; then 
    gribfile="$gribfile $file"
  fi
  fgdate=`advance_time $fgdate $LINT`
done
$WPS_DIR/link_grib.csh $gribfile
ln -sf $WPS_DIR/ungrib/Variable_Tables/Vtable.GFS Vtable
ln -fs $WPS_DIR/ungrib/src/ungrib.exe .
$SCRIPT_DIR/job_submit.sh 1 0 $HOSTPPN ./ungrib.exe >& ungrib.log
#./ungrib.exe >& ungrib.log
watch_log ungrib.log Successful 2 $rundir

#3. metgrid.exe --------------------------------------------------------------------
echo "  Running metgrid.exe..."
ln -fs $WPS_DIR/metgrid/METGRID.TBL.ARW METGRID.TBL
ln -fs $WPS_DIR/metgrid/src/metgrid.exe .
$SCRIPT_DIR/job_submit.sh $real_ntasks 0 $HOSTPPN ./metgrid.exe >& metgrid.log
#./metgrid.exe >& metgrid.log
watch_log metgrid.log Successful 2 $rundir
#mv met_em* $WORK_DIR/rc/$DATE/.

#4. real.exe ----------------------------------------------------------------------
echo "  Running real.exe..."
$SCRIPT_DIR/namelist_wrf.sh real > namelist.input
#ln -fs ../../../rc/$DATE/met_em* .
ln -fs $WRF_DIR/main/real.exe .
$SCRIPT_DIR/job_submit.sh $real_ntasks 0 $HOSTPPN ./real.exe >& real.log
#./real.exe >& real.log
watch_log rsl.error.0000 SUCCESS 2 $rundir

done


#nestdown from truth simulation
echo "  Running ndown.exe..."
$SCRIPT_DIR/namelist_wrf.sh ndown 2 > namelist.input
wdate=$start_date
while [[ $wdate -le `advance_time $start_date $run_minutes` ]]; do
  ln -fs $WORK/data/DYNAMO/3km_run_9km/wrfout_d01_`wrf_time_string $wdate`
  wdate=`advance_time $wdate $LBC_INTERVAL`
done
#cp -L $WORK_DIR/output/$DATE_START/wrfout_d01_`wrf_time_string $DATE` wrfndi_d02

if [[ $LBDATE != $DATE ]]; then
  cp $WORK_DIR/../1/output/$DATE_START/wrfout_d01_`wrf_time_string $DATE` wrfndi_d02
else
  cp wrfinput_d02 wrfndi_d02
fi

#ncl $SCRIPT_DIR/util_change_nc_att.ncl 'ncfile="wrfndi_d02"' 'attname="I_PARENT_START"' 'attvalue=10'
#ncl $SCRIPT_DIR/util_change_nc_att.ncl 'ncfile="wrfndi_d02"' 'attname="J_PARENT_START"' 'attvalue=10'

ln -fs $WRF_DIR/main/ndown.exe .
$SCRIPT_DIR/job_submit.sh $real_ntasks 0 $HOSTPPN ./ndown.exe >& ndown.log
#./ndown.exe >& ndown.log
watch_log rsl.error.0000 SUCCESS 2 $rundir

#ncl $SCRIPT_DIR/util_change_nc_att.ncl 'ncfile="wrfinput_d02"' 'attname="I_PARENT_START"' 'attvalue=1'
#ncl $SCRIPT_DIR/util_change_nc_att.ncl 'ncfile="wrfinput_d02"' 'attname="J_PARENT_START"' 'attvalue=1'

cp wrfinput_d02 $WORK_DIR/rc/$DATE/wrfinput_d01
cp wrfbdy_d02 $WORK_DIR/rc/$DATE/wrfbdy_d01
if [ $SST_UPDATE == 1 ]; then
  cp wrflowinp_d02 $WORK_DIR/rc/$DATE/wrflowinp_d01
fi

if $CLEAN; then rm -f *log.???? ; fi
echo complete > stat
