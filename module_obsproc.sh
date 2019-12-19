#!/bin/bash
. $CONFIG_FILE

rundir=$WORK_DIR/run/$DATE/obsproc

if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi

cd $rundir
if [[ `cat stat` == "complete" ]]; then exit; fi
echo running > stat
echo "  Running ObsProc..."

ln -fs $WRFDA_DIR/var/obsproc/obsproc.exe .
ln -fs $WRFDA_DIR/var/obsproc/obserr.txt .
echo > obs.raw

##### include NCAR_LITTLE_R (3-hourly) #####
if $INCLUDE_LITTLE_R; then
  rm -f datelist
  time_lag=1
  obs_interval=1
  for offset in `seq $((OBS_WIN_MIN/60-$time_lag)) $obs_interval $((OBS_WIN_MAX/60+$time_lag))`; do
    obsdate=`advance_time $DATE $((offset*60))`
    hh=`echo $obsdate |cut -c9-10`
    inc=`echo $hh%$obs_interval*60 |bc`
    if [[ $inc -lt $((obs_interval*60/2)) ]]; then
      obsdate=`advance_time $obsdate -$inc`
    else
      obsdate=`advance_time $obsdate $((obs_interval*60-inc))`
    fi
    echo $obsdate >> datelist
  done
  for d in `cat datelist |sort |uniq`; do

    #NCAR_LITTLE_R
    if [ -f $DATA_DIR/ncar_littler/${d:0:6}/obs.${d:0:10}.gz ]; then
      echo "ncar_littler"
      cp $DATA_DIR/ncar_littler/${d:0:6}/obs.${d:0:10}.gz .
      gunzip obs.${d:0:10}.gz
      cat obs.${d:0:10} >> obs.raw
      rm obs.${d:0:10}
    fi

    if [ -f $DATA_DIR/dropsonde/littler/${d:0:10} ]; then
      echo "dropsonde"
      cat $DATA_DIR/dropsonde/littler/${d:0:10} >> obs.raw
    fi
  done
fi

##### include MADIS data (hourly) data ######
if $INCLUDE_MADIS; then
  rm -f datelist
  time_lag=1
  obs_interval=1
  for offset in `seq $((OBS_WIN_MIN/60-$time_lag)) $obs_interval $((OBS_WIN_MAX/60+$time_lag))`; do
    obsdate=`advance_time $DATE $((offset*60))`
    hh=`echo $obsdate |cut -c9-10`
    inc=`echo $hh%$obs_interval*60 |bc`
    if [[ $inc -lt $((obs_interval*60/2)) ]]; then
      obsdate=`advance_time $obsdate -$inc`
    else
      obsdate=`advance_time $obsdate $((obs_interval*60-inc))`
    fi
    echo $obsdate >> datelist
  done
  for d in `cat datelist |sort |uniq`; do
    cat $DATA_DIR/madis/`echo $d |cut -c1-6`/madis_`echo $d |cut -c1-10`_littler >> obs.raw
  done
fi

##### include BUFR ADP Surface and Upperair (6-hourly) data #####
if $INCLUDE_BUFR; then
  bufr_sfc_decode_dir=$WORK/code/BUFR/bufr_decode_ADPsfc_littler/exe
  bufr_upa_decode_dir=$WORK/code/BUFR/bufr_decode_ADPupa_littler/exe
  bufr_dir=$WORK/data/bufr
  rm -f datelist
  time_lag=6
  obs_interval=6
  for offset in `seq $((OBS_WIN_MIN/60-$time_lag)) $obs_interval $((OBS_WIN_MAX/60+$time_lag))`; do
    obsdate=`advance_time $DATE $((offset*60))`
    hh=`echo $obsdate |cut -c9-10`
    inc=`echo $hh%$obs_interval*60 |bc`
    if [[ $inc -lt $((obs_interval*60/2)) ]]; then
      obsdate=`advance_time $obsdate -$inc`
    else
      obsdate=`advance_time $obsdate $((obs_interval*60-inc))`
    fi
    echo $obsdate >> datelist
  done
  for d in `cat datelist |sort |uniq`; do
    ccyymmdd=`echo $d |cut -c1-8`
    hh=`echo $d |cut -c9-10`
    $bufr_sfc_decode_dir/bufr_sfc2ob.x $bufr_dir/sfcobs.$ccyymmdd/gdas.adpsfc.t${hh}z.$ccyymmdd.bufr $ccyymmdd$hh >> decode_bufr.log 2>&1
    $bufr_sfc_decode_dir/bufr_ship2ob.x $bufr_dir/sfcobs.$ccyymmdd/gdas.sfcshp.t${hh}z.$ccyymmdd.bufr $ccyymmdd$hh >> decode_bufr.log 2>&1
    $bufr_upa_decode_dir/bufr_upa2ob.x $bufr_dir/upaobs.$ccyymmdd/gdas.adpupa.t${hh}z.$ccyymmdd.bufr $ccyymmdd$hh >> decode_bufr.log 2>&1
    $bufr_upa_decode_dir/bufr_aircar2ob.x $bufr_dir/upaobs.$ccyymmdd/gdas.aircar.t${hh}z.$ccyymmdd.bufr $ccyymmdd$hh >> decode_bufr.log 2>&1
    $bufr_upa_decode_dir/bufr_craft2ob.x $bufr_dir/upaobs.$ccyymmdd/gdas.aircft.t${hh}z.$ccyymmdd.bufr $ccyymmdd$hh >> decode_bufr.log 2>&1
    $bufr_upa_decode_dir/bufr_sat2ob.x $bufr_dir/upaobs.$ccyymmdd/gdas.satwnd.t${hh}z.$ccyymmdd.bufr $ccyymmdd$hh >> decode_bufr.log 2>&1
    rm -f files.txt
    for type in Airca Aircraft Satob Ship Surface Upper; do
      echo $type$ccyymmdd$hh.obs >> files.txt
    done
    $bufr_sfc_decode_dir/runob2lit_imd_obs.x files.txt $ccyymmdd$hh >> decode_bufr.log 2>&1
    cat OBS:$ccyymmdd$hh >> obs.raw
    cat SURFACE_OBS:$ccyymmdd$hh >> obs.raw
  done
  rm -f OBS:* SURFACE_OBS:* *.obs
fi

#####  START OBSPROC #####
for var_type in 3DVAR 4DVAR; do
  case $var_type in
    3DVAR)
      if ! $RUN_ENKF; then continue; fi
    ;;
    4DVAR)
      if ! $RUN_4DVAR; then continue; fi
    ;;
  esac
  echo > obsproc.log
  export use_for=$var_type
  $SCRIPT_DIR/namelist_obsproc.sh > namelist.obsproc
  ./obsproc.exe >& obsproc.log

  watch_log obsproc.log 99999 1 $rundir

  mv obs_gts_*.$var_type $WORK_DIR/obs/$DATE/.
done

###squeeze in HPI obs_gts file
mv $WORK_DIR/obs/$DATE/obs_gts_`wrf_time_string $DATE`.3DVAR tmpfile
total=`head -n1 tmpfile |cut -c8-14`
echo "`head -n1 tmpfile |cut -c1-7``printf '%7i' $((total+1))``head -n1 tmpfile |cut -c15-33`" > tmpfile1
head -n21 tmpfile |tail -n20 >> tmpfile1
tail -n3 $WORK/data/Patricia/HPI/HPI_obs_gts_`wrf_time_string $DATE`.3DVAR >> tmpfile1
cat tmpfile |awk 'NR>21{print}' >> tmpfile1
mv tmpfile1 $WORK_DIR/obs/$DATE/obs_gts_`wrf_time_string $DATE`.3DVAR

#if $CLEAN; then rm obs.raw; fi

echo complete > stat

