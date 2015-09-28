#!/bin/bash
. $CONFIG_FILE

id=$1

for i in 1; do
  dm=d`expr 100 + $i |cut -c2-`

  MP_PHYSICS=$(ncdump -h $WORK_DIR/fc/$DATE_START/wrfinput_${dm}_$id |grep :MP_PHYSICS |awk '{print $3}')
  CU_PHYSICS=$(ncdump -h $WORK_DIR/fc/$DATE_START/wrfinput_${dm}_$id |grep :CU_PHYSICS |awk '{print $3}')
  RA_SW_PHYSICS=$(ncdump -h $WORK_DIR/fc/$DATE_START/wrfinput_${dm}_$id |grep :RA_SW_PHYSICS |awk '{print $3}')
  RA_LW_PHYSICS=$(ncdump -h $WORK_DIR/fc/$DATE_START/wrfinput_${dm}_$id |grep :RA_LW_PHYSICS |awk '{print $3}')
  BL_PBL_PHYSICS=$(ncdump -h $WORK_DIR/fc/$DATE_START/wrfinput_${dm}_$id |grep :BL_PBL_PHYSICS |awk '{print $3}')
  SF_SFCLAY_PHYSICS=$(ncdump -h $WORK_DIR/fc/$DATE_START/wrfinput_${dm}_$id |grep :SF_SFCLAY_PHYSICS |awk '{print $3}')
  SF_SURFACE_PHYSICS=$(ncdump -h $WORK_DIR/fc/$DATE_START/wrfinput_${dm}_$id |grep :SF_SURFACE_PHYSICS |awk '{print $3}')

  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="MP_PHYSICS"' 'attvalue='$MP_PHYSICS 'ncfile="wrfinput_'$dm'"'
  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="CU_PHYSICS"' 'attvalue='$CU_PHYSICS 'ncfile="wrfinput_'$dm'"'
  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="RA_SW_PHYSICS"' 'attvalue='$RA_SW_PHYSICS 'ncfile="wrfinput_'$dm'"'
  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="RA_LW_PHYSICS"' 'attvalue='$RA_LW_PHYSICS 'ncfile="wrfinput_'$dm'"'
  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="BL_PBL_PHYSICS"' 'attvalue='$BL_PBL_PHYSICS 'ncfile="wrfinput_'$dm'"'
  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="SF_SFCLAY_PHYSICS"' 'attvalue='$SF_SFCLAY_PHYSICS 'ncfile="wrfinput_'$dm'"'
  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="SF_SURFACE_PHYSICS"' 'attvalue='$SF_SURFACE_PHYSICS 'ncfile="wrfinput_'$dm'"'
done
