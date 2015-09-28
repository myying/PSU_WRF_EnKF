#!/bin/bash
. $CONFIG_FILE

for i in `seq 1 $MAX_DOM`; do
  dm=d`expr 100 + $i |cut -c2-`
  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="MP_PHYSICS"' 'attvalue='${MP_PHYSICS[$i-1]} 'ncfile="wrfinput_'$dm'"'
  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="CU_PHYSICS"' 'attvalue='${CU_PHYSICS[$i-1]} 'ncfile="wrfinput_'$dm'"'
  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="RA_SW_PHYSICS"' 'attvalue='${RA_SW_PHYSICS[$i-1]} 'ncfile="wrfinput_'$dm'"'
  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="RA_LW_PHYSICS"' 'attvalue='${RA_LW_PHYSICS[$i-1]} 'ncfile="wrfinput_'$dm'"'
  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="BL_PBL_PHYSICS"' 'attvalue='${BL_PBL_PHYSICS[$i-1]} 'ncfile="wrfinput_'$dm'"'
  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="SF_SFCLAY_PHYSICS"' 'attvalue='${SF_SFCLAY_PHYSICS[$i-1]} 'ncfile="wrfinput_'$dm'"'
  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="SF_SURFACE_PHYSICS"' 'attvalue='${SF_SURFACE_PHYSICS[$i-1]} 'ncfile="wrfinput_'$dm'"'
done
