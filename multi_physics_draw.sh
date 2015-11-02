#!/bin/bash
. $CONFIG_FILE
nens=$1
dir=$2

#Define possible choices of physics schemes
#adjacent lines mean paired parameters (cannot mix and match)
BIN_MP_PHYSICS=(2 6 7 10)

BIN_CU_PHYSICS=(0)

BIN_RA_LW_PHYSICS=(3 1)
BIN_RA_SW_PHYSICS=(3 1)

BIN_BL_PBL_PHYSICS=(2)
BIN_SF_SFCLAY_PHYSICS=(2)

BIN_SF_SURFACE_PHYSICS=(1)

i=1
for n0 in `seq 1 ${#BIN_MP_PHYSICS[@]}`; do 
M_MP_PHYSICS=${BIN_MP_PHYSICS[$n0-1]}
for n1 in `seq 1 ${#BIN_CU_PHYSICS[@]}`; do 
M_CU_PHYSICS=${BIN_CU_PHYSICS[$n1-1]}
for n2 in `seq 1 ${#BIN_RA_SW_PHYSICS[@]}`; do 
M_RA_SW_PHYSICS=${BIN_RA_SW_PHYSICS[$n2-1]}
M_RA_LW_PHYSICS=${BIN_RA_LW_PHYSICS[$n2-1]}
for n3 in `seq 1 ${#BIN_BL_PBL_PHYSICS[@]}`; do
M_BL_PBL_PHYSICS=${BIN_BL_PBL_PHYSICS[$n3-1]}
M_SF_SFCLAY_PHYSICS=${BIN_SF_SFCLAY_PHYSICS[$n3-1]}
for n4 in `seq 1 ${#BIN_SF_SURFACE_PHYSICS[@]}`; do
M_SF_SURFACE_PHYSICS=${BIN_SF_SURFACE_PHYSICS[$n4-1]}
  id=`expr $i + 1000 |cut -c2-` 
  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="MP_PHYSICS"' 'attvalue='$M_MP_PHYSICS 'ncfile="'$dir'/wrfinput_d01_'$id'"'
  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="CU_PHYSICS"' 'attvalue='$M_CU_PHYSICS 'ncfile="'$dir'/wrfinput_d01_'$id'"'
  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="RA_SW_PHYSICS"' 'attvalue='$M_RA_SW_PHYSICS 'ncfile="'$dir'/wrfinput_d01_'$id'"'
  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="RA_LW_PHYSICS"' 'attvalue='$M_RA_LW_PHYSICS 'ncfile="'$dir'/wrfinput_d01_'$id'"'
  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="BL_PBL_PHYSICS"' 'attvalue='$M_BL_PBL_PHYSICS 'ncfile="'$dir'/wrfinput_d01_'$id'"'
  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="SF_SFCLAY_PHYSICS"' 'attvalue='$M_SF_SFCLAY_PHYSICS 'ncfile="'$dir'/wrfinput_d01_'$id'"'
  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="SF_SURFACE_PHYSICS"' 'attvalue='$M_SF_SURFACE_PHYSICS 'ncfile="'$dir'/wrfinput_d01_'$id'"'
  if [[ $i -ge $nens ]]; then exit; fi
  i=$((i+1))
done
done
done
done
done
