#!/bin/bash
. $CONFIG_FILE

#Define possible choices of physics schemes
#adjacent lines mean paired parameters (cannot mix and match)
BIN_MP_PHYSICS=(2 6 8 7 16)

BIN_CU_PHYSICS=(0)

BIN_RA_SW_PHYSICS=(4)
BIN_RA_LW_PHYSICS=(5)

BIN_BL_PBL_PHYSICS=(1) # 2 99)
BIN_SF_SFCLAY_PHYSICS=(1) # 2 1)

BIN_SF_SURFACE_PHYSICS=(1)

#draw a random physics option and modify wrfinput file accordingly
n=$((RANDOM%${#BIN_MP_PHYSICS[@]}))
export M_MP_PHYSICS=${BIN_MP_PHYSICS[$n]}
n=$((RANDOM%${#BIN_CU_PHYSICS[@]}))
export M_CU_PHYSICS=${BIN_CU_PHYSICS[$n]}
n=$((RANDOM%${#BIN_RA_SW_PHYSICS[@]}))
export M_RA_SW_PHYSICS=${BIN_RA_SW_PHYSICS[$n]}
export M_RA_LW_PHYSICS=${BIN_RA_LW_PHYSICS[$n]}
n=$((RANDOM%${#BIN_BL_PBL_PHYSICS[@]}))
export M_BL_PBL_PHYSICS=${BIN_BL_PBL_PHYSICS[$n]}
export M_SF_SFCLAY_PHYSICS=${BIN_SF_SFCLAY_PHYSICS[$n]}
n=$((RANDOM%${#BIN_SF_SURFACE_PHYSICS[@]}))
export M_SF_SURFACE_PHYSICS=${BIN_SF_SURFACE_PHYSICS[$n]}

for wrffile in $@; do
  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="MP_PHYSICS"' 'attvalue='$M_MP_PHYSICS 'ncfile="'$wrffile'"'
  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="CU_PHYSICS"' 'attvalue='$M_CU_PHYSICS 'ncfile="'$wrffile'"'
  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="RA_SW_PHYSICS"' 'attvalue='$M_RA_SW_PHYSICS 'ncfile="'$wrffile'"'
  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="RA_LW_PHYSICS"' 'attvalue='$M_RA_LW_PHYSICS 'ncfile="'$wrffile'"'
  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="BL_PBL_PHYSICS"' 'attvalue='$M_BL_PBL_PHYSICS 'ncfile="'$wrffile'"'
  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="SF_SFCLAY_PHYSICS"' 'attvalue='$M_SF_SFCLAY_PHYSICS 'ncfile="'$wrffile'"'
  ncl $SCRIPT_DIR/util_change_nc_att.ncl 'attname="SF_SURFACE_PHYSICS"' 'attvalue='$M_SF_SURFACE_PHYSICS 'ncfile="'$wrffile'"'
done
