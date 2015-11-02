#!/bin/bash
. $CONFIG_FILE
end_date=`advance_time $start_date $run_minutes`

if [ -f ij_parent_start ]; then
  i_parent_start=(`cat ij_parent_start |head -n1`)
  j_parent_start=(`cat ij_parent_start |tail -n1`)
else
  i_parent_start=(${I_PARENT_START[@]})
  j_parent_start=(${J_PARENT_START[@]})
fi

cat << EOF
&share
 wrf_core = 'ARW',
 max_dom = $MAX_DOM,
 start_year         = `for i in $(seq 1 $MAX_DOM); do printf ${start_date:0:4}, ; done`
 start_month        = `for i in $(seq 1 $MAX_DOM); do printf ${start_date:4:2}, ; done`
 start_day          = `for i in $(seq 1 $MAX_DOM); do printf ${start_date:6:2}, ; done`
 start_hour         = `for i in $(seq 1 $MAX_DOM); do printf ${start_date:8:2}, ; done`
 start_minute       = `for i in $(seq 1 $MAX_DOM); do printf ${start_date:10:2}, ; done`
 start_second       = `for i in $(seq 1 $MAX_DOM); do printf 00, ; done`
 end_year           = `for i in $(seq 1 $MAX_DOM); do printf ${end_date:0:4}, ; done`
 end_month          = `for i in $(seq 1 $MAX_DOM); do printf ${end_date:4:2}, ; done`
 end_day            = `for i in $(seq 1 $MAX_DOM); do printf ${end_date:6:2}, ; done`
 end_hour           = `for i in $(seq 1 $MAX_DOM); do printf ${end_date:8:2}, ; done`
 end_minute         = `for i in $(seq 1 $MAX_DOM); do printf ${end_date:10:2}, ; done`
 end_second         = `for i in $(seq 1 $MAX_DOM); do printf 00, ; done`
 interval_seconds = $((LBC_INTERVAL*60)),
 io_form_geogrid = 2,
 opt_output_from_geogrid_path = './',
 debug_level = 0,
/

&geogrid
 parent_id  = 0,`for i in $(seq 2 $MAX_DOM); do printf ${PARENT_ID[$i-1]}, ; done`
 parent_grid_ratio = 1,`for i in $(seq 2 $MAX_DOM); do printf ${GRID_RATIO[$i-1]}, ; done`
 i_parent_start = 1,`for i in $(seq 2 $MAX_DOM); do printf ${i_parent_start[$i-1]}, ; done`
 j_parent_start = 1,`for i in $(seq 2 $MAX_DOM); do printf ${j_parent_start[$i-1]}, ; done`
 geog_data_res  = `for i in $(seq 1 $MAX_DOM); do printf '30s', ; done`
 e_we       = `for i in $(seq 1 $MAX_DOM); do printf ${E_WE[$i-1]}, ; done`
 e_sn       = `for i in $(seq 1 $MAX_DOM); do printf ${E_SN[$i-1]}, ; done`
 dx = ${DX[0]},
 dy = ${DY[0]},
 map_proj = '$MAP_PROJ',
 ref_lat   = ${REF_LAT[0]},
 ref_lon   = ${REF_LON[0]},
 truelat1  = $TRUELAT1,
 truelat2  = $TRUELAT2,
 stand_lon = $STAND_LON,
 geog_data_path = '$GEOG_DIR',
 opt_geogrid_tbl_path = '$WPS_DIR/geogrid'
/

&ungrib
 out_format = 'WPS'
 prefix = 'FILE',
/

&metgrid
 fg_name = './FILE'
 io_form_metgrid = 2, 
 opt_output_from_metgrid_path = './',
 opt_metgrid_tbl_path         = '$WPS_DIR/metgrid',
/

EOF
