#!/bin/bash
#use_for = wrf:     normal use for wrf.exe
#	   wrfvar:  use with da_wrfvar.exe to generate perturbations / run 4dvar; 1 domain only.
#          ndown:   use for ndown.exe, parent-child domains.
#need to define: start_date, run_minutes, write_input (inputout_*)

. $CONFIG_FILE

end_date=`advance_time $start_date $run_minutes`
use_for=$1
idom=$2
thin_factor=${thin_factor:-1}
time_step_ratio=${time_step_ratio:-1}

if [ -f ij_parent_start ]; then
  i_parent_start=(`cat ij_parent_start |head -n1`)
  j_parent_start=(`cat ij_parent_start |tail -n1`)
else
  for n in `seq 1 $MAX_DOM`; do
    i_parent_start[$n-1]=${I_PARENT_START[$n-1]}
    j_parent_start[$n-1]=${J_PARENT_START[$n-1]}
  done
fi

domlist=`seq 1 $MAX_DOM`
if [ $use_for == "wrfvar" ]; then
  MAX_DOM=1
  domlist=$idom
fi
if [ $use_for == "ndown" ]; then
  MAX_DOM=2
  domlist="${PARENT_ID[$idom-1]} $idom"
fi

#=============TIME CONTROL PART=============
echo "&time_control"
cat << EOF
run_minutes        = $run_minutes,
start_year         = $(for i in $domlist; do printf ${start_date:0:4}, ; done)
start_month        = $(for i in $domlist; do printf ${start_date:4:2}, ; done)
start_day          = $(for i in $domlist; do printf ${start_date:6:2}, ; done)
start_hour         = $(for i in $domlist; do printf ${start_date:8:2}, ; done)
start_minute       = $(for i in $domlist; do printf ${start_date:10:2}, ; done)
start_second       = $(for i in $domlist; do printf 00, ; done)
end_year           = $(for i in $domlist; do printf ${end_date:0:4}, ; done)
end_month          = $(for i in $domlist; do printf ${end_date:4:2}, ; done)
end_day            = $(for i in $domlist; do printf ${end_date:6:2}, ; done)
end_hour           = $(for i in $domlist; do printf ${end_date:8:2}, ; done)
end_minute         = $(for i in $domlist; do printf ${end_date:10:2}, ; done)
end_second         = $(for i in $domlist; do printf 00, ; done)
input_from_file    = $(for i in $domlist; do printf .true., ; done)
interval_seconds   = $((LBC_INTERVAL*60)),
history_interval   = $(for i in $domlist; do printf $(min $wrfout_interval ${WRFOUT_INTERVAL[$i-1]}), ; done)
frames_per_outfile = $(for i in $domlist; do printf 1, ; done)
debug_level        = 0,
restart = .false.,
restart_interval = 1440,
EOF

if [ ! -z $inputout_interval ]; then
cat << EOF
input_outname="wrfinput_d<domain>_<date>",
write_input=true,
inputout_interval=$inputout_interval,
inputout_begin_m=$inputout_begin,
inputout_end_m=$inputout_end,
EOF
fi

if [[ $use_for == "ndown" ]]; then
  echo io_form_auxinput2=2,
fi

if [[ $SST_UPDATE == 1 ]]; then
  echo auxinput4_inname="wrflowinp_d<domain>",
  echo auxinput4_interval=$(for i in $domlist; do printf $CYCLE_PERIOD, ; done)
  echo io_form_auxinput4=2,
fi

echo "/"

#=============DOMAIN PART=============
echo "&domains"
echo max_dom    = $MAX_DOM,
if [ $use_for == "wrfvar" ]; then
  echo time_step  = `echo "${time_step:-${TIME_STEP[$idom-1]}}/$time_step_ratio" |bc`,
else
  echo time_step  = `echo "${TIME_STEP[0]}/$time_step_ratio" |bc`,
fi

echo e_we       = $(for i in $domlist; do printf `echo "(${E_WE[$i-1]} -1)/$thin_factor + 1" |bc`, ; done)
echo e_sn       = $(for i in $domlist; do printf `echo "(${E_SN[$i-1]} -1)/$thin_factor + 1" |bc`, ; done)
echo e_vert     = $(for i in $domlist; do printf ${E_VERT[$i-1]}, ; done)
echo dx         = $(for i in $domlist; do printf `echo "${DX[$i-1]} * $thin_factor" |bc`, ; done)
echo dy         = $(for i in $domlist; do printf `echo "${DY[$i-1]} * $thin_factor" |bc`, ; done)

if [[ $use_for == "ndown" ]]; then
cat << EOF
grid_id    = 1,2,
parent_id  = 0,1,
parent_grid_ratio = 1,${GRID_RATIO[$idom-1]},
i_parent_start = 1,${i_parent_start[$idom-1]},
j_parent_start = 1,${j_parent_start[$idom-1]},
EOF
else
cat << EOF
grid_id    = $(for i in $domlist; do printf $i, ; done)
parent_id  = 0,$(for i in $(seq 2 $MAX_DOM); do printf ${PARENT_ID[$i-1]}, ; done)
parent_grid_ratio = 1,$(for i in $(seq 2 $MAX_DOM); do printf ${GRID_RATIO[$i-1]}, ; done)
parent_time_step_ratio = 1,$(for i in $(seq 2 $MAX_DOM); do printf ${TIME_STEP_RATIO[$i-1]}, ; done)
i_parent_start = $(for i in $(seq 1 $MAX_DOM); do printf ${i_parent_start[$i-1]}, ; done)
j_parent_start = $(for i in $(seq 1 $MAX_DOM); do printf ${j_parent_start[$i-1]}, ; done)
EOF
fi

#Preset moves for following TC
if [ -f domain_moves ]; then
  cat domain_moves
fi

if $TWO_WAY_NESTING; then
  echo "feedback=1,"
else
  echo "feedback=0,"
fi

cat << EOF
eta_levels = 1.000000, 0.993042, 0.983211, 0.970571, 0.955211, 0.937239,
0.916793, 0.894028, 0.869104, 0.842199, 0.813512, 0.783267,
0.751688, 0.718973, 0.685325, 0.650946, 0.616034, 0.580793,
0.545438, 0.510175, 0.475184, 0.440634, 0.406689, 0.373503,
0.341211, 0.309918, 0.279723, 0.250727, 0.223030, 0.196731,
0.171918, 0.148674, 0.127087, 0.107244, 0.089261, 0.073324,
0.059519, 0.047598, 0.037267, 0.028358, 0.020713, 0.014185,
0.008638, 0.003947, 0.000000,
smooth_option=0,
num_metgrid_levels=${NUM_METGRID_LEVELS:-38},
p_top_requested=$P_TOP,
num_metgrid_soil_levels=4,
nproc_x=0,
EOF

#echo "eta_levels=1.000,0.99258,0.98275,0.96996,0.95372,0.93357,0.90913,0.87957,0.84531,0.80683,0.76467,0.7194,0.67163,0.62198,0.57108,0.51956,0.46803,0.4203,0.37613,0.33532,0.29764,0.2629,0.23092,0.20152,0.17452,0.14978,0.12714,0.10646,0.08761,0.07045,0.05466,0.03981,0.0258,0.01258,0.000,"
echo "/"

#=============PHYSICS PART=============
#if running multi-physics ensemble read the physics options from wrfinput files
#else get options from config file
echo "&physics"
if $MULTI_PHYS_ENS && [ $use_for == "wrf" ]; then
cat << EOF
mp_physics         = $(for i in $domlist; do printf $(ncdump -h wrfinput_d$(expr $i + 100 |cut -c2-) |grep :MP_PHYSICS |awk '{print $3}'), ; done)
ra_lw_physics      = $(for i in $domlist; do printf $(ncdump -h wrfinput_d$(expr $i + 100 |cut -c2-) |grep :RA_LW_PHYSICS |awk '{print $3}'), ; done)
ra_sw_physics      = $(for i in $domlist; do printf $(ncdump -h wrfinput_d$(expr $i + 100 |cut -c2-) |grep :RA_SW_PHYSICS |awk '{print $3}'), ; done)
sf_sfclay_physics  = $(for i in $domlist; do printf $(ncdump -h wrfinput_d$(expr $i + 100 |cut -c2-) |grep :SF_SFCLAY_PHYSICS |awk '{print $3}'), ; done)
sf_surface_physics = $(for i in $domlist; do printf $(ncdump -h wrfinput_d$(expr $i + 100 |cut -c2-) |grep :SF_SURFACE_PHYSICS |awk '{print $3}'), ; done)
bl_pbl_physics     = $(for i in $domlist; do printf $(ncdump -h wrfinput_d$(expr $i + 100 |cut -c2-) |grep :BL_PBL_PHYSICS |awk '{print $3}'), ; done)
cu_physics         = $(for i in $domlist; do printf $(ncdump -h wrfinput_d$(expr $i + 100 |cut -c2-) |grep :CU_PHYSICS |awk '{print $3}'), ; done)
EOF
else
cat << EOF
mp_physics         = $(for i in $domlist; do printf ${MP_PHYSICS[$i-1]}, ; done)
ra_lw_physics      = $(for i in $domlist; do printf ${RA_LW_PHYSICS[$i-1]}, ; done)
ra_sw_physics      = $(for i in $domlist; do printf ${RA_SW_PHYSICS[$i-1]}, ; done)
sf_sfclay_physics  = $(for i in $domlist; do printf ${SF_SFCLAY_PHYSICS[$i-1]}, ; done)
sf_surface_physics = $(for i in $domlist; do printf ${SF_SURFACE_PHYSICS[$i-1]}, ; done)
bl_pbl_physics     = $(for i in $domlist; do printf ${BL_PBL_PHYSICS[$i-1]}, ; done)
cu_physics         = $(for i in $domlist; do printf ${CU_PHYSICS[$i-1]}, ; done)
EOF
fi

cat << EOF
radt               = $(for i in $domlist; do printf `echo "${RADT[$i-1]} * $thin_factor" |bc`, ; done)
bldt               = $(for i in $domlist; do printf ${BLDT[$i-1]}, ; done)
cudt               = $(for i in $domlist; do printf ${CUDT[$i-1]}, ; done)
EOF

cat << EOF
mp_zero_out        = 2,
sst_update         = $SST_UPDATE,

EOF

#extra physics options here:
cat << EOF
levsiz = 59
paerlev = 29
cam_abs_dim1 = 4
cam_abs_dim2 = 45

 isfflx                              = 1,
 ifsnow                              = 1,
 icloud                              = 1,
 surface_input_source                = 1,

 num_soil_layers                     = 4,
EOF
echo "/"


#=============DYNAMICS PART=============
echo "&dynamics"
cat << EOF
 w_damping                           = 0,
 diff_opt                            = 2,
 km_opt                              = 4,
 diff_6th_opt                        = 0,      0,      0,
 diff_6th_factor                     = 0.12,   0.12,   0.12,
 base_temp                           = 290.
 damp_opt                            = 3,
 zdamp                               = 7000.,  7000.,  5000.,
 dampcoef                            = 0.1,    0.1,    0.2
 khdif                               = 0,      0,      0,
 kvdif                               = 0,      0,      0,
 non_hydrostatic                     = .true., .true., .true.,
 moist_adv_opt                       = 1,      1,      1,
 scalar_adv_opt                      = 1,      1,      1,
EOF
echo "/"

#=============BODY CONTROL PART=============
echo "&bdy_control"
cat << EOF
spec_bdy_width      = 5,
spec_zone           = 1,
relax_zone          = 4,
specified           = $(for i in $domlist; do if [ $i == 1 ]; then printf .true.,; else printf .false.,; fi; done)
nested              = $(for i in $domlist; do if [ $i == 1 ]; then printf .false.,; else printf .true.,; fi; done)
EOF
echo "/"

#=============OTHERS=============
cat << EOF
&noah_mp
/
&fdda
/
&scm
/
&grib2
/
&fire
/
&diags
/
&namelist_quilt
 nio_tasks_per_group = 0,
 nio_groups = 1,
/
&tc
/
&logging
/
&dfi_control
/
EOF
