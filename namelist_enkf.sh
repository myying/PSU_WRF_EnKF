#!/bin/bash
. $CONFIG_FILE
domain_id=$1
current_scale=$2
num_scales=$3
scale_roi_factor=(2.0 1.3 1.0)
thin_factor=(4 2 0)
if [ $num_scales == 1 ]; then
  srf=1
  thin=0
else
  srf=${scale_roi_factor[$current_scale-1]}
  thin=${thin_factor[$current_scale-1]}
fi

dx=`echo ${DX[$domain_id-1]}/1000 |bc -l`

##switch certain obs type off if OBSINT (obs interval) is set less frequent than CYCLE_PERIOD
offset=`echo "(${DATE:8:2}*60+${DATE:10:2})%${OBSINT_ATOVS:-$CYCLE_PERIOD}" |bc`
if [[ $offset != 0 ]]; then USE_ATOVS=false; fi

#if [ $DATE != $LBDATE ]; then USE_HURRICANE=false; fi
##This if statement swiths the radar rv data off for parent domains
##  the radar data is only assimilated for d03
#if [[ $domain_id != 3 ]]; then USE_RADAR_RV=false; fi

use_airborne_rv=${USE_AIRBORNE_RV[$domain_id-1]}

buffer=4 #buffer=0 if update_bc, buffer=spec_bdy_width-1 if bc is fixed as in perfect model case

cat << EOF
&enkf_parameter
numbers_en   = $NUM_ENS,
expername    = '$EXP_NAME',
enkfvar      = `for i in ${UPDATE_VAR[*]}; do printf "\'%-10s\', " $i; done`
updatevar    = `for i in ${UPDATE_VAR[*]}; do printf "\'%-10s\', " $i; done`
update_is    = `echo 1+$buffer |bc`,
update_ie    = `echo ${E_WE[$domain_id-1]}-$buffer |bc`,
update_js    = `echo 1+$buffer |bc`,
update_je    = `echo ${E_SN[$domain_id-1]}-$buffer |bc`,
update_ks    = 1,
update_ke    = `echo ${E_VERT[$domain_id-1]}-1 |bc`,
inflate      = ${INFLATION_COEF:-1.0},
relax_opt    = ${RELAX_OPT:-0},
relax_adaptive = .${RELAX_ADAPTIVE:-false}.,
mixing       = ${RELAXATION_COEF[$domain_id-1]:-0.0},
random_order = .false.,
print_detail = 0,
/

&parallel
manual_parallel = .true.,
nmcpu  = $NMCPU,
nicpu  = $NICPU,
njcpu  = $NJCPU,
/
EOF

if [ $num_scales == 1 ]; then
  cat << EOF
&multiscale
num_scales = 1,
krange = 0.00,
current_scale = 1,
run_alignment = .false.,
/
EOF
else
  cat << EOF
&multiscale
num_scales = $num_scales,
krange = $(for k in ${KRANGE[*]}; do printf '%5.2f, ' `echo $k/${DX[0]}*${DX[0]} |bc -l`; done)
current_scale = ${current_scale:-1},
run_alignment = .$RUN_ALIGNMENT.,
/
EOF
fi
cat << EOF
&osse
use_ideal_obs    = .false.,
gridobs_is   = 20,
gridobs_ie   = `echo ${E_WE[$domain_id-1]}-20 |bc`,
gridobs_js   = 20,
gridobs_je   = `echo ${E_SN[$domain_id-1]}-20 |bc`,
gridobs_ks   = 1,
gridobs_ke   = `echo ${E_VERT[$domain_id-1]}-1 |bc`,
gridobs_int_x= 40,
gridobs_int_k= 1,
use_simulated= .false.,
/

&hurricane_PI
use_hurricane_PI  = .$USE_HURRICANE.,
hroi_hurricane_PI = $(printf %.0f `echo $HROI_HURRICANE*$srf/$dx |bc -l`),
vroi_hurricane_PI = $VROI,
/

&surface_obs
use_surface      = .$USE_SURFOBS.,
datathin_surface = ${THIN_SURFACE:-0},
hroi_surface     = $(printf %.0f `echo $HROI_SFC*$srf/$dx |bc -l`),
vroi_surface     = $VROI,
/

&sounding_obs
use_sounding      = .$USE_SOUNDOBS.,
datathin_sounding = ${THIN_SOUNDING:-0},
datathin_sounding_vert = ${THIN_SOUNDING_VERT:-0},
hroi_sounding     = $(printf %.0f `echo ${HROI_SOUNDING:-$HROI_UPPER}*$srf/$dx |bc -l`),
vroi_sounding     = ${VROI_SOUNDING:-$VROI},
/

&profiler_obs
use_profiler      = .$USE_PROFILEROBS.,
datathin_profiler = ${THIN_PROFILER:-0},
datathin_profiler_vert = ${THIN_PROFILER_VERT:-0},
hroi_profiler     = $(printf %.0f `echo ${HROI_PROFILER:-$HROI_UPPER}*$srf/$dx |bc -l`),
vroi_profiler     = ${VROI_PROFILER:-$VROI},
/

&aircft_obs
use_aircft      = .$USE_AIREPOBS.,
datathin_aircft = ${THIN_AIRCFT:-0},
hroi_aircft     = $(printf %.0f `echo $HROI_UPPER*$srf/$dx |bc -l`),
vroi_aircft     = $VROI,
/

&metar_obs
use_metar      = .$USE_METAROBS.,
datathin_metar = ${THIN_METAR:-0},
hroi_metar     = $(printf %.0f `echo $HROI_SFC*$srf/$dx |bc -l`),
vroi_metar     = $VROI,
/

&sfcshp_obs
use_sfcshp      = .$USE_SHIPSOBS.,
datathin_sfcshp = ${THIN_SFCSHP:-0},
hroi_sfcshp     = $(printf %.0f `echo $HROI_SFC*$srf/$dx |bc -l`),
vroi_sfcshp     = $VROI,
/

&spssmi_obs
use_spssmi      = .$USE_SSMIOBS.,
datathin_spssmi = ${THIN_SPSSMI:-0},
hroi_spssmi     = $(printf %.0f `echo $HROI_UPPER*$srf/$dx |bc -l`),
vroi_spssmi     = $VROI,
/

&atovs_obs
use_atovs      = .$USE_ATOVS.,
datathin_atovs = ${THIN_ATOVS:-0},
datathin_atovs_vert = ${THIN_ATOVS_VERT:-0},
hroi_atovs     = $(printf %.0f `echo ${HROI_ATOVS:-$HROI_UPPER}*$srf/$dx |bc -l`),
vroi_atovs     = ${VROI_ATOVS:-$VROI},
/

&satwnd_obs
use_satwnd      = .$USE_GEOAMVOBS.,
datathin_satwnd = `expr $THIN_SATWND + $thin`,
hroi_satwnd     = $(printf %.0f `echo ${HROI_SATWND:-$HROI_UPPER}*$srf/$dx |bc -l`),
vroi_satwnd     = ${VROI_SATWND:-$VROI},
/

&seawind_obs
use_seawind      = .$USE_SEAWIND.,
datathin_seawind = ${THIN_SEAWIND:-0},
hroi_seawind     = $(printf %.0f `echo ${HROI_SEAWIND:-$HROI_UPPER}*$srf/$dx |bc -l`),
vroi_seawind     = ${VROI_SEAWIND:-$VROI},
/

&gpspw_obs
use_gpspw      = .$USE_GPSPWOBS.,
datathin_gpspw = ${THIN_GPSPW:-0},
hroi_gpspw     = $(printf %.0f `echo $HROI_SFC*$srf/$dx |bc -l`),
vroi_gpspw     = $VROI,
/

&radar_obs
radar_number   = 1,
use_radar_rf   = .$USE_RADAR_RF.,
use_radar_rv   = .$USE_RADAR_RV.,
datathin_radar = $THIN_RADAR,
hroi_radar     = $(printf %.0f `echo $HROI_RADAR*$srf/$dx |bc -l`),
vroi_radar     = $VROI_RADAR,
/

&airborne_radar
use_airborne_rf   = .$USE_AIRBORNE_RF.,
use_airborne_rv   = .$use_airborne_rv.,
datathin_airborne = `expr $THIN_RADAR + $thin`,
hroi_airborne     = $(printf %.0f `echo $HROI_RADAR*$srf/$dx |bc -l`),
vroi_airborne     = $VROI_RADAR,
/

&radiance
use_radiance      = .${USE_RADIANCE[$domain_id-1]:-false}.,
datathin_radiance = `expr ${THIN_RADIANCE[$domain_id-1]} + $thin`,
hroi_radiance     = $(printf %.0f `echo ${HROI_RADIANCE:-$HROI_UPPER}*$srf/$dx |bc -l`),
vroi_radiance     = ${VROI_RADIANCE:-99},
/
EOF

