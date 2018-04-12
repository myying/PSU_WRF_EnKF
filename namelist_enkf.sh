#!/bin/bash
. $CONFIG_FILE
domain_id=$1
dx=`echo ${DX[$domain_id-1]}/1000 |bc -l`

##switch certain obs type off if OBSINT (obs interval) is set less frequent than CYCLE_PERIOD
offset=`echo "(${DATE:8:2}*60+${DATE:10:2})%${OBSINT_ATOVS:-$CYCLE_PERIOD}" |bc`
if [ $offset != 0 ]; then USE_ATOVS=false; fi

##This if statement swiths the radar rv data off for parent domains
##  the radar data is only assimilated for d03
if [[ $domain_id != 3 ]]; then USE_RADAR_RV=false; fi

cat << EOF
&enkf_parameter
numbers_en   = $NUM_ENS, 
expername    = '$EXP_NAME',  
enkfvar      = 'U         ', 'V         ', 'W         ', 'T         ', 'QVAPOR    ', 'QCLOUD    ', 'QRAIN     ', 'QSNOW     ', 'QICE      ', 'QGRAUP    ', 'QHAIL     ', 'PH        ', 'MU        ', 'PSFC      ', 'P         ', 'PHB       ', 'PB        ', 'MUB       ',
EOF

if [ $minute_off == 0 ] || [ $minute_off == 180 ]; then
  echo "updatevar    = 'U         ', 'V         ', 'W         ', 'T         ', 'QVAPOR    ', 'QCLOUD    ', 'QRAIN     ', 'QSNOW     ', 'QICE      ', 'QGRAUP    ', 'QHAIL     ', 'PH        ', 'MU        ', 'PSFC      ', 'P         ',"
else
  echo "updatevar    = 'U         ', 'V         ', 'W         ', 'T         ', 'QVAPOR    ', 'QCLOUD    ', 'QRAIN     ', 'QSNOW     ', 'QICE      ', 'QGRAUP    ', 'QHAIL     ', 'PH        ', 'MU        ', 'PSFC      ', 'P         ',"
	#echo "updatevar    = 'QCLOUD    ', 'QRAIN     ', 'QSNOW     ', 'QICE      ', 'QGRAUP    ',"
fi

buffer=4 #buffer=0 if update_bc, buffer=spec_bdy_width-1 if bc is fixed as in perfect model case

cat << EOF
update_is    = `echo 1+$buffer |bc`,
update_ie    = `echo ${E_WE[$domain_id-1]}-$buffer |bc`,
update_js    = `echo 1+$buffer |bc`,
update_je    = `echo ${E_SN[$domain_id-1]}-$buffer |bc`,
update_ks    = 1,
update_ke    = `echo ${E_VERT[$domain_id-1]}-1 |bc`,
inflate      = $INFLATION_COEF,
relax_opt    = $RELAX_OPT,
relax_adaptive = .$RELAX_ADAPTIVE.,
mixing       = $RELAXATION_COEF,
random_order = .false.,
print_detail = 0,
/

&parallel
manual_parallel = .true.,
nmcpu  = $NMCPU,
nicpu  = $NICPU,
njcpu  = $NJCPU,
/

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
use_hurricane_PI  = .false.,
hroi_hurricane_PI = 60,
vroi_hurricane_PI = 35,
/

&surface_obs
use_surface      = .$USE_SURFOBS.,
datathin_surface = ${THIN_SURFACE:-0},
hroi_surface     = $(printf %.0f `echo $HROI_SFC/$dx |bc -l`),
vroi_surface     = $VROI,
/

&sounding_obs
use_sounding      = .$USE_SOUNDOBS.,
datathin_sounding = ${THIN_SOUNDING:-0},
datathin_sounding_vert = ${THIN_SOUNDING_VERT:-0},
hroi_sounding     = $(printf %.0f `echo ${HROI_SOUNDING:-$HROI_UPPER}/$dx |bc -l`),
vroi_sounding     = ${VROI_SOUNDING:-$VROI},
/

&profiler_obs
use_profiler      = .$USE_PROFILEROBS.,
datathin_profiler = ${THIN_PROFILER:-0},
datathin_profiler_vert = ${THIN_PROFILER_VERT:-0},
hroi_profiler     = $(printf %.0f `echo $HROI_PROFL/$dx |bc -l`),
vroi_profiler     = $VROI_PROFL,
/

&aircft_obs
use_aircft      = .$USE_AIREPOBS.,
datathin_aircft = ${THIN_AIRCFT:-0},
hroi_aircft     = $(printf %.0f `echo $HROI_UPPER/$dx |bc -l`),
vroi_aircft     = $VROI,
/

&metar_obs
use_metar      = .$USE_METAROBS.,
datathin_metar = ${THIN_METAR:-0},
hroi_metar     = $(printf %.0f `echo $HROI_SFC/$dx |bc -l`),
vroi_metar     = $VROI,
/

&sfcshp_obs
use_sfcshp      = .$USE_SHIPSOBS.,
datathin_sfcshp = ${THIN_SFCSHP:-0},
hroi_sfcshp     = $(printf %.0f `echo $HROI_SFC/$dx |bc -l`),
vroi_sfcshp     = $VROI,
/

&spssmi_obs
use_spssmi      = .$USE_SSMIOBS.,
datathin_spssmi = ${THIN_SPSSMI:-0},
hroi_spssmi     = $(printf %.0f `echo $HROI_UPPER/$dx |bc -l`),
vroi_spssmi     = $VROI,
/

&atovs_obs
use_atovs      = .$USE_ATOVS.,
datathin_atovs = ${THIN_ATOVS:-0},
datathin_atovs_vert = ${THIN_ATOVS_VERT:-0},
hroi_atovs     = $(printf %.0f `echo ${HROI_ATOVS:-$HROI_UPPER}/$dx |bc -l`),
vroi_atovs     = ${VROI_ATOVS:-$VROI},
/

&satwnd_obs
use_satwnd      = .$USE_GEOAMVOBS.,
datathin_satwnd = ${THIN_SATWND:-0},
hroi_satwnd     = $(printf %.0f `echo ${HROI_SATWND:-$HROI_UPPER}/$dx |bc -l`),
vroi_satwnd     = ${VROI_SATWND:-$VROI},
/

&seawind_obs
use_seawind      = .$USE_SEAWIND.,
datathin_seawind = ${THIN_SEAWIND:-0},
hroi_seawind     = $(printf %.0f `echo ${HROI_SEAWIND:-$HROI_UPPER}/$dx |bc -l`),
vroi_seawind     = ${VROI_SEAWIND:-$VROI},
/

&gpspw_obs
use_gpspw      = .$USE_GPSPWOBS.,
datathin_gpspw = ${THIN_GPSPW:-0},
hroi_gpspw     = $(printf %.0f `echo $HROI_SFC/$dx |bc -l`),
vroi_gpspw     = $VROI,
/

&radar_obs
radar_number   = 1,
use_radar_rf   = .$USE_RADAR_RF.,
use_radar_rv   = .$USE_RADAR_RV., 
datathin_radar = $THIN_RADAR,
hroi_radar     = $(printf %.0f `echo $HROI_RADAR/$dx |bc -l`),
vroi_radar     = $VROI_RADAR,
/

&airborne_radar   
use_airborne_rf   = .$USE_AIRBORNE_RF.,
use_airborne_rv   = .$USE_AIRBORNE_RV.,
datathin_airborne = $THIN_RADAR,
hroi_airborne     = $(printf %.0f `echo $HROI_RADAR/$dx |bc -l`),
vroi_airborne     = $VROI_RADAR,
/

&radiance
use_radiance      = .${USE_RADIANCE:-false}.,
datathin_radiance = ${THIN_RADIANCE:-0},
hroi_radiance     = $(printf %.0f `echo ${HROI_RADIANCE:-$HROI_UPPER}/$dx |bc -l`),
vroi_radiance     = ${VROI_RADIANCE:-$VROI},
/
EOF


