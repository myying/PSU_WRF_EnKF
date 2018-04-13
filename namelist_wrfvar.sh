#!/bin/bash
. $CONFIG_FILE

if [[ $analysis_type == "RANDOMCV" ]]; then
  VAR_SCALING1=`echo "$VAR_SCALING1*1" |bc -l`
  VAR_SCALING2=`echo "$VAR_SCALING2*1" |bc -l`
  VAR_SCALING3=`echo "$VAR_SCALING3*1" |bc -l`
  VAR_SCALING4=`echo "$VAR_SCALING4*1" |bc -l`
  VAR_SCALING5=`echo "$VAR_SCALING5*1" |bc -l`
fi

cat << EOF
&wrfvar1
var4d=${var4d:-false},
var4d_bin=$((MINUTES_PER_SLOT*60)),
var4d_lbc=$VAR4D_LBC,
multi_inc=${multi_inc:-0},
print_detail_radar=false,
print_detail_xa=false,
print_detail_xb=false,
print_detail_obs=false,
print_detail_grad=false,
print_detail_parallel=false,
/
&wrfvar2
/
EOF

echo "&wrfvar3"
echo ob_format=2,
if [ ! -z $num_fgat_time ]; then echo num_fgat_time=$num_fgat_time,; fi
if [ ! -z $tot_en_size ]; then echo tot_en_size=$tot_en_size,; fi
echo "/"

cat << EOF
&wrfvar4
use_synopobs=$USE_SYNOPOBS,
use_shipsobs=$USE_SHIPSOBS,
use_metarobs=$USE_METAROBS,
use_soundobs=$USE_SOUNDOBS,
use_pilotobs=$USE_PILOTOBS,
use_airepobs=$USE_AIREPOBS,
use_geoamvobs=$USE_GEOAMVOBS,
use_polaramvobs=$USE_POLARAMVOBS,
use_bogusobs=$USE_BOGUSOBS,
use_buoyobs=$USE_BUOYOBS,
use_profilerobs=$USE_PROFILEROBS,
use_satemobs=$USE_SATEMOBS,
use_gpspwobs=$USE_GPSPWOBS,
use_gpsrefobs=$USE_GPSREFOBS,
use_qscatobs=$USE_QSCATOBS,
use_radarobs=$USE_RADAROBS,
/
&wrfvar5
check_max_iv=true,
put_rand_seed=false,
/
&wrfvar6
max_ext_its=${max_ext_its:-1},
ntmax=$NTMAX,
orthonorm_gradient=true,
outer_loop_restart=${outer_loop_restart:-false},
/
&wrfvar7
cv_options=$CV_OPTIONS,
je_factor=$JE_FACTOR,
var_scaling1=$VAR_SCALING1,
var_scaling2=$VAR_SCALING2,
var_scaling3=$VAR_SCALING3,
var_scaling4=$VAR_SCALING4,
var_scaling5=$VAR_SCALING5,
len_scaling1=$LEN_SCALING1,
len_scaling2=$LEN_SCALING2,
len_scaling3=$LEN_SCALING3,
len_scaling4=$LEN_SCALING4,
len_scaling5=$LEN_SCALING5,
as1=`echo "$VAR_SCALING1*0.1" |bc`,`echo "$LEN_SCALING1*0.5" |bc`,`echo "$LEN_SCALING1*1.0" |bc`,
as2=`echo "$VAR_SCALING2*0.1" |bc`,`echo "$LEN_SCALING2*0.5" |bc`,`echo "$LEN_SCALING2*1.0" |bc`,
as3=`echo "$VAR_SCALING3*0.1" |bc`,`echo "$LEN_SCALING3*0.5" |bc`,`echo "$LEN_SCALING3*1.0" |bc`,
as4=`echo "$VAR_SCALING4*0.1" |bc`,`echo "$LEN_SCALING4*0.5" |bc`,`echo "$LEN_SCALING4*1.0" |bc`,
as5=`echo "$VAR_SCALING5*0.1" |bc`,`echo "$LEN_SCALING5*0.5" |bc`,`echo "$LEN_SCALING5*1.0" |bc`,
/
&wrfvar8
/
&wrfvar9
/
&wrfvar10
/
&wrfvar11
cv_options_hum=1,
check_rh=0,
calculate_cg_cost_fn=false,
/
&wrfvar12
/
&wrfvar13
/
&wrfvar14
/
&wrfvar15
/
&wrfvar16
alphacv_method=$ALPHACV_METHOD,
ensdim_alpha=${ensdim_alpha:-0},
alpha_truncation=0,
alpha_corr_type=3,
alpha_corr_scale=900,
alpha_std_dev=1.0,
/
&wrfvar17
analysis_type="$analysis_type",
/
&wrfvar18
analysis_date="`wrf_time_string $time_window_min`.0000"
/
&wrfvar19
/
&wrfvar20
/
&wrfvar21
time_window_min="`wrf_time_string $time_window_min`.0000",
/
&wrfvar22
time_window_max="`wrf_time_string $time_window_max`.0000",
/
&perturbation
jcdfi_use=true,
jcdfi_diag=1,
jcdfi_penalty=1000,
/
&namelist_quilt
/
EOF

