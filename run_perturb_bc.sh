#!/bin/bash
#BSUB -P UPSU0001
#BSUB -J perturb_bc
#BSUB -W 1:00
#BSUB -q small
#BSUB -n 16
#BSUB -R "span[ptile=4]"
#BSUB -o log.%J.out
#BSUB -e log.%J.err

. $WORK/PSU_WRF_EnKF/util.sh
DATE_START=201110120000
DATE_END=201111070000

#tid=0
#nt=$total_ntasks
for id in {003..060}; do
  echo member $id
  if [[ ! -d $id ]]; then mkdir $id; fi

  cd $id

  ln -fs $WORK/code/WRF_BC_v2.1_alltime/update_wrf_bc.exe .
  ln -fs ../../../../fc/wrfbdy_d01 wrfbdy_d01_real
  cp -L wrfbdy_d01_real wrfbdy_d01_update

  dd=`diff_time $DATE_START $DATE_END`

  for n_1 in `seq 1 $((dd/360+1))`; do
  cat > parame.in << EOF
&control_param
 wrf_3dvar_output_file = 'wrfinput_d01_update'
 wrf_bdy_file          = 'wrfbdy_d01_update'
 wrf_bdy_file_real     = 'wrfbdy_d01_real'
 wrf_input_from_si     = 'wrfinput_d01_real'
 wrf_input_from_si_randmean = 'random_mean'
 wrf_3dvar_random_draw = 'random_draw'
 cycling = .true.
 low_bdy_only = .false. 
 perturb_bdy = .true.
 n_1 = $n_1
/
EOF

  ln -fs ../../../../rc/$DATE_START/wrfinput_d01_`advance_time $DATE_START $((n_1*360-360))` wrfinput_d01_real
  ln -fs wrfinput_d01_real wrfinput_d01_update
  randnum=`expr $((RANDOM%99+1)) + 1000 |cut -c2-`
  echo $n_1 $randnum >> rand
  ln -fs ../../../../fc/$DATE_START/wrfinput_d01_$randnum random_draw
  ln -fs ../../../../fc/$DATE_START/wrfinput_d01 random_mean

#  $SCRIPT_DIR/job_submit.sh 1 $tid $HOSTPPN ./update_wrf_bc.exe >& update_wrf_bc.log &
  ./update_wrf_bc.exe >& update_wrf_bc.log 
  done

  mv wrfbdy_d01_update ../../../../fc/wrfbdy_d01_$id
  mv rand ../../../../fc/rand_$id

  cd .. 
done
