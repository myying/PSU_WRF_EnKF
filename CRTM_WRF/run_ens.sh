#!/bin/bash
#SBATCH -J CRTM_Met7_wind
#SBATCH -n 222
#SBATCH -p development
#SBATCH -t 2:00:00
. $WORK/PSU_WRF_EnKF/util.sh 
nens=60
t=201110130000
c=Meteosat7

while [[ $t -lt 201110131800 ]]; do 
echo $t

dir=/scratch/02135/yingyue/EnKF_OSSE/201110120000/$c/run/$t/enkf/d01
outdir=/work/02135/yingyue/DYNAMO/EnKF_OSSE/201110120000/$c
mkdir -p $outdir/$t

for i in `seq 1 $((nens+1))`; do
ibrun /work/02135/yingyue/code/CRTM/crtm_wrf/crtm.io.exe $dir/fort.`expr 80010 + $i` $outdir/$t/fort.`expr 80010 + $i`.bin > /dev/null
ibrun /work/02135/yingyue/code/CRTM/crtm_wrf/crtm.io.exe $dir/fort.`expr 90010 + $i` $outdir/$t/fort.`expr 90010 + $i`.bin > /dev/null
done

  t=`advance_time $t 180`
done
