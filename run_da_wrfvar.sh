#!/bin/bash
#BSUB -P UPSU0001
#BSUB -J da_wrfvar
#BSUB -W 2:00
#BSUB -q regular 
#BSUB -n 64
#BSUB -R "span[ptile=4]"
#BSUB -o log.%J.out
#BSUB -e log.%J.err

source /glade/u/apps/opt/lmod/4.2.1/init/bash
source ~/.bashrc_yy

for i in {073..100}; do 
cd /glade/scratch/fzhang/yy/DYNAMO/EnKF_OSSE/AtovsAmvMet7Cygnss/run/201110120000/perturb_ic/$i
mpirun.lsf ./da_wrfvar.exe >& da_wrfvar.log
mv wrfvar_output ../../../../fc/201110120000/wrfinput_d01_$i
done
