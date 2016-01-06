#!/bin/bash
#SBATCH -J CRTM
#SBATCH -n 222
#SBATCH -p normal
#SBATCH -t 8:00:00

cd /work/02135/yingyue/code/CRTM/crtm_wrf
ibrun ./crtm.exe  #> /dev/null
