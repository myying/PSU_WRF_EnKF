#!/bin/bash
#####header for jet######
#PBS -A hfip-psu
#PBS -N gen_be
#PBS -l walltime=0:30:00
#PBS -q debug 
#PBS -l partition=sjet
#PBS -l procs=256
#PBS -j oe
#PBS -d .

######header for stampede######
##SBATCH -J gen_be
##SBATCH -n 256
##SBATCH -p development 
##SBATCH -t 2:00:00

source ~/.bashrc

#load configuration files, functions, parameter
cd $WORK/DA
export CONFIG_FILE=$WORK/DA/config/BAMEX_genbe
. $CONFIG_FILE
. util.sh

####total_ntasks####
if [[ $HOSTTYPE == "stampede" ]]; then
  export total_ntasks=$SLURM_NTASKS
fi
if [[ $HOSTTYPE == "jet" ]]; then
  export total_ntasks=$PBS_NP
fi

#RUN GEN_BE_WRAPPER
#-----------------------------------------------------------------------
# Purpose: Calculates background error statistics for WRF-Var.
#-----------------------------------------------------------------------
#[1] Define job by overriding default environment variables:

export RUN_GEN_BE_STAGE0=true
export RUN_GEN_BE_STAGE1=true
export RUN_GEN_BE_STAGE2=true
export RUN_GEN_BE_STAGE2A=true
export RUN_GEN_BE_STAGE3=true
export RUN_GEN_BE_STAGE4=true
export RUN_GEN_BE_DIAGS=true
export RUN_GEN_BE_DIAGS_READ=true
export RUN_GEN_BE_MULTICOV=true

export WRFVAR_DIR=$WRFDA_DIR
export SCRIPTS_DIR=$SCRIPT_DIR

DATE_START=`advance_time $DATE_START $FORECAST_MINUTES`
DATE_END=`advance_time $DATE_END -$FORECAST_MINUTES`
export START_DATE=${DATE_START:0:10}
export END_DATE=${DATE_END:0:10}
export NUM_LEVELS=${E_VERT[0]}
export BIN_TYPE=5
#export DATA_ON_LEVELS=.true. # "False if fields projected onto modes."

export BE_METHOD=NMC
export FCST_RANGE=$((CYCLE_PERIOD/60))
#Example of changes required for "be_method=ENS":
#export BE_METHOD=ENS
#export NE=2 # 30

export FC_DIR=$WORK_DIR/output	# where wrf forecasts are
export RUN_DIR=$WORK_DIR
export DOMAIN=01
export FCST_RANGE1=$((FORECAST_MINUTES/60))
export FCST_RANGE2=$((CYCLE_PERIOD/60))
export INTERVAL=$((CYCLE_PERIOD/60))
export STRIDE=1
export USE_RFi=true		# use recursive filters?
export NUM_JOBS=$total_ntasks


#[2] Run gen_be:
if ${USE_RFi}; then
   $SCRIPT_DIR/gen_be/gen_be.ksh
else                          # loop over wavelet filter lengths:
   export DO_NORMALIZE=.false.	# normalize before wavelet transform?
   NEW_SUF=
   export RUN_DIR=${RUN_DIR}.
   for L in 7;do
      export WAVELET_NBAND=$L
      for N in C;do          # possible WAVELET_NAME values: B C D V
         export WAVELET_NAME=$N
         if [[ $WAVELET_NAME == B ]];then
            export ISTRT=18      
            export IINC=1
            export IFIN=$ISTRT
         elif [[ $WAVELET_NAME == C ]];then
            export ISTRT=30
            export IINC=6
            export IFIN=30
         elif [[ $WAVELET_NAME == D ]];then
            export ISTRT=6
            export IINC=2
            export IFIN=20
         elif [[ $WAVELET_NAME == V ]];then
            export ISTRT=24      
            export IINC=1
            export IFIN=$ISTRT
         fi
         for I in `seq ${ISTRT} ${IINC} ${IFIN}`; do
            export WAVELET_FILT_LEN=$I
            OLD_SUF=${NEW_SUF}
            NEW_SUF=${WAVELET_NBAND}${WAVELET_NAME}${WAVELET_FILT_LEN}n
            export RUN_DIR=${RUN_DIR%${OLD_SUF}}${NEW_SUF}
            $SCRIPT_DIR/gen_be/gen_be.ksh
         done
      done
   done
fi
