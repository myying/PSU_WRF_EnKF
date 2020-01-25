#!/bin/bash
#HOSTTYPE, HOSTPPN are defined in config file
#JOB_SUBMIT_MODE=1 :qsub run.sh as a single job in queue, once run, the resouces are
#                   scheduled according to work flow choreography (good for crowded queue)
#                   (could be wasting some resources if not choreographed carefully)
#               =2 :./run.sh run directly, each component submitted into queue separately
#                   no choreography needed, -o is discarded, jobs are schedule by host PBS. 
#                   (good for fast queue)

. $CONFIG_FILE

n=$1  # num of tasks job uses
o=$2  # offset location in task list, useful for several jobs to run together
ppn=$3  # proc per node for the job
exe=$4  # executable

###yellowstone/cheyenne NCAR CISL
if [[ $HOSTTYPE == "cheyenne" ]]; then

  if [ $JOB_SUBMIT_MODE == 1 ]; then
#    hosts=($LSB_MCPU_HOSTS)
#    rm -f nodefile_avail
#    for i in `seq 1 $((LSB_MAX_NUM_PROCESSORS/$HOSTPPN))`; do 
#      field1=${hosts[$i*2-2]}
#      field2=${hosts[$i*2-1]}
#      for c in `seq 1 $field2`; do
#        echo $field1 >> nodefile_avail
#      done
#    done
#    cat nodefile_avail |head -n$((o+$n)) |tail -n$n > nodefile 
#    mpirun -np $n -hostfile nodefile $exe
    mpiexec $exe

  elif [ $JOB_SUBMIT_MODE == 2 ]; then
    nodes=`echo "($n+$ppn-1)/$ppn" |bc`
    jobname=`basename $exe |awk -F. '{print $1}'`
    queue="regular"
    wtime="06:00:00"
    cat << EOF > run_$jobname.sh
#!/bin/bash
#PBS -A $HOSTACCOUNT
#PBS -N $jobname
#PBS -l walltime=$wtime
#PBS -q $queue
#PBS -l select=$nodes:ncpus=$ppn:mpiprocs=$ppn
#PBS -j oe
#PBS -o job_run.log

source ~/.bashrc
cd `pwd`
mpiexec $exe >& $jobname.log
EOF
    qsub run_$jobname.sh >& job_submit.log
    #wait for job to finish
    jobid=`cat job_submit.log |cut -c1-15`
    jobstat=1
    until [[ $jobstat == 0 ]]; do
      sleep 10m
      jobstat=`qstat -x |grep $jobid |awk '{if($5==R || $5==Q) print 1; else print 0;}'`
      #########BUGGY
      #########DOESNT WAIT BEFORE FINISH
    done
  fi

fi

###stampede TACC UTEXAS
if [[ $HOSTTYPE == "stampede2" ]]; then

  if [ $JOB_SUBMIT_MODE == 1 ]; then

		export SLURM_NTASKS=$((ppn*$SLURM_NNODES))
		export SLURM_NPROCS=$((ppn*$SLURM_NNODES))
		export SLURM_TACC_CORES=$((ppn*$SLURM_NNODES))
		export SLURM_TASKS_PER_NODE="$ppn(x$SLURM_NNODES)"
		ibrun -n $n -o $o $exe

  elif [ $JOB_SUBMIT_MODE == 2 ]; then
    nodes=`echo "($n+$ppn-1)/$ppn" |bc`
    jobname=`basename $exe |awk -F. '{print $1}'`
    queue="development"
    wtime="00:30:00"
    if [ $jobname == "wrf" ]; then
      queue="normal"
      wtime="00:30:00"
    fi
    if [ $jobname == "enkf" ]; then
      queue="normal"
      wtime="04:00:00"
    fi
    cat << EOF > run_$jobname.sh
#!/bin/bash
#SBATCH -A $HOSTACCOUNT
#SBATCH -J $jobname
#SBATCH -n $n -N $nodes
#SBATCH -p $queue
#SBATCH -t $wtime
#SBATCH -o job_run.log

source ~/.bashrc
cd `pwd`
ibrun $exe >& $jobname.log
EOF
    /bin/sbatch run_$jobname.sh >& job_submit.log
    #wait for job to finish
    jobid=`tail -n1 job_submit.log |cut -c20-`
    jobstat=1
    until [[ $jobstat == 0 ]]; do
      sleep 10m
      jobstat=`/bin/squeue |grep $jobid |awk '{if($5=="PD" || $5=="R") print 1; else print 0;}'`
      if [[ $jobstat == "" ]]; then
        jobstat=0
      fi
    done

  fi
fi

###jet NOAA
if [[ $HOSTTYPE == "jet" ]]; then

  if [ $JOB_SUBMIT_MODE == 1 ]; then
    rm -f nodefile_avail
    for i in `seq 1 $((PBS_NP/$HOSTPPN))`; do 
      cat $PBS_NODEFILE |head -n$((i*$HOSTPPN)) |tail -n$ppn >> nodefile_avail
    done
    cat nodefile_avail |head -n$((o+$n)) |tail -n$n > nodefile 
    mpiexec.mpirun_rsh -np $n -machinefile nodefile OMP_NUM_THREADS=1 $exe
 
  elif [ $JOB_SUBMIT_MODE == 2 ]; then
    nodes=`echo "($n+$ppn-1)/$ppn" |bc`
    jobname=`basename $exe |awk -F. '{print $1}'`
    if [ $jobname == "enkf" ]; then
      queue="debug"
      wtime="0:30:00"
    else
      queue="batch"
      wtime="4:00:00"
#      if [ $jobname == "wrf" ]; then
#        wtime="0:30:00"
#      fi
    fi
    cat << EOF > run_$jobname.sh
#!/bin/bash
#PBS -A hfip-psu
#PBS -N $jobname
#PBS -l walltime=$wtime
#PBS -q $queue
#PBS -l partition=sjet:vjet
#PBS -l nodes=$nodes:ppn=$ppn
#PBS -j oe
#PBS -d .
source ~/.bashrc_yingyue
cd `pwd`
mpiexec -np $n $exe >& $jobname.log
EOF
    qsub run_$jobname.sh >& job_submit.log
    #wait for job to finish
    jobid=`cat job_submit.log |awk -F. '{print $1}'`
    jobstat=1
    until [[ $jobstat == 0 ]]; do
      sleep 10m
      jobstat=`/apps/torque/default/bin/qstat |grep $jobid |awk '{if($5=="R" || $5=="Q") print 1; else print 0;}'`
    done
  fi
fi

###define your own mpiexec here if needed:
#if [[ $HOSTTYPE == "yourHPC" ]]; then
#  
#fi
