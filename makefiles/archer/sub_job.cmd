#!/bin/bash --login
#PBS -N sub_job
#PBS -l select=1
#PBS -l walltime=01:00:00
#PBS -A e01-WP2_Sheff14
#PBS -o sub_job.out
#PBS -e sub_job.err

# Make sure any symbolic links are resolved to absolute path
export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)               
  
# Change to the direcotry that the job was submitted from
cd $PBS_O_WORKDIR

# export PROCS=`qstat -u $USER | sed '1,5d' | awk '{print $6*24}'`

# Set the number of threads to 1
#   This prevents any system libraries from automatically 
#   using threading.

export OMP_NUM_THREADS=1

# Launch the parallel job
aprun -n 16 ./inc3d_archer
