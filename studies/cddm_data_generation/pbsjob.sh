#!/bin/bash
# Torque directives (#PBS) must always be at the start of a job script!

#PBS -N two_plasctics
#PBS -q mse
#PBS -l nodes=1:ppn=8
#PBS -o out.$PBS_JOBID
#PBS -e err.$PBS_JOBID
#PBS -M J.Yi@tudelft.nl


# Make sure I'm the only one that can read my output
umask 0077


PBS_ARRAYID=$(echo "${PBS_JOBID}" | sed -n 's/.*\[\([^]]*\)\].*/\1/p')

JOB_ID=$(echo "${PBS_JOBID}" | sed 's/\[[^][]*\]//g')
echo "Directory is $PBS_JOBID"

module load use.own
module load miniconda3
module load abaqus/2021
cd $PBS_O_WORKDIR

# Here is where the application is started on the node
# activating my conda environment:
source activate MLenv3


 # limiting number of threads
 OMP_NUM_THREADS=12
 export OMP_NUM_THREADS=12

if ! [ -n "${PBS_ARRAYID+1}" ]; then
  PBS_ARRAYID=0
fi

#Executing my python program
python main.py ++hpc.jobid=${PBS_ARRAYID}