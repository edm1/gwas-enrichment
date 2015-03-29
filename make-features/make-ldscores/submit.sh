#!/bin/sh 
#PBS -l nodes=2:ppn=16
#PBS -l walltime=00:02:00:00
#PBS -N n-ldsc
#PBS -m abe
##PBS -q testq

cd $PBS_O_WORKDIR

bash make-commands.sh | time parallel -j 12 --sshloginfile $PBS_NODEFILE --sshdelay 1
