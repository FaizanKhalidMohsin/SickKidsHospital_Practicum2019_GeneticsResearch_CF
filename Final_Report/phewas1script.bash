#!/bin/bash -x

#PBS -l mem=11g
#PBS -l vmem=10g
#PBS -l walltime=23:59:00
#PBS -l nodes=1:ppn=1
#PBS -o /hpf/largeprojects/cfcentre/strug/faizan/job_output/
#PBS -e /hpf/largeprojects/cfcentre/strug/faizan/job_output/
#PBS -d /hpf/largeprojects/cfcentre/strug/faizan/
#PBS -N phewas1run

hostname
date

echo "Working dir is ${PBS_O_WORKDIR}"
cd $PBS_O_WORKDIR


module load R/3.4.0
Rscript phewas1script2.R > phewas1script2.Rout



echo "Done"

date
