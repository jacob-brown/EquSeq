#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:mem=1gb

# load R 
echo "Loading R"
module load anaconda3/personal

# load fastqc
echo "Loading fastqc"
module load fastqc/0.11.5

# loop through files and generate report
echo "Generate reports"

# run control file
R --vanilla < $HOME/genomics/code/fastqc_report.R

# end
echo 'Report generation complete'
exit
