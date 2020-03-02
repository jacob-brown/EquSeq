#!/bin/bash
#PBS -l walltime=00:10:00
#PBS -l select=1:ncpus=1:mem=1gb


module load fastqc/0.11.5
#fastqc $EPHEMERAL/read1_trim.fq.gz

fastqc -d . -o $EPHEMERAL $EPHEMERAL/read1_trim.fq.gz

