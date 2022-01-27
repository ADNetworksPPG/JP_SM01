#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=04:00:00
#$ -l mem_free=4G
#$ -t 1-8

dataDir=/wynton/group/gladstone/biocore/projects/StephanieMiller-SM-01-behavior-motifs-sept-2021/data/
scriptDir=/wynton/group/gladstone/biocore/projects/StephanieMiller-SM-01-behavior-motifs-sept-2021/scripts/
resultsDir=/wynton/group/gladstone/biocore/projects/StephanieMiller-SM-01-behavior-motifs-sept-2021/results/
container=/wynton/group/gladstone/biocore/projects/StephanieMiller-SM-01-behavior-motifs-sept-2021/scripts/markovchain_latest.sif
export SINGULARITY_BINDPATH="$dataDir,$scriptDir, $resultsDir"

time_steps=(0 1 6 9 12 16 25 50 125 250)
singularity exec   $container Rscript test_markov_property_control_mice.R "${time_steps[$SGE_TASK_ID]}"

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID" 
