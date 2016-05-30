#!/bin/bash -l
#
#$ -N networkReliability
#$ -q low.q
#$ -pe threaded 1
#$ -w e
#$ -l h_vmem=5G
#$ -l h_rt=12:00:00
#$ -l s_rt=11:59:50
#$ -v SCENARIO_INDEX
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -R y
trap "echo recieved SIGUSR1;" SIGUSR1;
R CMD BATCH --no-save --no-restore -- jobScript.R jobScript.Rout.$JOB_ID
