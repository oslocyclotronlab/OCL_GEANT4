#!/bin/bash
#
# Laucher batch script file for TACC systems (like Frontera, Stampede2, etc.)
# Si Liu
# July 13, 2020
#
# Simple SLURM script for submitting multiple serial
# jobs (e.g. parametric studies) using a script wrapper
# to launch the jobs.
#
# To use, build the launcher executable and your
# serial application(s) and place them in your WORKDIR
# directory.  Then, edit the LAUNCHER_JOB_FILE to specify 
# each executable per process.
#-------------------------------------------------------
# 
#         <------ Setup Parameters ------>
#
#SBATCH -J resp_postproc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=bigmem
#SBATCH --mem-per-cpu=1GB
#SBATCH -o postproc.%j.out
#SBATCH -e postproc.%j.err
#SBATCH -t 01:30:00
#SBATCH -A nn9464k
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fabio.zeiser@fys.uio.no

#------------------------------------------------------

set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default

# setup working dir (remember to copy back later)
#workdir=$USERWORK/$SLURM_JOB_ID
#mkdir -p $workdir
#cp -r $SUBMITDIR/* $workdir
#cp -r /cluster/home/fabiobz/OCL_GEANT4/* $workdir

#module load launcher
#export LAUNCHER_DIR=/cluster/home/fabiobz/launcher

# USING SLURM; plugins defines SLURM env. vars.
#export LAUNCHER_RMI=SLURM
#export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins

#export LAUNCHER_WORKDIR=$workdir/OscarBuild/
export LAUNCHER_WORKDIR=/cluster/projects/nn9464k/fabio/OSCAR_response_results/641159/data
#export LOGDIR=$workdir/logs_postproc
#mkdir -p $LOGDIR

#export LAUNCHER_JOB_FILE=`pwd`/commands_postproc

#export LAUNCHER_SCHED=block

module load Python/3.8.2-GCCcore-9.3.0 ROOT/6.12.06-intel-2018a-Python-2.7.14 icc/2019.1.144-GCC-8.2.0-2.31.1 CMake/3.13.3-GCCcore-8.2.0
export LD_LIBRARY_PATH=/cluster/home/fabiobz/progs/xerces-c-3.2.3/install/lib:$LD_LIBRARY_PATH
source /cluster/home/fabiobz/progs/geant4.10.06.p02-install/bin/geant4.sh

module list 

#$LAUNCHER_DIR/paramrun
cd $LAUNCHER_WORKDIR
root -l export_hist_short.C
