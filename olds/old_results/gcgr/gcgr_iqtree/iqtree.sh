#!/bin/bash
#
# CompecTA (c) 2018
#
# makeblastdb_trial01 job submission script
#
# TODO:
#   - Set name of the job below changing "NAMD" value.
#   - Set the requested number of nodes (servers) with --nodes parameter.
#   - Set the requested number of tasks (cpu cores) with --ntasks parameter. (Total accross all nodes)
#   - Select the partition (queue) you want to run the job in:
#     - short : For jobs that have maximum run time of 120 mins. Has higher priority.
#     - mid   : For jobs that have maximum run time of 1 days. Lower priority than short.
#     - long  : For jobs that have maximum run time of 7 days. Lower priority than long.
#     - longer: For testing purposes, queue has 15 days limit but only 2 nodes.
#     - cuda  : For CUDA jobs. Solver that can utilize CUDA acceleration can use this queue. 7 days limit.
#   - Set the required time limit for the job with --time parameter.
#     - Acceptable time formats include "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"
#   - Put this script and all the input file under the same directory.
#   - Set the required parameters, input/output file names below.
#   - If you do not want mail please remove the line that has --mail-type and --mail-user. If you do want to get notification emails, set your email address.
#   - Put this script and all the input file under the same directory.
#   - Submit this file using:
#      sbatch slurm_submit.sh
#
# -= Resources =-
#
#SBATCH --job-name=gcgr_iqtree_boostrap_modelfinder
#SBATCH --account=investor
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --qos=mid_investor
#SBATCH --partition=mid_investor
#SBATCH --time=1-0:0
#SBATCH --output=%j-gcgr_iqtree_boostrap_modelfinder_1.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=omerfatihkonar@sabanciuniv.edu

module load iq-tree-2.0.6-gcc-9.2.0-fufn5ou

iqtree2 -s /cta/users/ofkonar/work/results/iqtree/gcgr_aligned.fasta -B 1000 -alrt 1000 -T AUTO --threads-max 10 --prefix gcgr_iqtree_boostrap_modelfinder_1 