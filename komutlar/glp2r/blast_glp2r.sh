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
#SBATCH --job-name=blast_glp2r
#SBATCH --account=investor
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=short_investor
#SBATCH --partition=short_investor
#SBATCH --time=120
#SBATCH --output=%j-blast_glp2r.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=omerfatihkonar@sabanciuniv.edu
module load blast-plus-2.9.0-gcc-9.2.0-k73hoo5
blastp -query /cta/users/ofkonar/work/fasta/glp2r_human.fasta -db /cta/users/ofkonar/work/uniprot_database/uniprot_db -out /cta/users/ofkonar/work/results/glp2r/glp2r_blast_results_csv.csv -outfmt "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sacc stitle"