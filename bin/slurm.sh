#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=20g
#SBATCH --time=1-00:00:00
#SBATCH --parsable
#SBATCH -J "ccbr-1326"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output "log/slurm_%j.log"
#SBATCH --output "log/slurm_%j.log"

module load nextflow
NXF_SINGULARITY_CACHEDIR=/data/CCBR_Pipeliner/SIFS

nextflow run main.nf -profile biowulf,slurm -resume
