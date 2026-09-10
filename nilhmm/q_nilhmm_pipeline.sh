#!/bin/bash
#SBATCH --job-name=nilhmm
#SBATCH --account=maize_cpu
#SBATCH --partition=compute
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# nilhmm Phase-1 head job. The Nextflow head runs here (1 cpu) and submits every
# process as its own Slurm job via the 'slurm' profile.
# Run from this directory:  sbatch q_nilhmm_pipeline.sh

source ~/.bashrc
conda activate /share/maize/frodrig4/conda/env/nextflow

nextflow run main.nf -profile slurm -resume
