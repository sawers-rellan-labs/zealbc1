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
# Nextflow env on PERSISTENT storage (the /share/.../env/nextflow one was wiped by cleanup:
# not a valid conda env, no nextflow binary, only Java 8). Create once with Java 17+:
#   conda create -p /rsstu/users/r/rrellan/BZea/envs/nextflow -c conda-forge -c bioconda 'nextflow>=26.04' 'openjdk>=17'
conda activate /rsstu/users/r/rrellan/BZea/envs/nextflow

# cache the per-label envs on persistent storage too, or cleanup wipes them
export NXF_CONDA_CACHEDIR=/rsstu/users/r/rrellan/BZea/envs/nf_cache

nextflow run main.nf -profile slurm -resume
