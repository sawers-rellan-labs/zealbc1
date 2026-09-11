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
# Requires a working nextflow env (>=26.04, JDK>=17). The one below was found broken
# (missing conda-meta + nextflow binary; cause unknown) — recreate if needed:
#   conda create -p /share/maize/frodrig4/conda/env/nextflow -c conda-forge -c bioconda 'nextflow>=26.04' 'openjdk>=17'
conda activate /share/maize/frodrig4/conda/env/nextflow

# cache the per-label envs Nextflow builds from envs/*.yml
export NXF_CONDA_CACHEDIR=/share/maize/frodrig4/conda/nf_cache

nextflow run main.nf -profile slurm -resume
