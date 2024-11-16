#!/bin/bash
#SBATCH -p long
#SBATCH --job-name=mouse_wt_long_dox_timecourse
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=krle4401@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=6gb
#SBATCH --time=10:00:00
#SBATCH --output=nextflow.out
#SBATCH --error=nextflow.err

pwd; hostname; date
echo "Here we go You've requested $SLURM_CPUS_ON_NODE core."

module load singularity/3.1.1

nextflow run nf-core/rnaseq -r 3.14.0 \
-resume \
-profile singularity \
--input samplesheet.csv \
--outdir /scratch/Shares/rinn/Kit/class_RNASeq/output_dir \
--reads /scratch/Shares/rinn/ML/mouse_wt_long_timecourse/*{_R1,_R2}.fastq.gz \
--fasta /scratch/Shares/rinn/ML/DESeq2/GENOMES/M25/GRCm38.p6.genome.fa \
--gtf /scratch/Shares/rinn/ML/DESeq2/GENOMES/M25/gencode.vM25.annotation.gtf \
--pseudo_aligner salmon \
--gencode \
--email krle4401@colorado.edu \
-c nextflow.config

date