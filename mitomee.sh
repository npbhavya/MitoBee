#!/bin/bash

#SBATCH --job-name mitomee_dev
#SBATCH --output dev-%j.error
#SBATCH --error dev-%j.out
#SBATCH -N 1
#SBATCH --cpus-per-task 16
#SBATCH --partition=cpu
#SBATCH --mem=150G
#SBATCH --time=8:00:00


mitomee run --input test-files/metagenomes --extn fastq.gz --sequencing paired  \
    --host_seq /test-files/am-dh4.fasta \
    --conda-frontend mamba --output output 

#MAGBuild run --input testReads/paired --extn fq.gz --sequencing paired \
#    --host_seq /users/bnalaga1/scratch/reference_db --conda-frontend mamba -k \
#    --output dev_test_run binning

#this is to fix the script from 6_bin_quality onwards
#MAGBuild run --input testReads/paired --extn fq.gz --sequencing paired \
#    --host_seq /users/bnalaga1/scratch/reference_db/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
#    --conda-frontend mamba --output output -k binning

