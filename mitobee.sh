#!/bin/bash

#SBATCH --job-name mitomee_dev
#SBATCH --output dev-%j.error
#SBATCH --error dev-%j.out
#SBATCH -N 1
#SBATCH --cpus-per-task 16
#SBATCH --partition=cpu
#SBATCH --mem=150G
#SBATCH --time=8:00:00


mitomee run --input test-files/metagenomes --extn fastq.gz --sequencing paired \
     --host_seq test-files/am-dh4.fasta \
     --conda-frontend mamba --output output

#build tree with host seq, assembled mitogenomes. 
#Also add other reference mitogenomes to  output/REPORTS/mitogenomes. No need to add the reference sequence, it will be included
#Then run the subcommand to build tree 

mitomee tree --input test-files/mitogenomes --extn fasta --host_seq test-files/am-dh4.fasta --output output -k all
