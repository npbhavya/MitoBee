[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![](https://img.shields.io/static/v1?label=CLI&message=Snaketool&color=blueviolet)](https://github.com/beardymcjohnface/Snaketool)
![GitHub last commit (branch)](https://img.shields.io/github/last-commit/npbhavya/MitoMee/main)


# MitoMee

## Snakemake workflow to get mitogenomes from metagenomic data

### Install 

**Source install**
Run the below commands:

    git clone https://github.com/npbhavya/MitoMee.git
    cd MitMee
    mamba create -y -n mitomee python=3.13
    conda activate mitomee
    pip install -e . 

**Downloading the test-files** 

The test files are large, so to download them, use `git-lfs`, Needed if you want to test the installation. 

**Once I have version release, I will upload them to conda and pip as well**

### Running the code

    mitomee run --input test-files/metagenomes --extn fastq.gz --sequencing paired \
         --host_seq test-files/am-dh4.fasta \
         --output output

Input files:
- Input directory with metagenomes
- Reference genome, include only one

Output files: Provide the output folder, contains subdirectories
- PROCESSING: Folder containing intermediate files
- RESULTS: Final results including the mitogenome fasta files from (hopefully) each metagenome sample \
      Also inlcudes the QC reports, to include stats on how many reads were processed, and not


### Build a mitogenome tree

    mitomee tree --input test-files/mitogenomes --extn fasta --host_seq test-files/am-dh4.fasta --output output -k all

Input files:
- Input directory containing all mitochondrial genomes \
  If you have other reference mitogenomes to `output/REPORTS/mitogenomes`. No need to add the reference sequence, it will be included
- Procvide the host reference mitogenome

Output files:
- RESULTS: the tree in nwk format
