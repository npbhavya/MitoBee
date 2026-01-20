"""
Separate Snakefile for building phylogenetic trees
- download and add any reference mitogenomes to "REPORTS/mitogenome/"
- build a phylogenetic tree using FastTree
"""
import yaml
import os
import glob

"""Parse config"""
configfile: os.path.join(workflow.basedir, "..", "config", "config.yaml")

"""
Input files and directories
"""
input_dir = config['args']['input']
extn=config['args']['extn']

# Pattern: all files with this extension
pattern = os.path.join(input_dir, f"*.{extn}")
file_paths = sorted(glob.glob(pattern))
# Extract sample names (strip extension only)
sample_names = [re.sub(rf"\.{re.escape(ext)}$", '', os.path.basename(fp)) for fp in file_paths]

print(f"Samples found: {sample_names}")


"""
Defining the targets dictionary
"""
dir_hostcleaned = os.path.join(dir_out, 'PROCESSING' ,'2_host_cleaned')
dir_reports = os.path.join(dir_out, 'REPORTS')

targets ={'host':[]}
#Definiting the targets
targets['host'].append(os.path.join(dir_reports, "mitogenome_consensus_summary.tsv"))


rule build_alignment_fasta:
    input:
        fasta=expand(os.path.join(dir_reports, "mitogenome", "{sample}_consensus.fasta"), sample=sample_names)
    output:
        final_fasta = os.path.join(dir_hostcleaned, "mitogenome", "final_mitogenome.aln")
    params:
        folder=os.path.join(dir_hostcleaned, "mitogenome")
    conda:
        os.path.join(dir_env, "mafft.yaml")
    shell:
        """
        set -euo pipefail
        if [ -f {output.final_fasta} ]; then
            echo "Final alignment fasta already exists. Skipping..."
            exit 0
        else
            cat {params.folder}/*_consensus.fasta > all_samples.fasta
            #snp-sites -o {output.final_fasta} all_samples.fasta
            mafft --auto all_samples_WRef.fasta > {output.final_fasta}
        fi
        """

rule phylo_tree:
    input:
        final_fasta = os.path.join(dir_hostcleaned, "mitogenome", "final_mitogenome.aln")
    output:
        tree = os.path.join(dir_reports, "mitogenome", "mitogenome_phylo_tree.nwk")
    conda:
        os.path.join(dir_env, "mafft.yaml")
    shell:
        """
        set -euo pipefail
        if [ -f {output.tree} ]; then
            echo "Phylogenetic tree already exists. Skipping..."
            exit 0
        else
            iqtree -s {input.final_fasta} -m GTR+G -bb 3000 -nt AUTO
        fi
        """