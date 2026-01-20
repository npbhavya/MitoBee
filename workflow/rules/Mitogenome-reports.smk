"""
Snakemake rules for generating mitogenome reports
"""
from glob import glob

rule mitogenome_summary:
    """
    Summarize mitogenome consensus FASTA:
    filename, header, length, GC content
    """
    input:
        fasta=os.path.join(dir_reports, "mitogenome", "{sample}_consensus.fasta")
    output:
        summary=os.path.join(dir_reports, "mitogenome", "{sample}_consensus.summary.tsv")
    params:
        sample="{sample}"
    shell:
        r"""
        awk -v fname="{params.sample}_consensus.fasta" '
            /^>/ {
                header = substr($0, 2)
                next
            }
            {
                seq = seq $0
            }
            END {
                len = length(seq)
                gc = gsub(/[GCgc]/, "", seq)
                gc_pct = (len > 0) ? (gc / len * 100) : 0
                printf "%s\t%s\t%d\t%.2f\n", fname, header, len, gc_pct
            }
        ' {input.fasta} > {output.summary}
        """

rule mitogenome_reports_aggregate:
    input:
        expand(os.path.join(dir_reports, "mitogenome", "{sample}_consensus.summary.tsv"), sample=sample_names)
    output:
        aggregate=os.path.join(dir_reports, "mitogenome_consensus_summary.tsv")
    shell:
        """
        echo -e "filename\theader\tlength\tGC_content" > {output.aggregate}
        cat {input} >> {output.aggregate}
        """