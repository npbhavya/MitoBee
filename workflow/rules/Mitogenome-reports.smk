"""
Snakemake rules for generating mitogenome reports
"""

rule mitogenome_output:
    """
    Generate a mitogenome report using provided data.
    """
    input:
        consensus_fasta = os.path.join(dir_hostcleaned, "mitogenome", "{sample}_consensus.fasta")
    output:
        report=os.path.join(dir_reports, "mitogenome_reports", "{sample}_consensus.fasta"))
    shell:
        """
        cp {input.consensus_fasta} {output.report}
        """


rule mitogenome_summary:
    """
    Summarize mitogenome consensus FASTA:
    filename, header, length, GC content
    """
    input:
        fasta=os.path.join(dir_reports, "mitogenome_reports", "{sample}_consensus.fasta")
    output:
        summary=os.path.join(dir_reports, "mitogenome_reports", "{sample}_consensus.summary.tsv")
    params:
        sample="{sample}"
    shell:
        r"""
        sample="{params.sample}"
        awk -v fname="{{sample}}_consensus.fasta" '
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
