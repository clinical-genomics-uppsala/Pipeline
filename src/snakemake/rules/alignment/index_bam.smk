_bam_index_input = "alignment/{sample}.{processed}.bam"
try:
    _bam_index_input = bam_index_input
except:
    pass

_bam_index_output = "alignment/{sample}.{processed}.bam.bai"
try:
    _bam_index_output = bam_index_out
except:
    pass

rule bam_index:
    input:
        _bam_index_input
    output:
        _bam_index_output
    params:
        ""
    wrapper:
        "0.31.1/bio/samtools/index"
