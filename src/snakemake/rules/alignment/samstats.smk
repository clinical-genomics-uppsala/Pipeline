_samstats_input = "alignment/{sample}.{part}.bam"
try:
    _samstats_input = samstats_input
except:
    pass

_samstats_output = "qc/samstats/{sample}.{part}.bam.html"
try:
    _samstats_output = samstats_output
except:
    pass

def get_region_file(wildcards):
    try:
        samples['region'][wildcards.sample]
    except:
        return ""

lambda wildcards: samples[''][wildcards.sample]

rule samtools_stats:
    input:
        _samstats_input
    output:
        _samstats_output
    params:
        extra="",                       # Optional: extra arguments.
        region=lambda wildcards: get_region_file(wildcards)      # Optional: region string.
    log:
        "logs/samtools_stats/{sample}.{part}.log"
    wrapper:
        "0.31.1/bio/samtools/stats"
