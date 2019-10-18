# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__="Patrik Smeds"
__copyright__ = "Copyright 2019, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "MIT"

_samtool_stats_input = "alignment/{sample}.{part}.bam"
try:
    _samtool_stats_input = samtool_stats_input
except:
    pass

_samtool_stats_output = "qc/samtool/stats/{sample}.{part}.bam.stats.txt"
try:
    _samtool_stats_output = samtool_stats_output
except:
    pass

def get_region_file(wildcards):
    try:
        samples['region'][wildcards.sample]
    except:
        return ""

lambda wildcards: samples[''][wildcards.sample]

rule samtool_stats:
    input:
        _samtool_stats_input
    output:
        _samtool_stats_output
    params:
        extra="",
        region=lambda wildcards: get_region_file(wildcards)
    log:
        "logs/samtools_stats/{sample}.{part}.log"
    wrapper:
        "0.31.1/bio/samtools/stats"
