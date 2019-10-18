# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__="Patrik Smeds"
__copyright__ = "Copyright 2019, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "MIT"

coveragebed_input = "alignment/{sample}.{part}.bam"
try:
    _coveragebed_input = coveragebed_input
except:
    pass

_coveragebed_output = "alignment/{sample}.{part}.coverage"
try:
    _coveragebed_output = coveragebed_output
except:
    pass

rule bedtools_coverage:
    input:
        b=coveragebed_input,
        a=lambda wildcards: samples['analyzable_region'][wildcards.sample]
    output:
        bam=_coveragebed_output
    log:
        "logs/bedtools/coveragebed/{sample}.{part}.log"
    params:
        extra=""
    threads: 16
    wrapper:
        "coverageBed-wrapper/bio/bedtools/coveragebed"

rule bedtools_coverage_d:
    input:
        b=coveragebed_input,
        a=lambda wildcards: samples['analyzable_region'][wildcards.sample]
    output:
        bam=_coveragebed_output + "_d"
    log:
        "logs/bedtools/coveragebed/{sample}.{part}.d.log"
    params:
        extra="-d"
    threads: 16
    wrapper:
        "coverageBed-wrapper/bio/bedtools/coveragebed"
