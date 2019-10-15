# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__="Patrik Smeds"
__copyright__ = "Copyright 2019, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "MIT"

"""
 Rule that uses primerclip to remove primer sequence from a bam file

 https://github.com/swiftbiosciences/primerclip

 Input, output and config
 ------------------------------------------------------------------------------
 Input variable: primerclip_input: optional
     Default:
         "alignment/{sample}.{part}.bam"
 Output variable:  _primerclip_output: optional
     Default:
         "alignment/{sample}.{part}.primerclip.bam"

 samples.tsv need to contain a column named master_primerclip_file containing
 path to a primerclip input file.

 Overriding input and output
 ------------------------------------------------------------------------------
 Required wildcards:
    sample
    part

 Override input format
 Ex
  primerclip_input = "alignment/{sample}.{unit}.merged.bam"

 Override output format
 Ex
  primerclip_output = "alignment/{sample}.{part}.final.gz"

"""

_primerclip_input = "alignment/{sample}.{part}.bam"
try:
    _primerclip_input = primerclip_input
except:
    pass

_primerclip_output = "alignment/{sample}.{part}.primerclip.bam"
try:
    _primerclip_output = primerclip_output
except:
    pass

rule sort_preprimerclip_queryname:
    input:
        _primerclip_input
    output:
        temp("alignment/.{sample}.{part}.tmp-qsorted.sam")
    log:
        "logs/umi/qsort/{sample}.{part}.log"
    params:
        sort_order="queryname",
        extra="VALIDATION_STRINGENCY=LENIENT"
    wrapper:
        "0.31.1/bio/picard/sortsam"

rule primerclip:
    input:
        sam="alignment/.{sample}.{part}.tmp-qsorted.sam",
        master_file=lambda wildcards: samples['master_primerclip_file'][wildcards.sample]
    output:
        sam=temp("alignment/.{sample}.{part}.tmp-qsorted-primerclip.sam")
    log:
        "logs/primerclip/{sample}.{part}.log"
    wrapper:
        "master/bio/primerclip"

rule primerclip_bam_generation:
    input:
        "alignment/.{sample}.{part}.tmp-qsorted-primerclip.sam"
    output:
        _primerclip_output#_primerclip_output[:]
    params:
        "-Sb"
    wrapper:
        "master/bio/samtools/view"
