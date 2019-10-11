# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__="Patrik Smeds"
__copyright__ = "Copyright 2019, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "MIT"

"""
 Rule that adds UMIs to a bam file.

 http://fulcrumgenomics.github.io/fgbio/tools/latest/AnnotateBamWithUmis.html

 Input, output and config
 ------------------------------------------------------------------------------
 Input variable:
        _annotate_bam_with_umi_input_bam: optional
        _annotate_bam_with_umi_input_umi: optional
     Default:
         bam="alignment/{sample}.{part}.bam",
         umi="umis/{sample}.{part}.UMIs.fastq"
 Output variable:  _annotate_bam_with_umi_output: optional
     Default:
         "alignment/{sample}.{part,[A-Za-z0-9]+-\d{4}}.umiAnnoBam.bam"

 Overriding input and output
 ------------------------------------------------------------------------------
 Required wildcards:
    sample
    part

 Override input format
 Ex
  annotate_bam_with_umi_input_bam = "alignment/{sample}.{part}.merged.bam"
  annotate_bam_with_umi_input_umi = "fastq/{sample}.{part}.fq2.fastq"

 Override output format
 Ex
  annotate_bam_with_umi_output = "alignment/{sample}.{part}.final.gz"

"""

_annotate_bam_with_umi_input_bam = "alignment/{sample}.{part}.bam"
try:
    _annotate_bam_with_umi_input_bam = annotate_bam_with_umi_input_bam
except:
    pass

_annotate_bam_with_umi_output = "alignment/{sample}.{part,[A-Za-z0-9]+-\d{4}}.umiAnnoBam.bam"
try:
  _annotate_bam_with_umi_output = annotate_bam_with_umi_output
except:
    pass

_annotate_bam_with_umi_input_umi = "umis/{sample}.{part}.UMIs.fastq"
try:
    _annotate_bam_with_umi_input_umi = annotate_bam_with_umi_input_umi
except:
    pass

rule annotate_bam_with_umi:
    input:
        bam=_annotate_bam_with_umi_input_bam,
        umi=_annotate_bam_with_umi_input_umi
    output:
        bam=_annotate_bam_with_umi_output
    log:
        "logs/umis/{sample}.{part}.fgbioAnnoBam.txt"
    wrapper:
        "0.31.1/bio/fgbio/annotatebamwithumis"
