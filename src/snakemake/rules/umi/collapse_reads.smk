# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__="Patrik Smeds"
__copyright__ = "Copyright 2019, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "MIT"

"""
 Rule uses umis stored in bam file to collapse reads.

 http://fulcrumgenomics.github.io/fgbio/tools/latest/AnnotateBamWithUmis.html

 Input, output and config
 ------------------------------------------------------------------------------
 Input variable: _collapse_reads_input: optional
     Default:
         "alignment/{sample}.{part}.sorted.bam"
 Output variable:  _collapse_reads_output: optional
     Default:
         "alignment/{sample}.{part,\d{4}}-consensus{num_support,\d+}.bam"

 Config dict keys: values
 'config: 'reference_genome'


 Overriding input and output
 ------------------------------------------------------------------------------
 Required wildcards:
    sample
    part
    num_support

 Override input format
 Ex
  collapse_reads_input = "alignment/{sample}.{part}.merged.bam"

 Override output format
 Ex
  collapse_reads_output = "alignment/{sample}.{part,\d{4}}.conse.{num_support,\d+}.bam"

"""

_collapse_reads_input = "alignment/{sample}.{part}.sorted.bam"
try:
    _collapse_reads_input = collapse_reads_input
except:
    pass

_collapse_reads_output = "alignment/{sample}.{part,\d{4}}-consensus{num_support,\d+}.bam"
try:
    _collapse_reads_output = collapse_reads_output
except:
    pass


rule pre_collapse_prep_revertsam:
    input:
        _collapse_reads_input
    output:
        temp("alignment/.{sample}.{part,\d{4}}.sanitised.bam")
    params:
        extra="SANITIZE=true REMOVE_DUPLICATE_INFORMATION=false REMOVE_ALIGNMENT_INFORMATION=false"
    wrapper:
        "0.31.1/bio/picard/revertsam"

rule pre_collapse_prep_set_mate_info:
    input:
        bam="alignment/.{sample}.{part}.sanitised.bam"
    output:
        bam=temp("alignment/.{sample}.{part,\d{4}}.sanitised.setmateinfo.bam")
    log:
        "logs/umis_collapse/{sample}.{part}.fgbioSetMQInfo.txt"
    wrapper:
        "0.31.1/bio/fgbio/setmateinformation"

rule pre_collapse_prep_sort_queryname:
    input:
      "alignment/.{sample}.{part}.sanitised.setmateinfo.bam"
    output:
        temp("alignment/.{sample}.{part,\d{4}}.fgbio.merged.setmateinfo.qsorted.bam")
    params:
        sort_order="queryname",
        extra="VALIDATION_STRINGENCY=LENIENT -Xms500m -Xmx20g"
    threads: 8
    wrapper:
        "0.31.1/bio/picard/sortsam"

rule pre_collapse_prep_sanitised_queryname:
    input:
      "alignment/.{sample}.{part}.fgbio.merged.setmateinfo.qsorted.bam"
    output:
        temp("alignment/.{sample}.{part,\d{4}}.fgbio.merged.setmateinfo.qsorted.sanitised.bam")
    log:
        "logs/umi_collapse/qsort/{sample}.{part}.log"
    params:
        sort_order="queryname",
        extra="VALIDATION_STRINGENCY=LENIENT -Xms500m -Xmx20g"
    threads: 8
    wrapper:
        "0.31.1/bio/picard/sortsam"

rule pre_collapse_prep_groupreads_by_umi:
    input:
        "alignment/.{sample}.{part}.fgbio.merged.setmateinfo.qsorted.sanitised.bam"
    output:
        temp("alignment/.{sample}.{part,\d{4}}.merged.sanitised.setmateinfo.groupreads.bam")
    log:
        "logs/umis_collapse/{sample}.{part}.merged.fgbioGroup.txt"
    params:
        extra="-s adjacency --edits 1 -Xms500m -Xmx64g"
    wrapper:
        "0.31.1/bio/fgbio/groupreadsbyumi"

rule collapse_reads_create_consensus_reads:
    input:
        "alignment/.{sample}.{part}.merged.sanitised.setmateinfo.groupreads.bam"
    output:
        temp("alignment/.{sample}.{part,\d{4}}.merged.sanitised.setmateinfo.groupread.consensus{num_support,\d+}.bam")
    log:
       "logs/umis_collapse/{sample}.{part}.merged.fgbioCMCR-{num_support}.txt"
    params:
        extra=lambda wildcards: "-M " + wildcards.num_support
    wrapper:
        "0.31.1/bio/fgbio/callmolecularconsensusreads"

rule create_reads_from_bam_consensus:
    input:
        "alignment/.{sample}.{part}.merged.sanitised.setmateinfo.groupread.consensus{num_support}.bam"
    output:
        fastq1=temp("fastq/{sample}.{part,\d{4}}.consensus{num_support,\d+}.R1.fq"),
        fastq2=temp("fastq/{sample}.{part}.consensus{num_support,\d+}.R2.fq")
    wrapper:
        "0.31.1/bio/picard/samtofastq"

rule bwa_concensus:
    input:
        reads=["fastq/{sample}.{part}.consensus{num_support}.R1.fq", "fastq/{sample}.{part}.consensus{num_support}.R2.fq"]
    output:
        _collapse_reads_output
    log:
        "logs/bwa_mem/{sample}.{part}.consensus{num_support}.log"
    threads: 3
    params:
        index=config['reference_genome'],
        extra=lambda wildcards: r"-M -R '@RG\tID:" + get_now() + "_" + wildcards.sample + r"\tSM:" + wildcards.sample + r"\tPL:illumina'",
        sort="samtools",
        sort_order="queryname",
        sort_extra="-@ 3"
    wrapper:
        "master/bio/bwa/mem"
