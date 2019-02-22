# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__="Patrik Smeds"
__copyright__ = "Copyright 2019, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "MIT"

"""
 Rule merge multiple bam files into one.

 Input, output and config
 ------------------------------------------------------------------------------
 Input variable: _merge_bam_input: optional
     Default:
         lambda wildcards: ["alignment/{}.{}.bam".format(wildcards.sample,unit_part) for unit_part in _get_split_and_unit_part_files(wildcards,units,config)]
 Output variable:  _merge_bam_output: optional
     Default:
        "alignment/{sample}.{part,[0]{4}}.sorted.bam"


 Overriding input and output
------------------------------------------------------------------------------
Required wildcards:
   sample
   part

 Override input format
 Ex
  merge_bam_input =  ["alignment/{sample}.L001.bam", "alignment/{sample}.L002.bam"]

 Override output format
 Ex
  merge_bam_input = "alignment/{sample}.{part,[0]{4}}.bam"

"""

def _get_split_and_unit_part_files(wildcards, units, config):
    num_splits = config.get("num_fastq_split", 1)
    if num_splits > 1:
        return [ unit + "-%04d" % part for part in range(0,num_splits) for unit in units.loc[wildcards.sample].index]
    else:
        return [ unit + "-0000" for unit in units.loc[wildcards.sample].index]


# ToDo See if part requirement can be remove.

_merge_bam_input = lambda wildcards: ["alignment/{}.{}.bam".format(wildcards.sample,unit_part) for unit_part in _get_split_and_unit_part_files(wildcards,units,config)]
try:
    _merge_bam_input = merge_bam_input
except:
    pass

_merge_bam_output = "alignment/{sample}.{part,[0]{4}}.sorted.bam"
try:
    _merge_bam_output = merge_bam_output
except:
    pass

rule merge_bam:
    input:
        _merge_bam_input
    output:
        temp("alignment/.tmp/{sample}.0000.merged.bam")
    threads: 8
    wrapper:
        "0.31.1/bio/samtools/merge"

rule coordinate_sort_mapped_reads_merged:
    input:
        "alignment/.tmp/{sample}.0000.merged.bam"
    output:
        _merge_bam_output
    params:
        sort_order="coordinate"
    wrapper:
        "0.31.1/bio/picard/sortsam"

rule create_bam_index_merged:
    input:
        _merge_bam_output
    output:
        _merge_bam_output + ".bai"
    params:
        ""
    wrapper:
        "0.31.1/bio/samtools/index"
