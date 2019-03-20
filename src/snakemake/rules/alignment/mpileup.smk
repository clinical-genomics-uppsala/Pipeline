# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__="Patrik Smeds"
__copyright__ = "Copyright 2019, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "MIT"

"""
 Rule used to create pileup files.

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
  merge_bam_input = "alignment/{sample}.{part}.bam"

"""

_mpileup_input = lambda wildcards: ["alignment/{}.{}.bam".format(wildcards.sample,unit_part) for unit_part in _get_split_and_unit_part_files(wildcards,units,config)]
try:
    _mpileup_input = mpileup_input
except:
    pass

_mpileup_output = "pileup/{sample}.{part}.mpileup.gz"
try:
    _mpileup_output = mpileup_output
except:
    pass

rule mpileup:
    input:
        bam=_mpileup_input,
        reference_genome=config['reference_genome']
    output:
        _mpileup_output
    params:
        extra="-d 100000000000 -ABQ0"
    threads: 8
    wrapper:
        "master/bio/samtools/mpileup"
