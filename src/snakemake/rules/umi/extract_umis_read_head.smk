# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__="Patrik Smeds"
__copyright__ = "Copyright 2019, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "MIT"

"""
 Rule that exract N first bases from a read in a fastq file

 Input, output and config
 ------------------------------------------------------------------------------
 Input variable: _extract_umis_read_head_input: optional
     Default:
         "fastq/{sample}.{part}.{read}.fastq.gz"
 Output variable:  _extract_umis_read_head_output: optional
     Default:
         temp("umis/{sample}.{part,[A-Za-z0-9]+-\d{4}}.{read}.UMIs{num_bases}.fastq")


 Overriding input and output
------------------------------------------------------------------------------
Required wildcards:
   sample
   part
   read
   num_bases

 Override input format
 Ex
  extract_umis_read_head_input = "fastq/{sample}.{part}.{read}.fq.gz"

 Override output format
 Ex
  extract_umis_read_head_output = "fastq/{sample}.{part}.{read}.umi.{num_bases}.fastq"

"""

_extract_umis_read_head_input = "fastq/{sample}.{part}.{read}.fastq.gz"
try:
    _extract_umis_read_head_input = extract_umis_read_head_input
except:
    pass

_extract_umis_read_head_output = temp("umis/{sample}.{part,[A-Za-z0-9]+-\d{4}}.{read}.UMIs.fastq")
try:
    _extract_umis_read_head_output = extract_umis_read_head_output
except:
    pass

rule extract_umis_read_head:
    input:
       _extract_umis_read_head_input
    output:
       _extract_umis_read_head_output
    params:
      trimmer=["CROP:%s" % config.get("fastq_extract_n_bases", 8)],
       extra=""
    threads: 8
    wrapper:
       "0.31.1/bio/trimmomatic/se"
