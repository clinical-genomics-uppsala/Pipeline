# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__="Patrik Smeds"
__copyright__ = "Copyright 2019, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "MIT"

"""
 Rule that removes N first bases from a read in a fastq file

 Input, output and config
 ------------------------------------------------------------------------------
 Input variable: _headcrop_read_input: optional
     Default:
         "fastq/{sample}.{part}.{read}.fastq.gz"
 Output variable:  _headcrop_read_output: optional
     Default:
         temp("fastq/{sample}.{part,[A-Za-z0-9]+-\d{4}}.{read}.headcrop{num_bases}.fastq.gz")


 Overriding input and output
 ------------------------------------------------------------------------------
 Required wildcards:
    sample
    part
    read
    num_bases

 Override input format
 Ex
  headcrop_read_input = "fastq/{sample}.{part}.{read}.fastq.gz"

 Override output format
 Ex
  headcrop_read_output = "fastq/{sample}.{part}.{read}.headcrop.{num_bases}.fastq.gz"
"""

_headcrop_read_input = "fastq/{sample}.{part}.{read}.fastq.gz"
try:
    _headcrop_read_input = headcrop_read_input
except:
    pass

_headcrop_read_output = temp("fastq/{sample}.{part,[A-Za-z0-9]+-\d{4}}.{read}.headcrop{num_bases,\d+}.fastq.gz")
try:
    _headcrop_read_output = headcrop_read_output
except:
    pass

rule headcrop_read:
   input:
        _headcrop_read_input
   output:
        _headcrop_read_output
   log:
      "logs/trimming/{sample}.{part}.{read}.headcro{num_bases}.trimmomatic.txt"
   params:
      trimmer= lambda wildcards: ["HEADCROP:%s" %  config.get('fastq_headcrop_n',wildcards.num_bases)],
      extra=""
   threads: 8
   wrapper:
      "0.31.1/bio/trimmomatic/se"
