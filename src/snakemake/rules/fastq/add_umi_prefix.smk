# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__="Patrik Smeds"
__copyright__ = "Copyright 2019, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "MIT"

_headcrop_read_input = "fastq/{sample}.{part}.{read}.fastq.gz"
try:
    _headcrop_read_input = headcrop_read_input
except:
    pass

_add_umi_prefix_output = temp("fastq/{sample}.{part,[A-Za-z0-9]+-\d{4}}.{read}.umi.fastq.gz")
try:
    _add_umi_prefix_output = add_umi_prefix_output
except:
    pass

rule add_umi_prefix:
   input:
        _headcrop_read_input
   output:
        _add_umi_prefix_output
   log:
      "logs/trimming/{sample}.{part}.{read}.headcro.trimmomatic.txt"
   shell:
      "awk \'{{}}"
