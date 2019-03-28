# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__="Patrik Smeds"
__copyright__ = "Copyright 2019, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "MIT"

"""

 Rule that performs alignment of reads using bwa mem.

 Input, output and config
 ------------------------------------------------------------------------------
 Input variable: _bwa_alignment_input: optional
     Default:
        "fastq/{sample}.{part,[A-Za-z0-9]+-\d{4}}.fq1.fastq.gz",
        "fastq/{sample}.{part,[A-Za-z0-9]+-\d{4}}.fq2.fastq.gz"
 Output variable:  bwa_mem_output: optional
     Default:
         "alignment/{sample}.{unit}.{part}.bam"

 Config dict keys: values
 'reference_genome' : "path/to/reference_genome": required

 Overriding input and output
 ------------------------------------------------------------------------------
 Required wildcards:
    sample
    part

 Override input format
 Ex
  bwa_alignment_input = [
                       "fastq/{sample}.{unit}.{part}.R1.cutadapt.fastq.gz",
                       "fastq/{sample}.{unit}.{part}.R2.cutadapt.fastq.gz"
                       ]

 Override output format
 Ex
   bwa_alignment_output = "alignment/{sample}.{part}.cutadapt.bam"

"""

def get_now():
    from datetime import datetime
    return datetime.now().strftime('%Y%m%d')

_bwa_mem_input = ["fastq/{sample}.{part,[A-Za-z0-9]+-\d{4}}.fq1.fastq.gz", "fastq/{sample}.{part,[A-Za-z0-9]+-\d{4}}.fq2.fastq.gz"]
try:
      _bwa_mem_input = bwa_mem_input
except:
      pass

_bwa_mem_output = temp("alignment/{sample}.{part,[A-Za-z0-9]+-\d{4}}.bam")
try:
    _bwa_mem_output = bwa_mem_output
except:
    pass

rule bwa_mem:
    input:
        reads=_bwa_mem_input
    output:
        _bwa_mem_output
    log:
        "logs/bwa_mem/{sample}.{part}.log"
    threads: 16
    params:
        index=config['reference_genome'],
        extra=lambda wildcards: r"-M -R '@RG\tID:%s_%s\tSM:%s\tPL:%s'" % (get_now(),wildcards.sample, wildcards.sample, "illumina",),
        sort="samtools",
        sort_order="coordinate",
        sort_extra="-@ 16"
    wrapper:
        "0.27.1/bio/bwa/mem"
