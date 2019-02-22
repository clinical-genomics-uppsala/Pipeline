# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__="Patrik Smeds"
__copyright__ = "Copyright 2019, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "MIT"

"""
 Rule that update read group information in bam file.

 Input, output and config
 ------------------------------------------------------------------------------
 Input variable: _add_replace_read_group_input: optional
     Default:
         "alignment/{sample}.{part}.rg.bam"
 Output variable:  _add_replace_read_group_output: optional
     Default:
         "alignment/{sample}.{part}.rg.bam"


Overriding input and output
------------------------------------------------------------------------------
Required wildcards:
   sample
   part

 Override input format
 Ex
  add_replace_read_group_input = "alignment/{sample}.{part}.rg.bam"

 Override output format
 Ex
  add_replace_read_group_output = add_replace_read_group_output

"""

def get_now():
    from datetime import datetime
    return datetime.now().strftime('%Y%m%d')

_add_replace_read_group_input = "alignment/{sample}.{part}.rg.bam"
try:
    _add_replace_read_group_input = add_replace_read_group_input
except:
      pass

_add_replace_read_group_output = "alignment/{sample}.{part}.rg.bam"
try:
    _add_replace_read_group_output = add_replace_read_group_output
except:
    pass

rule add_replace_read_group:
    input:
        _add_replace_read_group_input
    output:
        _add_replace_read_group_output
    log:
        "logs/alignment/{sample}.{part}.addReplaceReadgroup.log"
    params:
        lambda wildcards: "RGLB=lib1 RGID=%s_%s RGSM=%s RGPL=%s RGPU=%s CREATE_INDEX=TRUE SO=coordinate VALIDATION_STRINGENCY=STRICT" % (get_now(), wildcards.sample, wildcards.sample, "illumina","MiSeq")
    wrapper:
        "0.31.1/bio/picard/addorreplacereadgroups"
