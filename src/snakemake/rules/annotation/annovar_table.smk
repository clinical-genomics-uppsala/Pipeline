# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__="Patrik Smeds"
__copyright__ = "Copyright 2019, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "MIT"

"""
 Rule that uses annovar to annotate vcf file.

 https://github.com/charite/jannovar

 Input, output and config
 ------------------------------------------------------------------------------
 Input variable: _annovar_input: optional
     Default:
         "variants/{sample}.{part}.vcf"
 Output variable:  _annovar_output: optional
     Default:
         "variants/{sample}.{part}.annovar.vcf"

 Overriding input and output
 ------------------------------------------------------------------------------
 Required wildcards:
    sample
    part

 Override input format
 Ex
  annovar_input = "variants/{sample}.{unit}.vcf"

 Override output format
 Ex
  annovar_output = "variants/{sample}.{part}.an.vcf"

"""

_annovar_table_input = "variants/{sample}.{part}.vcf"
try:
    _annovar_table_input = annovar_table_input
except:
    pass

_annovar_table_output = "variants/{sample}.{part,\d+}.annovar.vcf"
try:
    _annovar_table_output = annovar_table_output
except:
    pass

rule annovar_annotate:
    input:
        _annovar_table_input
    output:
        _annovar_table_output
    log:
        "logs/annotate/qsort/{sample}.{part}.annovar.log"
    params:
        database_path = config['annovar']["database"],
        protocols = ",".join([protocol for protocol in config['annovar']['protocols']]),
        operations = ",".join([config['annovar']['protocols'][protocol]['operation'] for protocol in config['annovar']['protocols']]),
        arguments = ",".join([config['annovar']['protocols'][protocol]['argument'] for protocol in config['annovar']['protocols']]),
        build_version = config['annovar']['build_version'],
        extra=""         # optional parameters
    wrapper:
       "annovar-wrapper/bio/annovar/table_annovar"
