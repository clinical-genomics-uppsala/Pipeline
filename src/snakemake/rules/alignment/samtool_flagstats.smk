# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__="Patrik Smeds"
__copyright__ = "Copyright 2019, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "MIT"

_samtools_flagstat_input = "alignment/{sample}.{part}.bam"
try:
    _samtools_flagstat_input = samtools_flagstat_input
except:
    pass

_samtools_flagstat_output = "qc/samtool/flagstats/{sample}.{part}.bam.flagstats"
try:
    _samtools_flagstat_output = samtools_flagstat_output
except:
    pass

rule samtool_flagstat:
    input:
      _samtools_flagstat_input
    output:
      _samtools_flagstat_output
    wrapper:
        "0.31.1/bio/samtools/flagstat"
