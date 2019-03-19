# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__="Patrik Smeds"
__copyright__ = "Copyright 2019, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "MIT"

_multiqc_input = lambda wildcards: ["report/fastqc", "qc/samstats"]
try:
    _multiqc_input = multiqc_input
except:
    pass


_multiqc_output = "reports/multiqc.html"
try:
    _multiqc_output = multiqc_output
except:
    pass

from src.lib.data.report.wp1 import generate_filtered_mutations

rule multiqc:
  input:
      _multiqc_input
  output:
      _multiqc_output
  log:
      "logs/multiqc.log"
  params:
      ""
  wrapper:
      "0.31.1/bio/multiqc"
