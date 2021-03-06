# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__="Patrik Smeds"
__copyright__ = "Copyright 2019, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "MIT"

_fastqc_input = lambda wildcards: units.loc[(path.split(wildcards.sample)[-1], wildcards.unit), [wildcards.read]].dropna()[0]
try:
    _fastqc_input = fastqc_input
except:
    pass

_fastqc_output_html = "qc/fastqc/{sample}.{unit}.{read}_fastqc.html"
try:
    _fastqc_output_html = fastqc_output_html
except:
    pass

_fastqc_output_zip = "qc/fastqc/{sample}.{unit}.{read}_fastqc.zip"
try:
    _fastqc_output_html = fastqc_output_zip
except:
    pass

rule fastqc:
    input:
        _fastqc_input
    output:
        html=_fastqc_output_html,
        zip=_fastqc_output_zip
    params: ""
    log:
        "logs/fastqc/{sample}.{unit}.{read}.log"
    wrapper:
        "master/bio/fastqc"
