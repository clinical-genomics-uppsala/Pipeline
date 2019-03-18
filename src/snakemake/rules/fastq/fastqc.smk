# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


_fastqc_input = [lambda wildcards: units.loc[(path.split(wildcards.sample)[-1], wildcards.unit), ["fq1"]].dropna()[0], lambda wildcards: units.loc[(path.split(wildcards.sample)[-1], wildcards.unit), ["fq2"]].dropna()[0]]
try:
    _fastqc_input = fastqc_input
except:
    pass

_fastqc_output_html = "qc/fastqc/{sample}.{unit}.html"
try:
    _fastqc_output_html = fastqc_output_html
except:
    pass

_fastqc_output_zip = "qc/fastqc/{sample}.{unit}.zip"
try:
    _fastqc_output_html = fastqc_output_zip
except:
    pass

rule import_read:
    input:
        _fastqc_input
    output:
        html=_fastqc_output_html,
        zip=_fastqc_output_zip
    params: ""
    log:
        "logs/fastqc/{sample}.{unit}.log"
    wrapper:
        "0.31.1/bio/fastqc"
