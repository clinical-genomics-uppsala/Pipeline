# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__="Patrik Smeds"
__copyright__ = "Copyright 2019, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "MIT"

_extract_umis_read_queryname_input = "fastq/{sample}.{part}.{read}.fastq.gz"
try:
    _extract_umis_read_queryname_input = extract_umis_read_queryname_input
except:
    pass

_extract_umis_read_queryname_output = temp("umis/{sample}.{part,[A-Za-z0-9]+-\d{4}}.{read}.UMIs.fastq")
try:
    _extract_umis_read_queryname_output = extract_umis_read_queryname_output
except:
    pass

rule extract_umis_read_queryname:
    input:
       _extract_umis_read_queryname_input
    output:
       _extract_umis_read_queryname_output
    shell:
        "zcat {input} | awk \'BEGIN{{FS=\" \"}}{{split($1,a,\":\");print($1);print(a[length(a)]);print(\"+\");print(a[length(a)]);getline;getline;getline}}\' > {output}"
