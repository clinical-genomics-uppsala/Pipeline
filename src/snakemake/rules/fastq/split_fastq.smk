# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__="Patrik Smeds"
__copyright__ = "Copyright 2019, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "MIT"

#
# Rule that split a fastq file into num_fastq_split parts, and compress the
# resulting files. The fastq file will only be copied if num_split is set to a
# value below 2
#
# Input, output and config
# ------------------------------------------------------------------------------
# Input variable: _split_fastq_input: optional
#     Default:
#       lambda wildcards: "fastq/‰s.%s.%s.fastq.gz" % (wildcards.sample, wildcards.unit)
# Output variable:  _split_fastq_output: optional
#     Default:
#         temp("fastq/{sample}.{part,[A-Za-z0-9]+-\d{4}}.{read}.fastq.gz")
#
# Config dict keys: values
# 'num_fastq_split' : 5 : default 1
#
#
# Overriding input and output
# ------------------------------------------------------------------------------
# Required wildcards:
#    sample
#    part
#    read
#
# Override input format
# Ex
#  split_fastq_input = "fastq_split/{sample}.{unit}.{read}.fastq.gz"
#
# Override output format
# Ex
#  split_fastq_output = "fastq_split/{sample}.{part}.{read}.fastq.gz"
#

def _cgu_get_num_splits(config):
    """
        extract number of times a fastq file should be split.
    """
    return int(config.get("num_fastq_split",1))

_split_fastq_input = lambda wildcards: units.loc[(wildcards.sample, wildcards.unit),wildcards.read]
try:
    _split_fastq_input = split_fastq_input
except:
    pass

localrules: link_fastq

if _cgu_get_num_splits(config) > 1:
    _split_fastq_output = temp("fastq/{sample}.{part,[A-Za-z0-9]+-\d{4}}.{read}.fastq.gz")
    try:
        _split_fastq_output = split_fastq_output
    except:
      pass

    rule count_lines_in_fastq:
        input:
            _split_fastq_input
        output:
          temp("fastq/{sample}.{unit}.{read}.var")
        run:
          import subprocess, os
          lines = int(float(subprocess.run("gunzip -c " + input[0] + " |  wc -l | awk '{print $1/4}'", stdout=subprocess.PIPE,shell=True).stdout.decode('utf-8').rstrip("\n")))
          storage.store(wildcards.sample + "." + wildcards.unit + "." + wildcards.read + ".var",str(lines))
          shell("echo 'reads: '" + str(lines) + "'' > "  + output[0])

    rule split_fastq:
        input:
          _split_fastq_input,
          "fastq/{sample}.{unit}.{read}.var"
        output:
          temp(['fastq/{sample}.{unit}-%04d.{read}.fastq' % num for num in range(0,_cgu_get_num_splits(config))])
        params:
          output_prefix=lambda wildcards: "fastq/" + wildcards.sample + "." + wildcards.unit + "-",
          output_suffix=lambda wildcards: "." + wildcards.read + ".fastq"
        run:
          import math
          num_reads = int(storage.fetch(wildcards.sample + "." + wildcards.unit + "." + wildcards.read + ".var"))
          num_split = _cgu_get_num_splits(config)
          lines_per_file = 4*math.ceil(num_reads / num_split)
          shell('gunzip -c {input[0]} | awk \'BEGIN{{ file = 0; filename = sprintf("{params.output_prefix}%.04d{params.output_suffix}", file) }}{{ print > filename}} NR % {lines_per_file} == 0 {{ close(filename); file = file + 1; filename = sprintf("{params.output_prefix}%.04d{params.output_suffix}",file)}}\'')
          num_files_generated = 4*math.floor(num_reads / lines_per_file)
          while num_files_generated < num_split:
            shell("touch {params.output_prefix}%04d{params.output_suffix}" % num_split)
            num_split -= 1

    rule compress_split_fastq:
        input:
            "fastq/{sample}.{part}.{read}.fastq"
        output:
            _split_fastq_output
        run:
            shell("gzip -c {input} > {output}")
else:
    _split_fastq_output = temp("fastq/{sample}.{unit,[A-Za-z0-9]+}-0000.{read}.fastq.gz")
    try:
        _split_fastq_output = split_fastq_output
    except:
        pass

    rule link_fastq:
        input:
            _split_fastq_input
        output:
            _split_fastq_output
        run:
            shell("ln -s {input} {output}")
