# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__="Patrik Smeds"
__copyright__ = "Copyright 2019, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "MIT"

"""

 Rule that performs gatk3 indel realignment



"""

_gatk3_indelrealignment_input = "alignment/{sample}.{part}.bam"
try:
      _gatk3_indelrealignment_input = gatk3_indelrealignment_input
except:
      pass

_gatk3_indelrealignment_output = "alignment/{sample}.{part}.bam"
try:
      _gatk3_indelrealignment_output = gatk3_indelrealignment_output
except:
      pass

rule gatk3_realignertargetcreator:
    input:
        bam=_gatk3_indelrealignment_input,
        ref=config['reference_genome'],
        known=config['known_sites']
    output:
        temp("alignment/{sample}.{part}.indelrealignment.intervals")
    log:
        "logs/gatk/indelrealignment/{sample}.{part}.intervalscreator.log"
    params:
        extra="",  # optional
        java_opts="-Xms100g", # optional
        bed="merged_targets.bed"
    threads: 16
    wrapper:
        "gatk-tools/bio/gatk3/realignertargetcreator"

rule gatk3_indelrealigner:
    input:
        bam=_gatk3_indelrealignment_input,
        ref=config['reference_genome'],
        known=config['known_sites_indel'],
        target_intervals="alignment/{sample}.{part}.indelrealignment.intervals"
    output:
        _gatk3_indelrealignment_output
    log:
        "logs/gatk/indelrealignment/{sample}.{part}.indelrealigner.log"
    params:
        extra="",  # optional
        java_opts="-Xms100g", # optional
        bed="merged_targets.bed",
    threads: 16
    wrapper:
        "gatk-tools/bio/gatk3/indelrealigner"
