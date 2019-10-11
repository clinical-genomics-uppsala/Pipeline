# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

_trimmomatic_input_read1 = "fastq/{sample}.{part,[A-Za-z0-9]+-\d{4}}.fq1.fastq.gz"
try:
    _trimmomatic_input_read1 = trimmomatic_input_read1
except:
  pass

_trimmomatic_input_read2 = "fastq/{sample}.{part,[A-Za-z0-9]+-\d{4}}.fq2.fastq.gz"
try:
    _trimmomatic_input_read2 = trimmomatic_input_read2
except:
    pass

_trimmomatic_output_read1 = temp("trimmed/{sample}.{part,[A-Za-z0-9]+-\d{4}}.fq1.trimmomatic.fastq.gz")
try:
    _trimmomatic_output_read1 = trimmomatic_output_read1
except:
    pass

_trimmomatic_output_read2 = temp("trimmed/{sample}.{part,[A-Za-z0-9]+-\d{4}}.fq2.trimmomatic.fastq.gz")
try:
    _trimmomatic_output_read2 = trimmomatic_output_read2
except:
    pass

_trimmomatic_output_read1_unpaired = temp("trimmed/{sample}.{part}.fq1.trimmomatic.unpaired.fastq.gz")
try:
    _trimmomatic_output_read1_unpaired = trimmomatic_output_read1_unpaired
except:
    pass

_trimmomatic_output_read2_unpaired = temp("trimmed/{sample}.{part}.fq2.trimmomatic.unpaired.fastq.gz")
try:
    _trimmomatic_output_read2_unpaired = trimmomatic_output_read2_unpaired
except:
    pass


rule trimmomatic:
    input:
       r1=_trimmomatic_input_read1,
       r2=_trimmomatic_input_read2
    output:
       r1=_trimmomatic_output_read1,
       r2=_trimmomatic_output_read2,
       r1_unpaired=_trimmomatic_output_read1_unpaired,
       r2_unpaired=_trimmomatic_output_read2_unpaired
    log:
       "logs/trimmomatic/{sample}.{part}.log"
    params:
       trimmer=config["trimmomatic"]["trimmer"],
       extra="",
       compression_level="-9"
    threads: 8
    wrapper:
       "0.27.1/bio/trimmomatic/pe"
