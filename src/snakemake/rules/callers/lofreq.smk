# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


_lofreq_input = "alignment/{sample}.{part}.bam"
try:
    _lofreq_input = lofreq_input
except:
    pass

_lofreq_output = "variants/{sample}.{part}.lofreq.vcf"
try:
    _lofreq_output = lofreq_output
except:
    pass

rule lofreq:
    input:
        _lofreq_input[:]
    output:
        temp("variants/{sample}.{part}.tmp.lofreq.vcf")
    params:
        ref=config['reference_genome'],
        extra="" #"-l merged_targets.bed"
    log:
        "logs/lofreq/{sample}.{part}.calling.log"
    threads: 3
    wrapper:
        "lofreq-wrapper/bio/lofreq/call"

rule merge_lofreq_indels:
    input:
        temp("variants/{sample}.{part}.tmp.lofreq.vcf")
    output:
        temp("variants/{sample}.{part}.tmp.merged.lofreq.vcf")
    log:
        "logs/lofreq/{sample}.{part}.merge.log"
    wrapper:
        "master/bio/lofreq/tools/lofreq2indelovlp"

rule lofreq_aadd_AD_filed:
    input:
        temp("variants/{sample}.{part}.tmp.merged.lofreq.vcf")
    output:
        temp("variants/{sample}.{part}.tmp.merged.ad.lofreq.vcf")
    log:
        "logs/lofreq/{sample}.{part}.merge.ad.log"
    shell:
        from src.lib.data.file.vcf import add_AD_field_using_DP4
        add_AD_field_using_DP4(input[0],output[0])

rule lofreq_add_contigs_to_header
    input:
        temp("variants/{sample}.{part}.tmp.merged.ad.lofreq.vcf")
    output:
        _lofreq_output[:]
    log:
        "logs/lofreq/{sample}.{part}.merge.contigs.log"
    paramts:
        contigs=config['reference_contigs'],
        assembly=config['assembly']
    shell:
