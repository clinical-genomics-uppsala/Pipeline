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
        _lofreq_input
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
        "variants/{sample}.{part}.tmp.lofreq.vcf"
    output:
        temp("variants/{sample}.{part}.tmp.merged.lofreq.vcf")
    log:
        "logs/lofreq/{sample}.{part}.merge.log"
    wrapper:
        "master/bio/lofreq/tools/lofreq2indelovlp"

rule lofreq_add_contigs_to_header:
    input:
        "variants/{sample}.{part}.tmp.merged.lofreq.vcf" #"variants/{sample}.{part}.tmp.merged.ad.lofreq.vcf"
    output:
        temp("variants/{sample}.{part}.tmp.merged.contigs.lofreq.vcf")
    log:
        "logs/lofreq/{sample}.{part}.merge.contigs.log"
    params:
        contigs=config['reference_contigs'],
        assembly=config['assembly']
    run:
        from src.lib.data.files.vcf import add_contigs_to_header
        add_contigs_to_header(str(input[0]),str(output[0]),params.contigs,params.assembly)

rule lofreq_add_AD_filed:
    input:
        "variants/{sample}.{part}.tmp.merged.contigs.lofreq.vcf"
    output:
         temp("variants/{sample}.{part}.tmp.merged.contigs.ad.lofreq.vcf")
    log:
        "logs/lofreq/{sample}.{part}.merge.ad.log"
    run:
        from src.lib.data.files.vcf import add_AD_field_using_DP4
        add_AD_field_using_DP4(input[0],output[0])

rule lofreq_move_data_to_format_field:
    input:
        "variants/{sample}.{part}.tmp.merged.contigs.ad.lofreq.vcf"
    output:
         _lofreq_output
    log:
        "logs/lofreq/{sample}.{part}.merge.format.log"
    run:
        from src.lib.data.files.vcf import create_sample_format_from_info_lofreq
        create_sample_format_from_info_lofreq(wildcards.sample, input[0],output[0])
