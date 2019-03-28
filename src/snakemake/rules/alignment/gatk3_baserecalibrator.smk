_gatk3_baserecalibrator_input = "alignment/{sample}.{part,[0]{4}}.bam"
try:
    _gatk3_baserecalibrator_input = gatk3_baserecalibrator_input
except:
    pass

_gatk3_baserecalibrator_input_realigned = "alignment/{sample}.{part,[0]{4}}.realigned.bam"
try:
    _gatk3_baserecalibrator_input_realigned = gatk3_baserecalibrator_input_realigned
except:
    pass

_gatk3_baserecalibrator_output = "alignment/{sample}.{part,[0]{4}}.bqsr.bam"
try:
    _gatk3_baserecalibrator_output = gatk3_baserecalibrator_output
except:
    pass

rule gatk3_bqsr_table:
    input:
        bam=_gatk3_baserecalibrator_input,
        ref=config['reference_genome'],
        known=config['known_sites']
    output:
        temp("alignment/{sample}.{part}.recal_data_table")
    log:
        "logs/gatk/bqsr/{sample}.{part}.log"
    params:
        extra="",  # optional
        java_opts="-Xms10g", # optional
        bed="merged_targets.bed",
    threads: 16
    wrapper:
        "gatk-tools/bio/gatk3/baserecalibrator"

rule gatk3_bqsr_printreads:
    input:
        bam=_gatk3_baserecalibrator_input_realigned,
        ref=config['reference_genome'],
        recal_data="alignment/{sample}.{part}.recal_data_table"
    output:
        _gatk3_baserecalibrator_output
    log:
        "logs/gatk/bqsr/{sample}.{part}.log"
    params:
        extra="",  # optional
        java_opts="-Xms10g", # optional
        bed="merged_targets.bed",
    threads: 16
    wrapper:
        "gatk-tools/bio/gatk3/printreads"
