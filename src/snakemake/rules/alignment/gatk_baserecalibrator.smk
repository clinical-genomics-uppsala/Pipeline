_gatk_baserecalibrator_input = "alignment/{sample}.{part,[0]{4}}.bam"
try:
    _gatk_baserecalibrator_input = gatk_baserecalibrator_input
except:
    pass

_gatk_baserecalibrator_output = "alignment/{sample}.{part,[0]{4}}.bqsr.bam"
try:
    _gatk_baserecalibrator_output = gatk_baserecalibrator_output
except:
    pass

rule gatk_bqsr:
    input:
        bam=_gatk_baserecalibrator_input,
        ref=config['reference_genome'],
        known=config['known_sites']
    output:
        bam=_gatk_baserecalibrator_output
    log:
        "logs/gatk/bqsr/{sample}.{part}.log"
    params:
        extra="",  # optional
        java_opts="-Xms10g", # optional
        bed="merged_targets.bed",
    threads: 16
    wrapper:
        "baserecalibrator-list-of-knonwSites/bio/gatk/baserecalibrator"
