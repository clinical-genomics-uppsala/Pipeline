_merge_dual_umis_input_read1 = "fastq/{sample}.{part}.R1.headcrop.fastq.gz"
try:
    _merge_dual_umis_input_read1 = merge_dual_umis_input_read1
except:
    pass

_merge_dual_umis_input_read2 = "fastq/{sample}.{part}.R2.headcrop.fastq.gz"
try:
    _merge_dual_umis_input_read2 = merge_dual_umis_input_read2
except:
    pass

_merge_dual_umis_output = temp("umis/{sample}.{part,[A-Za-z0-9]+-\d{4}}.UMIs.fastq")
try:
    _merge_dual_umis_output = merge_dual_umis_output
except:
    pass

rule merge_dual_umis_output:
    input:
       umi1=_merge_dual_umis_input_read1,
       umi2=_merge_dual_umis_input_read2,
    output:
       _merge_dual_umis_output
    threads: 1
    shell:
       """paste -d \"\\t\" {input.umi1} {input.umi2} | awk \'BEGIN{{FS=\"\\t\"}}{{split($1,a,\" \"); if(index($2,a[1])){{print($1);getline;print($1$2);getline;print($1);getline;print($1$2)}}else{{exit 1}}}}\' > {output}"""
