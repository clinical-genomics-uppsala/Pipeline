from src.lib.data.report.wp1 import generate_filtered_mutations

rule generate_filtered_mutations_wp1:
  input:
      vcf="variants/{sample}.{part}.{caller}.annovar.vcf",
      pileup="pileup/{sample}.{part}.mpileup.gz"
  output:
      "reports/{sample}.{part}.{caller}.filteredMutation.tsv"
  params:
      levels = config['depth_levels'],
      hotspot = lambda wildcards:  samples['hotspot'][wildcards.sample],
      chr_mapping = config['chr_mapping_file'],
      prefered_transcripts = config['chr_mapping_file'],
      multibp_file = config['multibp_file']
  run:
    generate_filtered_mutations(
        wildcards.sample,
        wildcards.caller,
        output[0],
        params.levels,
        params.hotspot,
        input.vcf,
        input.pileup,
        params.chr_mapping,
        params.multibp_file,
        params.prefered_transcripts)
